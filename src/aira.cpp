#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include <libgen.h>
#include <time.h>

#include <atomic>
#include <iostream>
#include <list>
#include <mutex>
#include <thread>

#include "problem.h"
#include "env.h"
#include "symgroup.h"
#include "result.h"
#include "solutions.h"

#define ERR_CPLEX -1

#define DEBUG

/* Number of IPs we've solved */
std::atomic<int> ipcount;

/* The filename of the problem */
char *lpfn;

/* Solve CLMOIP and return status */
int solve(Env & e, Problem & p, int * result, double * rhs, int thread_id, std::atomic<int> & ipcount);

/* Optimise!
 * first_result is the result of the optimisation with no constraints on
 * objective values.
 * p is the problem (class)
 * thread_id is the index into S_n of the order in which we optimise each
 * objective.
 **/
void optimise(int thread_id, Problem & p, Solutions &all, std::mutex & solutionMutex);

bool problems_equal(const Result * a, const Result * b, int objcnt);

int main (int argc, char *argv[])
{

  int status = 0; /* Operation status */

  /* LP File Processing  */
  int from, to; /* First row, last row */
  int cur_numrows, cur_numcols, cur_numnz; /* # rows in current problem, # cols, # nonzero values*/
  int rmatspace; /* Length of rmatind & rmatval */
  int nzcnt, surplus; /* Nonzero row value count, array length check */
  int *rmatbeg, *rmatind; /* Row indices, column indices of rmatval */
  double *rmatval; /* Constraint RHS, relaxation ip, nozero row values */

  Problem p;
  Env e;

  /* LP and log file */
  FILE *outfp;
  char *outfn, *tmpstr;


  /* Timing */
  clock_t starttime, endtime;
  time_t startelapsed, endelapsed;
  double cpu_time_used, elapsedtime;
  int solcount;

  int ipcount_nonatomic;

  if (argc < 3) {
    printf("\nUsage: aira <numthreads> <filename.lp>\n\n");
    exit(1);
  }

  int num_threads = atoi(argv[1]);

  /* Last arg is the lp file */
  lpfn = argv[argc-1];

  /* Output file: <filename>.out */
  tmpstr = (char *) malloc(strlen(lpfn) * sizeof(char));
  outfn = (char *) malloc(strlen(lpfn) * sizeof(char));
  tmpstr = strrchr(lpfn, '.');
  strncpy(outfn, lpfn, tmpstr-lpfn);
  outfn[tmpstr-lpfn] = '\0';
  strcat(outfn, ".out");
  outfp = fopen(outfn, "w+");

  /* Initialize the CPLEX environment */
  e.env = CPXopenCPLEX (&status);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, 1);

  if (e.env == NULL) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (e.env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
    return -ERR_CPLEX;
  }

 status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
 if (status) {
   fprintf (stderr,
       "Failure to turn off screen indicator, error %d.\n", status);
    return -ERR_CPLEX;
 }

  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, lpfn);

  if (e.lp == NULL) {
    fprintf (stderr, "Failed to create LP.\n");
    return -ERR_CPLEX;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, lpfn, NULL);
  if (status) {
    fprintf (stderr, "Failed to read and copy the problem data.\n");
    return -ERR_CPLEX;
  }

  /* Get last rhs and determine the number of objectives.*/
  cur_numcols = CPXgetnumcols(e.env, e.lp);
  cur_numrows = CPXgetnumrows(e.env, e.lp);
  cur_numnz = CPXgetnumnz(e.env, e.lp);

  p.rhs = new double[cur_numrows];

  /* Get RHS of last objective - this tells us the objective count */
  status = CPXgetrhs (e.env, e.lp, p.rhs, cur_numrows-1, cur_numrows-1);

  if (status) {
    fprintf (stderr, "Failed to get RHS.\n");
    return -ERR_CPLEX;
  }

  p.objcnt = static_cast<int>(p.rhs[0]);

  /* Create a pair of multidimensional arrays to store the objective
   * coefficients and their indices indices first */
  p.objind = new int*[p.objcnt];
  for(int j = 0; j < p.objcnt; j++) {
    p.objind[j] = new int[cur_numcols];
    for(int i = 0; i < cur_numcols; i++){
      p.objind[j][i] = i;
    }
  }
  /* Now coefficients */
  p.objcoef = new double*[p.objcnt];
  for(int j = 0; j < p.objcnt; j++) {
    p.objcoef[j] = new double[cur_numcols];
    memset(p.objcoef[j], 0, cur_numcols * sizeof(double));
  }

  /* Parse out the objectives working backwards from the last constraint */
  rmatbeg = new int[cur_numrows];
  rmatind = new int[cur_numnz];
  rmatval = new double[cur_numnz];
  rmatspace = cur_numnz;

  status = CPXgetrows (e.env, e.lp, &nzcnt, rmatbeg, rmatind, rmatval,
                      rmatspace, &surplus, cur_numrows-p.objcnt, cur_numrows-1);

  if (status) {
    fprintf (stderr, "Couldn't get rows.\n");
    return -ERR_CPLEX;
  }

  for (int j = 0; j < p.objcnt; j++) {
    from = rmatbeg[j];
    if (j == p.objcnt-1) {
      to = nzcnt-1;
    }
    else {
      to = rmatbeg[(j+1)] - 1;
    }
    for (int k = from; k <= to; k++){
      p.objcoef[j][rmatind[k]] = rmatval[k];
    }
  }
  delete[] rmatbeg;
  delete[] rmatind;
  delete[] rmatval;

  /* Setup problem for solving */

  /* Resize rhs to fit only objective function constraints */
  delete[] p.rhs;
  p.rhs = new double[p.objcnt];

  /* Get objective sense */
  int cpx_sense = CPXgetobjsen(e.env, e.lp);
  p.objsen = (cpx_sense == CPX_MIN ? MIN : MAX);
  /* Set objective constraint sense and RHS */
  p.consense = new char[p.objcnt];
  for (int j = 0; j < p.objcnt; j++) {
    if (p.objsen == MIN) {
      p.consense[j] = 'L'; /* Set sense to <= */
      p.rhs[j] = CPX_INFBOUND;
    }
    else {
      p.consense[j] = 'G'; /* Set sense to >= */
      p.rhs[j] = -CPX_INFBOUND;
    }
  }

  /* Specify index of objective constraints */
  p.conind = new int[p.objcnt];
  for (int j = 0, k = p.objcnt; j < p.objcnt; j++, k--) {
    p.conind[j] = cur_numrows-k;
  }

  /* Set sense of objective constraints */
  status = CPXchgsense (e.env, e.lp, p.objcnt, p.conind, p.consense);
  if (status) {
    fprintf (stderr, "Failed to change constraint sense\n");
    return -ERR_CPLEX;
  }

  /* Set rhs of objective constraints */
  status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, p.rhs);
  if (status) {
    fprintf (stderr, "Failed to change constraint rhs\n");
    return -ERR_CPLEX;
  }

  fprintf(outfp, "\nUsing improved algorithm\n");

  /* Start the timer */
  starttime = clock();
  startelapsed = time(NULL);


  if (num_threads > S[p.objcnt].size())
    num_threads = S[p.objcnt].size();

  std::list<std::thread> threads;
  Solutions all(p.objcnt);
  std::mutex solutionMutex;
  for(int t = 0; t < num_threads; ++t) {
    threads.emplace_back(optimise, t, std::ref(p), std::ref(all), std::ref(solutionMutex));
  }
  for (auto& thread: threads)
    thread.join();
  //optimise(0 /* thread_id */, p /* Problem */, result /* first_result */);

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  endelapsed = time (NULL);
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  elapsedtime=(double) endelapsed - startelapsed;

  //list = g_slist_sort(list, (GCompareFunc)icmp);
  //list = g_slist_reverse(list);
  all.sort();
  solcount = 0;

  /* List of things we've printed */
  std::list<const Result *> printed;
  for (const Result* soln: all) {
    if (soln->infeasible)
      continue;
    bool printMe = true;
    for (auto& done: printed) {
      if (problems_equal(done, soln, p.objcnt)) {
        printMe = false;
        break;
      }
    }
    if (printMe) {
      for (int i = 0; i < p.objcnt; i++) {
        fprintf(outfp,"%d \t", soln->result[i]);
      }
      fprintf(outfp,"\n");
      printed.push_back(soln);
      solcount++;
    }
  }

  fprintf(outfp,"\n---\n%8.4f CPU seconds\n", cpu_time_used);
  fprintf(outfp,"%8.4f elapsed seconds\n", elapsedtime);
  ipcount_nonatomic = ipcount;
  fprintf(outfp,"%8d IPs solved\n", ipcount_nonatomic);
  fprintf(outfp,"%8d Solutions found\n", solcount);

  /* Free up memory as necessary. */
  if ( e.lp != NULL ) {
     status = CPXfreeprob (e.env, &e.lp);
     if ( status ) {
        fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
     }
  }
  if ( e.env != NULL ) {
     status = CPXcloseCPLEX (&e.env);

     if ( status ) {
        char  errmsg[CPXMESSAGEBUFSIZE];
        fprintf (stderr, "Could not close CPLEX environment.\n");
        CPXgeterrorstring (e.env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
     }
  }


  if (outfp != NULL) fclose(outfp);

  return (status);
}

/* Solve CLMOIP and return solution status */
int solve(Env & e, Problem & p, int * result, double * rhs, int thread_id, std::atomic<int> & ipcount) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];
  const int *perm = S[p.objcnt][thread_id];

  //for(int i = 0; i < p.objcnt; ++i)
  //  srhs[i] = rhs[perm[i]];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);


  // TODO Permutation applies here.
  for (int j_preimage = 0; j_preimage < p.objcnt; j_preimage++) {
    int j = perm[j_preimage];
    status = CPXchgobj(e.env, e.lp, cur_numcols, p.objind[j], p.objcoef[j]);
    if (status) {
       fprintf (stderr, "Failed to set objective.\n");
    }

    status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, srhs);
    if (status) {
      fprintf (stderr, "Failed to change constraint srhs\n");
    }

    /* solve for current objective*/
    status = CPXmipopt (e.env, e.lp);
    if (status) {
       fprintf (stderr, "Failed to optimize LP.\n");
    }

    // This is shared across threads, but it's an atomic integer (so read/writes
    // are atomic) and thus we don't need a lock/mutex
    ipcount++;

    solnstat = CPXgetstat (e.env, e.lp);
    if (solnstat == CPXMIP_INFEASIBLE) {
       break;
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
       fprintf (stderr, "Failed to obtain objective value.\n");
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  return solnstat;
}


void optimise(int thread_id, Problem & p, Solutions & all, std::mutex &solutionMutex) {
  Env e;
  Solutions s(p.objcnt);
  int status;
  /* Initialize the CPLEX environment */
  e.env = CPXopenCPLEX (&status);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, 1);

  if (e.env == NULL) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (e.env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
  }

 status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
 if (status) {
   fprintf (stderr,
       "Failure to turn off screen indicator, error %d.\n", status);
 }

  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, lpfn);


  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, lpfn, NULL);

  /* Set sense of objective constraints */
  status = CPXchgsense (e.env, e.lp, p.objcnt, p.conind, p.consense);

  /* Set rhs of objective constraints */
  status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, p.rhs);

  // Have some permutation that maps preimage_i to i
  // with a different permutation for each thread?
  // with a different lead objective for each thread?
  int infcnt;
  bool inflast;
  bool infeasible;
  int * max, * min,  * result;
  double * rhs, * ipr;

  result = new int[p.objcnt];
  int solnstat = solve(e, p, result, p.rhs, thread_id, ipcount);

  /* Need to add a result to the list here*/
  s.insert(p.rhs, result, solnstat == CPXMIP_INFEASIBLE);
#ifdef DEBUG
    std::cout << "Thread " << thread_id << " with constraints ∞,∞,∞ found ";
        for(int i = 0; i < p.objcnt; ++i) {
          std::cout << result[i] << ",";
        }
        std::cout << std::endl;
#endif
  rhs = new double[p.objcnt];
  min = new int[p.objcnt];
  max = new int[p.objcnt];
  ipr = new double[p.objcnt];


  for (int j = 0; j < p.objcnt; j++) {
    rhs[j] = p.rhs[j];
    min[j] = max[j] = result[j];
  }
  const int* perm = S[p.objcnt][thread_id];
  for (int objective_counter = 1; objective_counter < p.objcnt; objective_counter++) {
    int objective = perm[objective_counter];
#ifdef DEBUG
    std::cout << "Thread " << thread_id << " optimising obj # " << objective
       << std::endl;
#endif
    int depth_level = 1; /* Track current "recursion" depth */
    int depth = perm[depth_level]; /* Track depth objective */
    int onwalk = false; /* Are we on the move? */
    infcnt = 0; /* Infeasible count*/
    inflast = false; /* Last iteration infeasible?*/

    /* Set all contraints back to infinity*/
    for (int j = 0; j < p.objcnt; j++) {
      if (p.objsen == MIN) {
        rhs[j] = CPX_INFBOUND;
      }
      else {
        rhs[j] = -CPX_INFBOUND;
      }
    }

    /* Set rhs of current depth */
    if (p.objsen == MIN) {
      rhs[objective] = max[objective]-1;
    }
    else {
      rhs[objective] = min[objective]+1;
    }
    max[objective] = (int) -CPX_INFBOUND;
    min[objective] = (int) CPX_INFBOUND;
    while (infcnt < objective_counter) {
      bool relaxed;
      int solnstat;
      /* Look for possible relaxations to the current problem*/
#ifdef DEBUG
      std::cout << "Thread " << thread_id << " with constraints ";
      for(int i = 0; i < p.objcnt; ++i) {
        if (rhs[i] > 1e19)
          std::cout << "∞";
        else if (rhs[i] < -1e19)
          std::cout << "-∞";
        else
          std::cout << rhs[i];
        std::cout << ",";
      }
      std::cout << " found ";
#endif
      const Result *relaxation;

      relaxation = s.find(rhs, p.objsen);
      relaxed = (relaxation != nullptr);
      if (relaxed) {
        infeasible = relaxation->infeasible;
        ipr = relaxation->ip;
        result = relaxation->result;
#ifdef DEBUG
        std::cout << "relaxation ";
#endif
      } else {
        /* Solve in the absence of a relaxation*/
        solnstat = solve(e, p, result, rhs, thread_id, ipcount);
        infeasible = (solnstat == CPXMIP_INFEASIBLE);
      }
      if (infeasible) {
#ifdef DEBUG
        std::cout << "infeasible." << std::endl;
#endif
        infcnt++;
        inflast = true;
      }
      else {
#ifdef DEBUG
        for(int i = 0; i < p.objcnt; ++i) {
          std::cout << result[i] << ",";
        }
        std::cout << std::endl;
#endif
        infcnt = 0;
        inflast = false;
        /* Update maxima */
        for (int j = 0; j < p.objcnt; j++) {
          if (result[j] > max[j]) {
            max[j] = result[j];
          }
        }
        /* Update minima */
        for (int j = 0; j < p.objcnt; j++) {
          if (result[j] < min[j]) {
            min[j] = result[j];
          }
        }
      }
      /* Store result */
      if (!relaxed) {
        s.insert(rhs, result, infeasible);
      }

      if (infcnt == objective_counter-1) {
        /* Set all contraints back to infinity */
        for (int j = 0; j < p.objcnt; j++) {
          if (p.objsen == MIN) {
            rhs[j] = CPX_INFBOUND;
          }
          else {
            rhs[j] = -CPX_INFBOUND;
          }
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
        if (p.objsen == MIN) {
          rhs[objective] = max[objective]-1;
          max[objective] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[objective] = min[objective]+1;
          min[objective] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth_level = 1;
        depth = perm[depth_level];
        onwalk = false;
      }
      else if (inflast && infcnt != objective_counter) {
        rhs[depth] = CPX_INFBOUND;
        depth_level++;
        depth = perm[depth_level];
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = true;
      }
      else if (!onwalk && infcnt != 1) {
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      }
      else if (onwalk && infcnt != 1)  {
        depth_level = 1;
        depth = perm[depth_level];
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = false;
      }
    }
  }
  solutionMutex.lock();
  all.merge(s);
  solutionMutex.unlock();
}


bool problems_equal(const Result * a, const Result * b, int objcnt) {
  // Are they both feasible, or both infeasible
  if (a->infeasible != b->infeasible) {
    return false;
  }
  // If one is infeasible, both are
  if (a->infeasible) {
    return true;
  }
  // Have to check actual result/objective values
  for (int i = 0; i < objcnt; i++) {
    if (a->result[i] != b->result[i]) {
      return false;
    }
  }
  return true;
}
