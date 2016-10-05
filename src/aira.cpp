#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include <libgen.h>
#include <time.h>

#include <list>

/* List node. */
class Problem {
  public:
    double *ip;
    int *result;
    int infeasible;
};

/* List of possible relaxations */
std::list<Problem*> list;

/* List of things we've printed */
std::list<Problem *> printed;

/* CPLEX environment and problem */
CPXENVptr env = NULL;
CPXLPptr  lp = NULL;

/* Objective count */
int objcnt;

/* Solve CLMOIP and return status */
int solve(int objcnt, int **objind, double **objcoef, int *conind, double *rhs, int *result, int *ipcount);

/* Add new solution/problem to list of possible relaxations */
void add_to_list(double *ip, int *result, bool infeasible);

/* Search list for relaxations.
 * ip: Current problem
 * ipr: Place holder for the relaxation to the current problem
 */
int getrelaxation(double *ip, double *ipr, int *result, bool *infeasible, int objsen);

/* Compare two solutions (from struct Problem)
 * returns negative value if a < b; zero if a = b; positive value if a > b.
 * */
bool icmp(const Problem * a, const Problem * b);

bool problems_equal(const Problem * a, const Problem * b);

int main (int argc, char *argv[])
{

  int status = 0; /* Operation status */
  int solnstat; /* Solution status (from CPLEX) */

  /* LP File Processing  */
  int from, to; /* First row, last row */
  int cur_numrows, cur_numcols, cur_numnz; /* # rows in current problem, # cols, # nonzero values*/
  int rmatspace, cur_objsen; /* Length of rmatind & rmatval, objective sense of current objective */
  int nzcnt, surplus; /* Nonzero row value count, array length check */
  int *rmatbeg, *rmatind; /* Row indices, column indices of rmatval */
  double *rhs, *ipr, *rmatval; /* Constraint RHS, relaxation ip, nozero row values */
  int **objind; /* Objective indices */
  double **objcoef; /* Objective coefficients */
  int *conind; /* Constraint indices */
  char *consense; /* Constraint senses */

  /* Algorithm */
  int infcnt; /* Infeasibility count */
  bool *infeasible; /* Infeasibility flag for relaxation */
  int inflast, relaxed; /* Last iteration was feasible?, did we find a relaxation? */
  int *result, *max, *min; /* Current result, max result, min result */

  /* LP and log file */
  FILE *outfp;
  char *lpfn, *outfn, *tmpstr;


  /* Timing */
  clock_t starttime, endtime;
  time_t startelapsed, endelapsed;
  double cpu_time_used, elapsedtime;
  int *ipcount;
  int solcount;

  if (argc < 2) {
    printf("\nUsage: aira <filename.lp>\n\n");
    exit(1);
  }

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
  env = CPXopenCPLEX (&status);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(env, CPXPARAM_Threads, 1);

  if (env == NULL) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
    goto TERMINATE;
  }

 status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
 if (status) {
   fprintf (stderr,
       "Failure to turn off screen indicator, error %d.\n", status);
   goto TERMINATE;
 }

  /* Create the problem, using the filename as the problem name */
  lp = CPXcreateprob(env, &status, lpfn);

  if (lp == NULL) {
    fprintf (stderr, "Failed to create LP.\n");
    goto TERMINATE;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(env, lp, lpfn, NULL);
  if (status) {
    fprintf (stderr, "Failed to read and copy the problem data.\n");
    goto TERMINATE;
  }

  /* Get last rhs and determine the number of objectives.*/
  cur_numcols = CPXgetnumcols(env, lp);
  cur_numrows = CPXgetnumrows(env, lp);
  cur_numnz = CPXgetnumnz(env, lp);

  rhs = new double[cur_numrows];
  if (rhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Get RHS of last objective - this tells us the objective count */
  status = CPXgetrhs (env, lp, rhs, cur_numrows-1, cur_numrows-1);

  if (status) {
    fprintf (stderr, "Failed to get RHS.\n");
    goto TERMINATE;
  }

  objcnt = static_cast<int>(rhs[0]);

  /* Create a pair of multidimensional arrays to store the objective
   * coefficients and their indices indices first */
  objind = new int*[objcnt];
  if (objind == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for(int j = 0; j < objcnt; j++) {
    objind[j] = (int *) malloc(cur_numcols * sizeof(int));
    if (objind[j] == NULL) {
      fprintf(stderr, "Error! memory is not available\n");
      goto TERMINATE;
    }
    for(int i = 0; i < cur_numcols; i++){
      objind[j][i] = i;
    }
  }
  /* Now coefficients */
  objcoef = static_cast<double **>(malloc(objcnt * sizeof(double *)));
  if (objcoef == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for(int j = 0; j < objcnt; j++) {
    objcoef[j] = (double *) malloc(cur_numcols * sizeof(double));
    if (objcoef[j] == NULL) {
      fprintf(stderr, "Error! memory is not available\n");
      goto TERMINATE;
    }
    memset(objcoef[j], 0, cur_numcols * sizeof(double));
  }

  /* Parse out the objectives working backwards from the last constraint */
  rmatbeg = (int *) malloc(cur_numrows * sizeof(int));
  rmatind = (int *) malloc(cur_numnz * sizeof(int));
  rmatval = (double *) malloc(cur_numnz * sizeof(double));
  if (rmatbeg == NULL || rmatind == NULL || rmatval == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  rmatspace = cur_numnz;

  status = CPXgetrows (env, lp, &nzcnt, rmatbeg, rmatind, rmatval,
                      rmatspace, &surplus, cur_numrows-objcnt, cur_numrows-1);

  if (status) {
    fprintf (stderr, "Couldn't get rows.\n");
    goto TERMINATE;
  }

  for (int j = 0; j < objcnt; j++) {
    from = rmatbeg[j];
    if (j == objcnt-1) {
      to = nzcnt-1;
    }
    else {
      to = rmatbeg[(j+1)] - 1;
    }
    for (int k = from; k <= to; k++){
      objcoef[j][rmatind[k]] = rmatval[k];
    }
  }

  /* Setup problem for solving */

  /* Resize rhs to fit only objective function constraints */
  delete[] rhs;
  rhs = new double[objcnt];
  if (rhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Get objective sense */
  cur_objsen = CPXgetobjsen(env, lp);

  /* Set objective constraint sense and RHS */
  consense = (char *) malloc(objcnt * sizeof(char));
  if (consense == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for (int j = 0; j < objcnt; j++) {
    if (cur_objsen == CPX_MIN) {
      consense[j] = 'L'; /* Set sense to <= */
      rhs[j] = CPX_INFBOUND;
    }
    else {
      consense[j] = 'G'; /* Set sense to >= */
      rhs[j] = -CPX_INFBOUND;
    }
  }

  /* Specify index of objective constraints */
  conind = (int *) malloc(objcnt * sizeof(int));
  if (conind == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for (int j = 0, k = objcnt; j < objcnt; j++, k--) {
    conind[j] = cur_numrows-k;
  }

  /* Set sense of objective constraints */
  status = CPXchgsense (env, lp, objcnt, conind, consense);
  if (status) {
     fprintf (stderr, "Failed to change constraint sense\n");
     goto TERMINATE;
  }

  /* Set rhs of objective constraints */
  status = CPXchgrhs (env, lp, objcnt, conind, rhs);
  if (status) {
     fprintf (stderr, "Failed to change constraint rhs\n");
     goto TERMINATE;
  }

  /* Allocate an array to store the result */
  result = (int *) malloc(objcnt * sizeof(int));
  if (result == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Allocate an array to track the maximum values */
  max = (int *) malloc(objcnt * sizeof(int));
  if (max == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  /* Allocate an array to track the minimum values */
  min = (int *) malloc(objcnt * sizeof(int));
  if (min == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  ipcount = (int *) malloc(sizeof(int));
  if (ipcount == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  *ipcount = 0;

  fprintf(outfp, "\nUsing improved algorithm\n");
  ipr = (double *) malloc(objcnt * sizeof(double));
  if (ipr == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Start the timer */
  starttime = clock();
  startelapsed = time(NULL);

  solnstat = solve(objcnt, objind, objcoef, conind, rhs, result, ipcount);
  memcpy(max, result, objcnt * sizeof(int));
  memcpy(min, result, objcnt * sizeof(int));

  infeasible = new bool;
  *infeasible = false;

  /* Need to add a result to the list here*/
  add_to_list(rhs, result, *infeasible);

  // Have some permutation that maps preimage_i to i
  // with a different permutation for each thread?
  // with a different lead objective for each thread?
  for (int i = 1; i < objcnt; i++) {
    int depth = 1; /* Track current "recursion" depth */
    int onwalk = false; /* Are we on the move? */
    infcnt = 0; /* Infeasible count*/
    inflast = false; /* Last iteration infeasible?*/

    /* Set all contraints back to infinity*/
    for (int j = 0; j < objcnt; j++) {
      if (cur_objsen == CPX_MIN) {
        rhs[j] = CPX_INFBOUND;
      }
      else {
        rhs[j] = -CPX_INFBOUND;
      }
    }

    /* Set rhs of current depth */
    if (cur_objsen == CPX_MIN) {
      rhs[i] = max[i]-1;
    }
    else {
      rhs[i] = min[i]+1;
    }
    max[i] = (int) -CPX_INFBOUND;
    min[i] = (int) CPX_INFBOUND;

    while (infcnt < i) {
      relaxed = false;

      /* Look for possible relaxations to the current problem*/
      solnstat = getrelaxation(rhs, ipr, result, infeasible, cur_objsen);
      relaxed = solnstat ? true : false;

      if(!relaxed) {
        /* Solve in the absence of a relaxation*/
        solnstat = solve(objcnt, objind, objcoef, conind, rhs, result, ipcount);
      }

      if (solnstat == CPXMIP_INFEASIBLE || *infeasible) {
        infcnt++;
        inflast = true;
      }
      else {
        infcnt = 0;
        inflast = false;
        /* Update maxima */
        for (int j = 1; j < objcnt; j++) {
          if (result[j] > max[j]) {
            max[j] = result[j];
          }
        }
        /* Update minima */
        for (int j = 1; j<objcnt; j++) {
          if (result[j] < min[j]) {
            min[j] = result[j];
          }
        }
      }
      /* Store result */
      if (!relaxed) {
        add_to_list(rhs, result, *infeasible);
      }

      if (infcnt == i-1) {
        /* Set all contraints back to infinity */
        for (int j = 0; j < objcnt; j++) {
          if (cur_objsen == CPX_MIN) {
            rhs[j] = CPX_INFBOUND;
          }
          else {
            rhs[j] = -CPX_INFBOUND;
          }
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
        if (cur_objsen == CPX_MIN) {
          rhs[i] = max[i]-1;
          max[i] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[i] = min[i]+1;
          min[i] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth = 1;
        onwalk = false;
      }
      else if (inflast && infcnt != i) {
        rhs[depth] = CPX_INFBOUND;
        depth++;
        if (cur_objsen == CPX_MIN) {
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
        if (cur_objsen == CPX_MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      }
      else if (onwalk && infcnt != 1)  {
        depth = 1;
        if (cur_objsen == CPX_MIN) {
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

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  endelapsed = time (NULL);
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  elapsedtime=(double) endelapsed - startelapsed;

  //list = g_slist_sort(list, (GCompareFunc)icmp);
  //list = g_slist_reverse(list);
  list.sort(icmp);
  solcount = 0;

  for (auto soln: list) {
    bool printMe = true;
    for (auto done: printed) {
      if (problems_equal(done, soln)) {
        printMe = false;
        break;
      }
    }
    if (printMe) {
      for (int i = 0; i < objcnt; i++) {
        fprintf(outfp,"%d \t", soln->result[i]);
      }
      fprintf(outfp,"\n");
      printed.push_back(soln);
      solcount++;
    }
  }

  fprintf(outfp,"\n---\n%8.4f CPU seconds\n", cpu_time_used);
  fprintf(outfp,"%8.4f elapsed seconds\n", elapsedtime);
  fprintf(outfp,"%8d IPs solved\n", *ipcount);
  fprintf(outfp,"%8d Solutions found\n", solcount);

  TERMINATE:
   /* Free up memory as necessary. */
   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }
   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }
   if (rhs != NULL)
     delete[] rhs;
   if (rmatbeg != NULL) free(rmatbeg);
   if (rmatind != NULL) free(rmatind);
   if (rmatval != NULL) free(rmatval);
   if (objind !=NULL) {
    for(int j = 0; j < objcnt; j++) {
      if (objind[j] != NULL) free(objind[j]);
    }
   }
   if (objind != NULL)
     delete[] objind;
   if (objcoef !=NULL) {
    for(int j = 0; j < objcnt; j++) {
      if (objcoef[j] != NULL) free(objcoef[j]);
    }
   }
   if (objcoef != NULL) free(objcoef);
   if (conind != NULL) free(conind);
   if (consense != NULL) free(consense);
   if (outfp != NULL) fclose(outfp);

   for(auto p: list)
     delete p;

   return (status);
}

/* Solve CLMOIP and return solution status */
int solve(int objcnt, int **objind, double **objcoef, int *conind, double *rhs,
    int *result, int *ipcount) {

  int j;
  int cur_numcols, status, solnstat;
  double objval, *srhs;
  srhs = (double *) malloc(objcnt * sizeof(double));
  if (srhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    return -1;
  }

  memcpy(srhs, rhs, objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(env, lp);

  // TODO Permutation applies here.
  for (int j = 0; j < objcnt;j++) {
    status = CPXchgobj(env, lp, cur_numcols, objind[j], objcoef[j]);
    if (status) {
       fprintf (stderr, "Failed to set objective.\n");
    }

    status = CPXchgrhs (env, lp, objcnt, conind, srhs);
    if (status) {
      fprintf (stderr, "Failed to change constraint srhs\n");
    }

    /* solve for current objective*/
    status = CPXmipopt (env, lp);
    if (status) {
       fprintf (stderr, "Failed to optimize LP.\n");
    }

    (*ipcount)++;

    solnstat = CPXgetstat (env, lp);
    if (solnstat == CPXMIP_INFEASIBLE) {
       break;
    }
    status = CPXgetobjval (env, lp, &objval);
    if ( status ) {
       fprintf (stderr, "Failed to obtain objective value.\n");
    }

    result[j] = round(objval);
    srhs[j] = round(objval);
  }

  return solnstat;
}

/* Add new solution/problem to list of possible relaxations */
void add_to_list(double *ip, int *result, bool infeasible) {
  Problem* p = new Problem;
  p->ip = new double[objcnt];
  memcpy(p->ip, ip, objcnt * sizeof(double));
  if (infeasible) {
    p->infeasible = true;
    p->result = NULL;
  }
  else {
    p->infeasible = false;
    p->result = new int[objcnt];
    memcpy(p->result, result, objcnt * sizeof(int));
  }
  list.push_back(p);
}

/* Search list for relaxations.
 * ip: Current problem
 * ipr: Place holder for the relaxation to the current problem
 */
int getrelaxation(double *ip, double *ipr, int *result, bool *infeasible, int objsen) {

  int t1,t2,t3,i;
  for (auto prob: list) {

    t1 = true;
    t2 = false;
    t3 = true;

    if(objsen == CPX_MIN) {
      for (i=0; i<objcnt; i++){
        /* t1: All values of candidate must be >= than ip */
        if (prob->ip[i] < ip[i]) {
          t1 = false;
          break;
        }
        /* t2: At least one inequality must be strict */
        if (prob->ip[i] > ip[i]) {
          t2 = true;
        }
        /* t3: All values of candidate result must be <= than ip */
        if (!prob->infeasible) {
          if (prob->result[i] > ip[i]) {
            t3 = false;
            break;
          }
        }
      }
    }
    else {
      for (i=0; i<objcnt; i++){
        /* t1: All values of candidate must be <= than ip */
        if (prob->ip[i] > ip[i]) {
          t1 = false;
          break;
        }
        /* t2: At least one inequality must be strict */
        if (prob->ip[i] < ip[i]) {
          t2 = true;
        }
        /* t3: All values of candidate result must be >= than ip */
        if (!prob->infeasible) {
          if (prob->result[i] < ip[i]) {
            t3 = false;
            break;
          }
        }
      }
    }
    /* If all conditions are met copy problem & solution and return */
    if (t1 && t2 && t3) {
      memcpy(ipr, prob->ip, objcnt * sizeof(double));
      *infeasible = prob->infeasible;
      if (!*infeasible) {
        memcpy(result, prob->result, objcnt * sizeof(int));
      }
      return 1;
    }
  }
  return 0;
}

/* Compare two solutions (from struct Problem)
 * returns negative value if a < b; zero if a = b; positive value if a > b.
 * True if a ≤ b, False otherwise.
 * WARNING: This comparison function has been reversed to save a call to
 * list::reverse. */
bool icmp(const Problem * a, const Problem * b) {
  int i;
  for (i=0;i<objcnt;i++) {
    if (a->result[i] < b->result[i]) {
      return false;
    }
    else if (a->result[i] > b->result[i]) {
      return true;
    }
  }
  return false;
}

bool problems_equal(const Problem * a, const Problem * b) {
  for (int i=0;i<objcnt;i++) {
    if (a->result[i] < b->result[i]) {
      return false;
    }
    else if (a->result[i] > b->result[i]) {
      return false;
    }
  }
  return true;

}
