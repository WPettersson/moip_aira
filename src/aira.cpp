#include <atomic>
#include <condition_variable>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <mutex>
#include <thread>

#include <ilcplex/cplex.h>

#include <boost/program_options.hpp>

#include "problem.h"
#include "env.h"
#include "symgroup.h"
#include "result.h"
#include "solutions.h"

#define ERR_CPLEX -1

//#define DEBUG
//#define FINETIMING

#ifdef DEBUG
std::mutex debug_mutex;
#endif

enum Status { RUNNING, DONE };

/* The threads need to synchronise the limits on objectives 3 through n */
struct Locking_Vars {
  std::mutex status_mutex;
  std::mutex ready_mutex;
  std::atomic<Status> thread_status;
  std::condition_variable cv;
  std::condition_variable ready_cv;
} ;

int num_threads;
int cplex_threads;

/* perms[i] is the permutation (as indexed by sym_group.cpp) that will be
  * used by thread i. */
int * perms = nullptr;

double midpoint = 0;
bool split;
/* Number of IPs we've solved */
std::atomic<int> ipcount;


/* Solve CLMOIP and return status */
int solve(Env & e, Problem & p, int * result, double * rhs, int perm_id);

int read_lp_problem(Env& e, Problem& p, bool store_objectives);
int read_mop_problem(Env& e, Problem& p, bool store_objectives);

/* Optimise!
 * first_result is the result of the optimisation with no constraints on
 * objective values.
 * p is the problem (class)
 * thread_id is the index into S_n of the order in which we optimise each
 * objective.
 **/
void optimise(int thread_id, const char * fname, Solutions &all,
    std::mutex & solutionMutex, Locking_Vars *lv,
    std::atomic<double> *shared_limits, std::atomic<double> *global_limits,
    std::list<int*> * my_feasibles, std::list<int*> * partner_feasibles);

bool problems_equal(const Result * a, const Result * b, int objcnt);

namespace po = boost::program_options;

int main (int argc, char *argv[])
{

  int status = 0; /* Operation status */

  std::string permString;
  Env e;

  std::string pFilename, outputFilename;

  /* Timing */
  clock_t starttime, endtime;
  double cpu_time_used, elapsedtime, startelapsed;
  int solcount;

  po::variables_map v;
  po::options_description opt("Options for aira");
  opt.add_options()
    ("help,h", "Show this help.")
    ("lp,p",
      po::value<std::string>(&pFilename),
     "The LP file to solve. Required.")
    ("output,o",
      po::value<std::string>(&outputFilename),
     "The output file. Optional.")
    ("split,s",
     po::bool_switch(&split),
     "Split the range of the first objective into one strip per thread\n"
     "Optional, default to False.")
    ("threads,t",
      po::value<int>(&num_threads)->default_value(1),
     "Number of threads to use internally. Optional, default to 1.")
    ("cplex_threads,c",
        po::value<int>(&cplex_threads)->default_value(1),
     "Number of threads to allocate to CPLEX.\n"
     "Note that each internal thread calls CPLEX, so the total number of "
     "threads used is threads*cplex_threads.\n"
     "Optional, defaults to 1.")
    ("perms", po::value<std::string>(&permString),
     "The permutations (as indexed by sym_group.cpp) to be used by each "
     "the threads. These must be entered as a comma separated list\n"
     "  e.g. 0,1,4,5,8,9\n"
     "Note that to ensure correctness, for even i we must have\n"
     "  permutation[i] being even\n"
     "and\n"
     "  permutation[i+1] = permutation[i] + 1\n")
  ;

  po::store(po::parse_command_line(argc, argv, opt), v);
  po::notify(v);

  if (v.count("help")) {
    // usage();
    std::cout << opt << std::endl;
    return(1);
  }


  if (v.count("lp") == 0) {
    // usage();
    std::cerr << "Error: You must pass in a problem file. All other parameters are optional." << std::endl;
    std::cout << opt << std::endl;
    return(1);
  }

  Problem p(pFilename.c_str(), cplex_threads);



  std::ofstream outFile;

  if (v.count("output")) {
    outFile.open(outputFilename);
  } else {
    /* Output file: <filename>.out */
    size_t last = pFilename.find_last_of(".");
    outFile.open(pFilename.substr(0, last).append(".out"));
  }



  /* Initialize the CPLEX environment */
  e.env = CPXopenCPLEX (&status);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, p.cplex_threads);

  if (e.env == NULL) {
    std::cerr << "Could not open CPLEX environment." << std::endl;
    return -ERR_CPLEX;
  }

  status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
  if (status) {
    std::cerr << "Failure to turn off screen indicator, error." << std::endl ;
      return -ERR_CPLEX;
  }


  {
    int file_status;
    if (p.filetype == LP) {
      file_status = read_lp_problem(e, p, true /* store_objective*/);
    } else if (p.filetype == MOP) {
      file_status = read_mop_problem(e, p, true /* store_objective*/);
    } else {
      std::cerr << "Unknown filetype" << std::endl;
      return -1;
    }
    if (file_status != 0) {
      std::cerr << "Error reading file" << std::endl;
      return -1;
    }
  }
  /* Free up memory as all we needed was a count of the number of objectives. */
  if ( e.lp != NULL ) {
     status = CPXfreeprob (e.env, &e.lp);
     if ( status ) {
       std::cerr << "CPXfreeprob failed." << std::endl;
     }
  }
  if ( e.env != NULL ) {
     status = CPXcloseCPLEX (&e.env);

     if ( status ) {
       std::cerr << "Could not close CPLEX environment." << std::endl;
     }
  }

  // Assign permutations to threads
  if (v.count("perms") == 1) {
    perms = new int[num_threads];
    int index = 0;
    std::string::iterator it = permString.begin();
    while (it != permString.end()) {
      if (index >= num_threads) {
        std::cerr << "Error in permutation string - too many permutations "
            "for " << num_threads << " threads." << std::endl;
        return -1;
      }
      // Move past commas
      while (*it == ',')
        it++;
      // Check to see if we ended on a comma
      if (it == permString.end())
        break;
      std::string token(it, permString.end());
      // Find next comma
      size_t next_comma = token.find(',');
      // If no more commas, this is the last entry
      if (std::string::npos == next_comma) {
        it = permString.end();
        perms[index] = std::atoi(token.c_str());
      } else {
        it += next_comma;
        token.resize(next_comma);
        perms[index] = std::atoi(token.c_str());
      }

      // Make sure the permutation is valid.
      if (perms[index] >= S[p.objcnt].size()) {
        std::cerr << "Invalid permutation : " << perms[index] << " is too big" <<
          "(" << S[p.objcnt].size() << ")." << std::endl;
        return -1;
      }
      if ((index % 2 == 0)) {
        if ((perms[index] % 2) != 0) {
          std::cerr << "Invalid permutation string at " << perms[index] <<
            " - must be even." << std::endl;
          return -1;
        }
      } else {
        if (perms[index] != perms[index-1]+1) {
          std::cerr << "Invalid permutation string at " << perms[index] <<
            " - must equal one more than previous permutation (" <<
            perms[index-1] << ")." << std::endl;
          return -1;
        }
      }
      index++;
    }
    if (index != num_threads) {
      std::cerr << "Error in permutation string - not enough permutations "
          "for " << num_threads << " threads." << std::endl;
      return -1;
    }
  }

  outFile << std::endl << "Using improved algorithm" << std::endl;

  /* Start the timer */
  starttime = clock();
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  startelapsed = start.tv_sec + start.tv_nsec/1e9;


  if (num_threads > S[p.objcnt].size())
    num_threads = S[p.objcnt].size();

  std::list<std::thread> threads;
  Solutions all(p.objcnt);
  std::mutex solutionMutex;

  std::atomic<double> *global_limits = new std::atomic<double>[p.objcnt];

  double lim = (p.objsen == MIN) ? CPX_INFBOUND : -CPX_INFBOUND;
  for (int i = 0; i < p.objcnt; ++i) {
    global_limits[i] = lim;
  }

  std::list<Locking_Vars*> locking_var_list;
  for (int t = 0; t < num_threads; t += 2) {
    std::atomic<double> *shared_limits = new std::atomic<double>[p.objcnt];
    Locking_Vars *lv = new Locking_Vars;
    lv->thread_status = RUNNING;
    locking_var_list.push_back(lv);
    for (int i = 0; i < p.objcnt; ++i) {
      shared_limits[i] = lim;
    }
    std::list<int *> *t1_solns = new std::list<int *>;
    std::list<int *> *t2_solns = new std::list<int *>;
    threads.emplace_back(optimise,
        t, pFilename.c_str(), std::ref(all), std::ref(solutionMutex), lv,
        shared_limits, global_limits, t1_solns, t2_solns);
    // Only launch second thread in pair if we have an even number of threads
    if ((num_threads % 2) == 0) {
      threads.emplace_back(optimise,
          t+1, pFilename.c_str(), std::ref(all), std::ref(solutionMutex), lv,
          shared_limits, global_limits, t2_solns, t1_solns);
    }
  }
  for (auto& thread: threads)
    thread.join();

  for (auto lv: locking_var_list)
    delete lv;

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  clock_gettime(CLOCK_MONOTONIC, &start);
  elapsedtime = (start.tv_sec + start.tv_nsec/1e9 - startelapsed);

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
        outFile << soln->result[i] << "\t";
      }
      outFile << std::endl;
      printed.push_back(soln);
      solcount++;
    }
  }

  constexpr int width = 8;
  constexpr int precision = 3;
  outFile << std::endl << "---" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << cpu_time_used << " CPU seconds" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << elapsedtime << " elapsed seconds" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << ipcount << " IPs solved" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << solcount << " Solutions found" << std::endl;


  outFile.close();

  return (status);
}

/* Solve CLMOIP and return solution status */
int solve(Env & e, Problem & p, int * result, double * rhs, int perm_id) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];
  const int *perm = S[p.objcnt][perm_id];

  //for(int i = 0; i < p.objcnt; ++i)
  //  srhs[i] = rhs[perm[i]];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);

  // TODO Permutation applies here.
  for (int j_preimage = 0; j_preimage < p.objcnt; j_preimage++) {
    int j = perm[j_preimage];
    status = CPXchgobj(e.env, e.lp, cur_numcols, p.objind[j], p.objcoef[j]);
    if (status) {
      std::cerr << "Failed to set objective." << std::endl;
    }

    status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, srhs);
    if (status) {
      std::cerr << "Failed to change constraint srhs" << std::endl;
    }

    /* solve for current objective*/
    status = CPXmipopt (e.env, e.lp);
    if (status) {
      std::cerr << "Failed to optimize LP." << std::endl;
    }

    // This is shared across threads, but it's an atomic integer (so read/writes
    // are atomic) and thus we don't need a lock/mutex
    ipcount++;

    solnstat = CPXgetstat (e.env, e.lp);
    if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
       break;
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
      std::cerr << "Failed to obtain objective value." << std::endl;
      exit(0);
    }
    if ( objval > 1/p.mip_tolerance ) {
      while (objval > 1/p.mip_tolerance) {
        p.mip_tolerance /= 10;
      }
      CPXsetdblparam(e.env, CPXPARAM_MIP_Tolerances_MIPGap, p.mip_tolerance);
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  delete[] srhs;

  return solnstat;
}


void optimise(int thread_id, const char * pFilename, Solutions & all,
    std::mutex &solutionMutex, Locking_Vars *lv,
    std::atomic<double> *shared_limits, std::atomic<double> *global_limits,
    std::list<int *> * my_feasibles, std::list<int *> * partner_feasibles) {
  Env e;
  Problem p(pFilename, cplex_threads);
  int status;
#ifdef FINETIMING
  double cplex_time = 0;
  double wait_time = 0;
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  double total_time = start.tv_sec + start.tv_nsec/1e9;
#endif
  /* Initialize the CPLEX environment */
  e.env = CPXopenCPLEX (&status);
  if (status != 0) {
    std::cerr << "Could not open CPLEX environment." << std::endl;
  }

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, p.cplex_threads);

  if (e.env == NULL) {
    std::cerr << "Could not open CPLEX environment." << std::endl;
  }

  status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
  if (status) {
    std::cerr << "Failure to turn off screen indicator." << std::endl;
  }


  if (p.filetype == LP) {
    read_lp_problem(e, p, true /* store_objective*/);
  } else if (p.filetype == MOP) {
    read_mop_problem(e, p, true /* store_objective*/);
  }
  Solutions s(p.objcnt);

  int perm_id = thread_id;
  if (nullptr != perms) {
    perm_id = perms[thread_id];
  }
  const int* perm = S[p.objcnt][perm_id];

  int infcnt;
  bool inflast;
  bool infeasible;
  int * max, * min,  * result, *resultStore;
  double * rhs;

  result = resultStore = new int[p.objcnt];
  rhs = new double[p.objcnt];
  for(int i = 0; i < p.objcnt; ++i) {
    rhs[i] = p.rhs[i];
  }
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  double starttime = (start.tv_sec + start.tv_nsec/1e9);
#endif
  int solnstat = solve(e, p, result, rhs, perm_id);
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
#ifdef DEBUG
  debug_mutex.lock();
  std::cout << "Thread " << thread_id << " with constraints ∞* found ";
  for(int i = 0; i < p.objcnt; ++i) {
    std::cout << result[i] << ",";
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  // Share the midpoint if splitting
  if (split && num_threads > 1) {
    {
      std::unique_lock<std::mutex> lk(lv->status_mutex);
      if (midpoint == 0) {
        // First in, wait.
        midpoint += ((double)result[0])/2;
        lv->cv.wait(lk);
      } else {
        // Second in, update and keep going
        midpoint += ((double)result[0])/2;
      }
    }
    lv->cv.notify_all();
  }

  /* Need to add a result to the list here*/
  s.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
  if (solnstat != CPXMIP_INFEASIBLE) {
    int *objectives = new int[p.objcnt];
    for (int i = 0; i < p.objcnt; ++i) {
      objectives[i] = result[i];
    }
    partner_feasibles->push_back(objectives);
  }
  min = new int[p.objcnt];
  max = new int[p.objcnt];

#ifdef DEBUG
    debug_mutex.lock();
    std::cout << "Thread " << thread_id << " using permutation ";
    for (int i = 0; i < p.objcnt; ++i) {
      std::cout << perm[i] <<  ", ";
    }
    std::cout << std::endl;
    debug_mutex.unlock();
#endif
  std::atomic<double> &my_limit = shared_limits[perm[0]];
  std::atomic<double> &partner_limit = shared_limits[perm[1]];

  if ((solnstat != CPXMIP_INFEASIBLE) && (p.objcnt > 1)) {
    partner_limit = result[perm[1]];
  }
  for (int j = 0; j < p.objcnt; j++) {
    min[j] = max[j] = result[j];
  }
  for (int objective_counter = 1; objective_counter < p.objcnt; objective_counter++) {
    int objective = perm[objective_counter];
    int depth_level = 1; /* Track current "recursion" depth */
    int depth = perm[depth_level]; /* Track depth objective */
    int onwalk = false; /* Are we on the move? */
    infcnt = 0; /* Infeasible count*/
    inflast = false; /* Last iteration infeasible?*/

    /* Set all constraints back to infinity*/
    for (int j = 0; j < p.objcnt; j++) {
      if (split && (thread_id == 0) && (j == 0))
        rhs[j] = midpoint;
      else
        rhs[j] = global_limits[j];
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
      const Result *relaxation;

#ifdef DEBUG_SOLUTION_SEARCH
      debug_mutex.lock();
      std::cout << "Thread " << thread_id;
      debug_mutex.unlock();
#endif
      relaxation = s.find(rhs, p.objsen);
      relaxed = (relaxation != nullptr);
      if (relaxed) {
        infeasible = relaxation->infeasible;
        result = relaxation->result;
      } else {
        /* Solve in the absence of a relaxation*/
        result = resultStore;
#ifdef FINETIMING
        clock_gettime(CLOCK_MONOTONIC, &start);
        double starttime = (start.tv_sec + start.tv_nsec/1e9);
#endif
        solnstat = solve(e, p, result, rhs, perm_id);
#ifdef FINETIMING
        clock_gettime(CLOCK_MONOTONIC, &start);
        cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
        infeasible = ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD));
        /* Store result */
        s.insert(rhs, result, infeasible);
      }
      if (split) {
        if ((thread_id == 1) && !infeasible) {
          // check if we cross midpoint if we are thread 1
          if (p.objsen == MIN) {
            if (result[0] < midpoint) {
              infeasible = true;
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << thread_id << " reached midpoint";
              std::cout << " " << midpoint << std::endl;
              debug_mutex.unlock();
#endif
            }
          } else {
            if (result[0] > midpoint) {
              infeasible = true;
            }
          }
        }
      } else {
        // Not splitting.
        /* We want to keep the actual objective vector, and share it with our
        * partner as a relaxation. */
        if (!infeasible) {
          int *objectives = new int[p.objcnt];
          for (int i = 0; i < p.objcnt; ++i) {
            objectives[i] = result[i];
          }
          partner_feasibles->push_back(objectives);
        }
        if (!infeasible && (p.objcnt > 1) && (infcnt == 0)) {
          partner_limit = rhs[perm[1]];
        //  rhs[perm[0]] = my_limit;
        }
      }
#ifdef DEBUG
      debug_mutex.lock();
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
      if (relaxed)
        std::cout << "relaxation ";
      if (infeasible) {
        std::cout << "infeasible." << std::endl;
      } else {
        for(int i = 0; i < p.objcnt; ++i) {
          std::cout << result[i] << ",";
        }
        std::cout << std::endl;
      }
      if (!infeasible && (num_threads > 1) && (infcnt == 0)) {
        if (p.objsen == MIN) {
          if (result[perm[0]] >= my_limit) {
            std::cout << "Thread " << thread_id << " result found by partner, bailing." << std::endl;
          }
        } else {
          if (result[perm[0]] <= my_limit) {
            std::cout << "Thread " << thread_id << " result found by partner, bailing." << std::endl;
          }
        }
      }
      debug_mutex.unlock();
#endif
      if (!infeasible && (num_threads > 1) && (infcnt == 0)) {
        if (p.objsen == MIN) {
          if (result[perm[0]] >= my_limit) {
            // Pretend infeasible to backtrack properly
            infeasible = true;
          }
        } else {
          if (result[perm[0]] <= my_limit) {
            // Pretend infeasible to backtrack properly
            infeasible = true;
          }
        }
        // Duplicate code as we are marking this result infeasible
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
      if (infeasible) {
        infcnt++;
        inflast = true;
      } else {
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

      if (infeasible && (infcnt == 1) && (num_threads > 1)) { // Wait/share results of 2-objective problem
#ifdef DEBUG
        debug_mutex.lock();
        std::cout << "Thread " << thread_id <<  " done" << std::endl;
        debug_mutex.unlock();
#endif
        bool wait = false;
        {
          std::unique_lock<std::mutex> lk(lv->status_mutex);
          if (lv->thread_status == RUNNING) {
            lv->thread_status = DONE; // This thread was first in.
            for(int i = 2; i < p.objcnt; ++i) {
              if (p.objsen == MIN) {
                shared_limits[perm[i]] = max[perm[i]];
              } else {
                shared_limits[perm[i]] = min[perm[i]];
              }
            }
            wait = true;
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            double start_wait = start.tv_sec + start.tv_nsec/1e9;
#endif
            lv->cv.wait(lk); // Wait for partner to update limits.
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            wait_time += (start.tv_sec + start.tv_nsec/1e9) - start_wait;
#endif
            for(int i = 2; i < p.objcnt; ++i) {
              if (p.objsen == MIN) {
                if (shared_limits[perm[i]] > max[perm[i]]) {
                  max[perm[i]] = shared_limits[perm[i]];
                }
              } else {
                if (shared_limits[perm[i]] < min[perm[i]]) {
                  min[perm[i]] = shared_limits[perm[i]];
                }
              }
            }
            /* From each feasible result that our partner calculated, we can
            * "make up" an IP problem which would have that result as it's
            * solution. Since we know that the list from our partner, combined
            * with our own results, is the complete list of solutions, this is
            * doable. */
            if (!my_feasibles->empty()) {
              double *lp = new double[p.objcnt];
              for( int i = 2; i < p.objcnt; ++i) {
                lp[perm[i]] = rhs[perm[i]];
              }
              lp[perm[0]] = global_limits[perm[0]];
              int *res = my_feasibles->front();
              my_feasibles->pop_front();
              if (p.objsen == MIN) {
                lp[perm[1]] = res[perm[1]] - 1;
              } else {
                lp[perm[1]] = res[perm[1]] + 1;
              }
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << thread_id << " adding relaxation ";
              for (int j = 0; j < p.objcnt; ++j) {
                if (lp[j] > 1e19)
                  std::cout << "+∞,";
                else if (lp[j] < -1e19)
                  std::cout << "-∞,";
                else
                  std::cout << lp[j] << ",";
              }
              std::cout << " with result infeasible";
              std::cout << std::endl;
              debug_mutex.unlock();
#endif
              s.insert(lp, res, true);
              int *last = res;
              while (! my_feasibles->empty()) {
                res = my_feasibles->front();
                my_feasibles->pop_front();
                if (p.objsen == MIN) {
                  lp[perm[1]] = res[perm[1]] - 1;
                } else {
                  lp[perm[1]] = res[perm[1]] + 1;
                }
#ifdef DEBUG
                debug_mutex.lock();
                std::cout << "Thread " << thread_id << " adding relaxation ";
                for (int j = 0; j < p.objcnt; ++j) {
                  if (lp[j] > 1e19)
                    std::cout << "+∞,";
                  else if (lp[j] < -1e19)
                    std::cout << "-∞,";
                  else
                    std::cout << lp[j] << ",";
                }
                std::cout << " with result ";
                for (int j = 0; j < p.objcnt; ++j) {
                  if (last[j] > 1e19)
                    std::cout << "+∞,";
                  else if (last[j] < -1e19)
                    std::cout << "-∞,";
                  else
                    std::cout << last[j] << ",";
                }
                std::cout << std::endl;
                debug_mutex.unlock();
#endif
                s.insert(lp, last, false);
                delete[] last;
                last = res;
              }
              delete[] last;
              delete[] lp;
            }
          } else if (lv->thread_status == DONE) {
            // This thread can't be first in, so must be last.
            for(int i = 2; i < p.objcnt; ++i) {
              if (p.objsen == MIN) {
                if (shared_limits[perm[i]] < max[perm[i]]) {
                  shared_limits[perm[i]] = max[perm[i]];
                } else {
                  max[perm[i]] = shared_limits[perm[i]];
                }
              } else {
                if (shared_limits[perm[i]] > min[perm[i]]) {
                  shared_limits[perm[i]] = min[perm[i]];
                } else {
                  min[perm[i]] = shared_limits[perm[i]];
                }
              }
            }
            /* From each feasible result that our partner calculated, we can
            * "make up" an IP problem which would have that result as it's
            * solution. Since we know that the list from our partner, combined
            * with our own results, is the complete list of solutions, this is
            * doable. */
            if (! my_feasibles->empty()) {
              double *lp = new double[p.objcnt];
              for( int i = 2; i < p.objcnt; ++i) {
                lp[perm[i]] = rhs[perm[i]];
              }
              lp[perm[0]] = global_limits[perm[0]];
              int *res = my_feasibles->front();
              my_feasibles->pop_front();
              if (p.objsen == MIN) {
                lp[perm[1]] = res[perm[1]] - 1;
              } else {
                lp[perm[1]] = res[perm[1]] + 1;
              }
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << thread_id << " adding relaxation ";
              for (int j = 0; j < p.objcnt; ++j) {
                if (lp[j] > 1e19)
                  std::cout << "+∞,";
                else if (lp[j] < -1e19)
                  std::cout << "-∞,";
                else
                  std::cout << lp[j] << ",";
              }
              std::cout << " with result infeasible";
              std::cout << std::endl;
              debug_mutex.unlock();
#endif
              s.insert(lp, res, true);
              int *last = res;
              while (! my_feasibles->empty()) {
                res = my_feasibles->front();
                my_feasibles->pop_front();
                if (p.objsen == MIN) {
                  lp[perm[1]] = res[perm[1]] - 1;
                } else {
                  lp[perm[1]] = res[perm[1]] + 1;
                }
#ifdef DEBUG
                debug_mutex.lock();
                std::cout << "Thread " << thread_id << " adding relaxation ";
                for (int j = 0; j < p.objcnt; ++j) {
                  if (lp[j] > 1e19)
                    std::cout << "+∞,";
                  else if (lp[j] < -1e19)
                    std::cout << "-∞,";
                  else
                    std::cout << lp[j] << ",";
                }
                std::cout << " with result ";
                for (int j = 0; j < p.objcnt; ++j) {
                  if (last[j] > 1e19)
                    std::cout << "+∞,";
                  else if (last[j] < -1e19)
                    std::cout << "-∞,";
                  else
                    std::cout << last[j] << ",";
                }
                std::cout << std::endl;
                debug_mutex.unlock();
#endif
                s.insert(lp, last, false);
                delete[] last;
                last = res;
              }
              delete[] last;
              delete[] lp;
            }

            /* We do modify partner_limit here, which normally we wouldn't do.
            * However, we need to know that partner_limit is reset before this
            * thread continues, and as there is no more waits, this is the
            * easiest way of achieving this. Note that we reset limits before
            * notifying the other thread, so the other thread will not start the
            * next run until limits are reset too. */

            // Note the typecasting to avoid the deleted copy constructor.
            my_limit = (double)global_limits[perm[0]];
            partner_limit = (double)global_limits[perm[1]];
            lv->thread_status = RUNNING;
          }
        }
        if(wait) {
          // Make sure that the second thread has finished its tasks.
          // Note that we have to wait for this as otherwise we might notify
          // the second thread before it is ready to be notified (causing 2nd
          // thread to get stuck forever).
          {
            std::unique_lock<std::mutex> ready_lk(lv->ready_mutex);
          }
          // Tell second thread to run.
          lv->ready_cv.notify_all();
        } else {
          std::unique_lock<std::mutex> ready_lk(lv->ready_mutex);
          // Notify first thread to keep going
          lv->cv.notify_all();
#ifdef FINETIMING
          clock_gettime(CLOCK_MONOTONIC, &start);
          double start_wait = start.tv_sec + start.tv_nsec/1e9;
#endif
          // Wait for first thread to be ready
          lv->ready_cv.wait(ready_lk);
#ifdef FINETIMING
          clock_gettime(CLOCK_MONOTONIC, &start);
          wait_time += (start.tv_sec + start.tv_nsec/1e9) - start_wait;
#endif
        }
      }

      if ((p.objcnt > 2) && (infcnt == objective_counter) && (infcnt == p.objcnt - 2)) {
        if (p.objsen == MIN) {
          global_limits[perm[p.objcnt-1]] = max[perm[p.objcnt-1]]-1;
        } else {
          global_limits[perm[p.objcnt-1]] = min[perm[p.objcnt-1]]+1;
        }
      }
      if (infcnt == objective_counter-1) {
        if ((p.objcnt > 2) && (objective_counter == p.objcnt - 1)) {
          if (p.objsen == MIN) {
            global_limits[objective] = max[objective]-1;
          } else {
            global_limits[objective] = min[objective]+1;
          }
        }
        /* Set all contraints back to infinity */
        for (int j = 0; j < p.objcnt; j++) {
          rhs[j] = global_limits[j];
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
        if (p.objsen == MIN) {
          rhs[objective] = max[objective]-1;
          max[objective] = (int) -CPX_INFBOUND;
        } else {
          rhs[objective] = min[objective]+1;
          min[objective] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth_level = 1;
        depth = perm[depth_level];
        onwalk = false;
      } else if (inflast && infcnt != objective_counter) {
        rhs[depth] = CPX_INFBOUND;
        depth_level++;
        depth = perm[depth_level];
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = true;
      } else if (!onwalk && infcnt != 1) {
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      } else if (onwalk && infcnt != 1)  {
        depth_level = 1;
        depth = perm[depth_level];
        if (p.objsen == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = false;
      }
    }
  }
  solutionMutex.lock();
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  total_time = start.tv_sec + start.tv_nsec/1e9 - total_time;
  std::cout << "Thread " << thread_id << " used " << cplex_time << "s in cplex";
  std::cout << ", waited for " << wait_time << "s";
  std::cout << " and " << total_time << "s overall." << std::endl;
#endif
  all.merge(s);
  solutionMutex.unlock();
  delete[] resultStore;
  delete[] rhs;
  delete[] min;
  delete[] max;
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

int read_lp_problem(Env& e, Problem& p, bool store_objectives) {
  int status;
  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, p.filename());

  if (e.lp == NULL) {
    std::cerr << "Failed to create LP." << std::endl;
    return -ERR_CPLEX;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, p.filename(), NULL);
  if (status) {
    std::cerr << "Failed to read and copy the problem data." << std::endl;
    return -ERR_CPLEX;
  }

  // If we aren't storing objectives, we are now done.
  if (!store_objectives) {
    return 0;
  }

  /* Get last rhs and determine the number of objectives.*/
  int cur_numcols = CPXgetnumcols(e.env, e.lp);
  int cur_numrows = CPXgetnumrows(e.env, e.lp);
  int cur_numnz = CPXgetnumnz(e.env, e.lp);

  p.rhs = new double[cur_numrows];

  /* Get RHS of last objective - this tells us the objective count */
  status = CPXgetrhs (e.env, e.lp, p.rhs, cur_numrows-1, cur_numrows-1);

  if (status) {
    std::cerr << "Failed to get RHS." << std::endl;
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
    for (int i = 0; i < cur_numcols; ++i)
      p.objcoef[j][i] = 0;
  }

  /* Parse out the objectives working backwards from the last constraint */
  int * rmatbeg = new int[cur_numrows];
  int * rmatind = new int[cur_numnz];
  double * rmatval = new double[cur_numnz];
  int rmatspace = cur_numnz;

  int nzcnt, surplus;
  status = CPXgetrows (e.env, e.lp, &nzcnt, rmatbeg, rmatind, rmatval,
                      rmatspace, &surplus, cur_numrows-p.objcnt, cur_numrows-1);

  if (status) {
    std::cerr << "Couldn't get rows." << std::endl;
    return -ERR_CPLEX;
  }

  for (int j = 0; j < p.objcnt; j++) {
    int to;
    int from = rmatbeg[j];
    if (j == p.objcnt-1) {
      to = nzcnt-1;
    }
    else {
      to = rmatbeg[(j+1)] - 1;
    }
    for (int k = from; k <= to; k++) {
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
    std::cerr << "Failed to change constraint sense" << std::endl;
    return -ERR_CPLEX;
  }

  /* Set rhs of objective constraints */
  status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, p.rhs);
  if (status) {
    std::cerr << "Failed to change constraint rhs" << std::endl;
    return -ERR_CPLEX;
  }
  return 0;
}


int read_mop_problem(Env& e, Problem& p, bool store_objectives) {
  int status;
  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, p.filename());

  if (e.lp == NULL) {
    std::cerr << "Failed to create problem." << std::endl;
    return -ERR_CPLEX;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, p.filename(), NULL);
  if (status) {
    std::cerr << "Failed to read and copy the problem data." << std::endl;
    return -ERR_CPLEX;
  }


  size_t cur_numcols = CPXgetnumcols(e.env, e.lp);
  size_t cur_numrows = CPXgetnumrows(e.env, e.lp);


  char** colNames = new char*[cur_numcols];
  size_t sz = cur_numcols * 1024;
  char* store = new char[sz];
  int surplus;
  status = CPXgetcolname(e.env, e.lp, colNames, store, sz, &surplus, 0, cur_numcols-1);
  if (status) {
    std::cerr << "Error retrieving column names." << std::endl;
    return -ERR_CPLEX;
  }


  // Iterate through file until we get to ROWS
  // count number of objectives
  // store names
  // break when we get COLUMNS
  // read a b c
  // if b == objective_name[i] p.objcoeff[i][ name_to_index[a]] = int(c)
  std::fstream mop_file(p.filename(), std::ios::in);
  std::string line;
  // Find "ROWS" line
  while (std::getline( mop_file, line)) {
    if ("ROWS" == line) {
      break;
    }
  }
  // Read objective names
  std::vector<std::string> objNames;
  while (std::getline( mop_file, line)) {
    // Objective functions have N as the second character
    if ('N' != line[1]) {
      break;
    }
    std::stringstream ss(line);
    std::string n, name;
    ss >> n; // Skip "N"
    ss >> name;
    objNames.push_back(name);
  }
  if (store_objectives) {
    p.objcnt = static_cast<int>(objNames.size());
  }
  int new_nzcnt = p.objcnt * cur_numcols;
  double * newRowMatVal = new double[new_nzcnt] {0};
  int * newRowMatInd = new int[new_nzcnt];

  // We will store each coefficient, in order, so set matrix up
  for(int i = 0; i < p.objcnt; ++i) {
    for(size_t j = 0; j < cur_numcols; ++j) {
      newRowMatInd[i*cur_numcols + j] = j;
    }
  }
  int * newRowMatBeg = new int[objNames.size()];
  for(int i = 0; i < p.objcnt; ++i) {
    newRowMatBeg[i] = i*cur_numcols;
  }

  if (store_objectives) {
    /* Create a pair of multidimensional arrays to store the objective
    * coefficients and their indices indices first */
    p.objind = new int*[p.objcnt];
    for(int j = 0; j < p.objcnt; j++) {
      p.objind[j] = new int[cur_numcols];
      for(size_t i = 0; i < cur_numcols; i++){
        p.objind[j][i] = i;
      }
    }
    /* Now coefficients */
    p.objcoef = new double*[p.objcnt];
    for(int j = 0; j < p.objcnt; j++) {
      p.objcoef[j] = new double[cur_numcols];
      for (size_t i = 0; i < cur_numcols; ++i) {
        p.objcoef[j][i] = 0;
      }
    }
  }
  // Find columns
  while (std::getline( mop_file, line)) {
    if ("COLUMNS" == line) {
      break;
    }
  }

  /* Parse out the objectives working backwards from the last constraint */
  while (std::getline( mop_file, line)) {
    std::stringstream ss(line);
    std::string col, obj;
    signed int val;
    ss >> col;
    ss >> obj;
    ss >> val;
    size_t colInd;
    int objInd;
    for(colInd = 0; colInd < cur_numcols; ++colInd) {
      if ( col.compare(colNames[colInd]) == 0) {
        break;
      }
    }
    // Skip if we don't recognise because why not?
    if (colInd == cur_numcols)
      continue;

    for(objInd = 0; objInd < p.objcnt; ++objInd) {
      if (objNames[objInd] == obj) {
        break;
      }
    }
    if (objInd == p.objcnt) {
      // Just an inequality, which we've already read.
      continue;
    }
    if (store_objectives) {
      p.objcoef[objInd][colInd] = val;
    }
    newRowMatVal[objInd * cur_numcols + colInd] = val;
  }

  // Now we need to set up the RHS.
  p.rhs = new double[p.objcnt];

  if (status) {
    std::cerr << "Failed to get RHS." << std::endl;
    return -ERR_CPLEX;
  }
  if (store_objectives) {
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
  }

  delete[] colNames;
  colNames = new char*[objNames.size()];
  for(size_t i = 0; i < objNames.size(); ++i) {
    colNames[i] = new char[objNames[i].size()+1];
    for(size_t j = 0; j < objNames[i].size(); j++) {
      colNames[i][j] = objNames[i][j];
    }
    colNames[i][objNames.size()] = '\0';
  }
  CPXaddrows(e.env, e.lp, 0 /*ccnt*/, objNames.size(), new_nzcnt, p.rhs, p.consense, newRowMatBeg, newRowMatInd, newRowMatVal, colNames, NULL /*rowname*/);

  p.conind = new int[p.objcnt];
  /* Specify index of objective constraints */
  for (int j = 0; j < p.objcnt; j++) {
    p.conind[j] = cur_numrows+j;
  }

  /* Set sense of objective constraints */
  status = CPXchgsense (e.env, e.lp, p.objcnt, p.conind, p.consense);
  if (status) {
    std::cerr << "Failed to change constraint sense" << std::endl;
    return -ERR_CPLEX;
  }

  delete[] newRowMatVal;
  delete[] newRowMatInd;
  delete[] newRowMatBeg;
  return 0;
}
