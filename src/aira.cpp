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

#include "cluster.h"
#include "problem.h"
#include "env.h"
#include "lockingvars.h"
#include "symgroup.h"
#include "result.h"
#include "solutions.h"
#include "errors.h"
#include "thread.h"

//#define DEBUG
//#define FINETIMING

#ifdef DEBUG
std::mutex debug_mutex;
#endif

int num_threads;
int cplex_threads;

/**
 * Are we splitting up the range of values an objective can take, such that
 * individual threads only find solutions in their own pre-determined range?
 */
bool split;

/**
 * Track how many individual IPs we solve.
 */
std::atomic<int> ipcount;


/**
 * Solves a CLMOIP and returns the status and result.
 * Variables:
 *
 * \param e An Env object holding an environment suitable for solving IPs
 * \param p The Problem object to be solved
 * \param result An array of ints of size p.objcnt where the result will be
 * stored
 * \param rhs An array of doubles holding the right hand side of the problem
 * \param perm_id The index of the permutation which denotes the hierarchy or
 * ordering of the objectives, from most significant to least.
 */
int solve(Env & e, Problem & p, int * result, double * rhs, Thread * t);
int solve(Env & e, Problem & p, int * result, double * rhs, const int * perm);

/**
 * Optimise!
 * This function runs through the selected algorithm of Pettersson and
 * Ozlen[1]. It is designed to be run multi-threaded, hence using locking
 * mechanisms.
 *
 * \param thread_id An identifier for the thread running this function
 * \pram fname The file name of the file holding the problem description
 * \param all A Solutions object into which all solutions will be placed
 * \param lv A Locking_Vars struct holding all synchronisation particulars for
 * synchronisting two paired threads.
 * \param shared_limits An array of atomic ints, used to share limits between
 * two paired threads.
 * \param global_limits An array of atomic ints, used to share limits amongst
 * all running threads.
 * \param my_feasibles
 * \param partner_feasibles
 * first_result is the result of the optimisation with no constraints on
 * objective values.
 * p is the problem (class)
 * thread_id is the index into S_n of the order in which we optimise each
 * objective.
 **/
template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Thread *t);
//void optimise(int thread_id, const char * fname, Solutions &all, Locking_Vars
//    *lv, std::atomic<double> *shared_limits, std::atomic<double>
//    *global_limits, std::list<int*> * my_feasibles, std::list<int*> *
//    partner_feasibles, double split_start, double split_stop);

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
     "  e.g. 0,1,4,5,8,9\n")
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

  Problem p(pFilename.c_str(), e);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, cplex_threads);

  if (e.env == NULL) {
    std::cerr << "Could not open CPLEX environment." << std::endl;
    return -ERR_CPLEX;
  }

  status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
  if (status) {
    std::cerr << "Failure to turn off screen indicator, error." << std::endl ;
      return -ERR_CPLEX;
  }



  /**
  * perms[i] is the permutation (as indexed by sym_group.cpp) that will be
  * used by thread i. perms remain a nullptr if no permutation is specified.
  */
  int * perms = nullptr;
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

  Solutions all(p.objcnt);

  std::atomic<double> *global_limits = new std::atomic<double>[p.objcnt];

  double lim = (p.objsen == MIN) ? CPX_INFBOUND : -CPX_INFBOUND;
  for (int i = 0; i < p.objcnt; ++i) {
    global_limits[i] = lim;
  }

  std::list<Locking_Vars*> locking_var_list;
  double start_point, stop_point;
  if (split) {
    double * rhs = new double[p.objcnt];
    for (int i=0; i < p.objcnt; ++i) {
      rhs[i] = p.rhs[i];
    }
    int * result = new int[p.objcnt];
    int solnstat = solve(e, p, result, rhs, S[p.objcnt][0]);
    if ( solnstat == CPXMIP_INFEASIBLE) {
      std::cout << "Infeasible" << std::endl;
      return 0;
    }
    all.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
#ifdef DEBUG
    debug_mutex.lock();
    std::cout << "Main thread with constraints ";
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
    for(int i = 0; i < p.objcnt; ++i) {
      std::cout << result[i] << ",";
    }
    std::cout << std::endl;
    debug_mutex.unlock();
#endif
    start_point = result[p.objcnt-1];
    if (p.objsen == MIN) {
      start_point += 1;
    } else {
      start_point -= 1;
    }

    solnstat = solve(e, p, result, rhs, S[p.objcnt][S[p.objcnt].size()-1]);
    if ( solnstat == CPXMIP_INFEASIBLE) {
      std::cout << "Infeasible" << std::endl;
      return 0;
    }
    all.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
#ifdef DEBUG
    debug_mutex.lock();
    std::cout << "Main thread with constraints ";
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
    for(int i = 0; i < p.objcnt; ++i) {
      std::cout << result[i] << ",";
    }
    std::cout << std::endl;
    debug_mutex.unlock();
#endif
    stop_point = result[p.objcnt-1];
    if (p.objsen == MIN) {
      stop_point -= 1;
    } else {
      stop_point += 1;
    }

    delete[] rhs;
    delete[] result;
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

  std::list<Thread*> threads;
  if (split) {
    double split_start = start_point;
    double split_stop;
    double step_size = (stop_point - start_point)/num_threads;
    for (int t = 0; t < num_threads; ++t) {
      split_stop = split_start + step_size;
      threads.push_back(new Thread(t, p.objcnt, split_start, split_stop));
      split_start = split_stop;
    }
  } else {
    // Not splitting.

    // We can share relaxations and limits if there is something to share them with,
    // and if the next thread is using an appropriate permutation
//    if ((t+1 < num_threads) && (!split) &&
//        (
//         (perms == nullptr) || // No perms specified is ok for sharing.
//         ((perms[t] % 2 == 0) && (perms[t] + 1 == perms[t+1]))
//        )
//      ) { // Can share limits+relaxations
//      shared_limits = new std::atomic<double>[p.objcnt];
//      for (int i = 0; i < p.objcnt; ++i) {
//        shared_limits[i] = lim;
//      }
//    }
//    Locking_Vars *lv = new Locking_Vars;
//    lv->num_running_threads = (t+1) == num_threads ? 1 : 2;
//    locking_var_list.push_back(lv);
    int * ordering = new int[p.objcnt];
    int ** share_from = new int*[p.objcnt] {nullptr};
    int ** share_to = new int*[p.objcnt] {nullptr};
    Locking_Vars ** lvs = new Locking_Vars*[p.objcnt] {nullptr};
    Cluster c(num_threads, p.objcnt, p.objsen, false /* spread_threads */, p.objcnt, ordering, share_from,
        share_to, threads, lvs);
    delete[] ordering;
    delete[] share_to;
    delete[] share_from;
  }
  std::list<std::thread> threadList;
  if (p.objsen == MIN) {
    for(Thread* thread: threads) {
      threadList.emplace_back(optimise<MIN>, pFilename.c_str(),
          std::ref(all), thread);
    }
  } else {
    for(Thread* & thread: threads) {
      threadList.emplace_back(optimise<MAX>, pFilename.c_str(),
          std::ref(all), thread);
    }
  }
  for (auto& thread: threadList)
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
int solve(Env & e, Problem & p, int * result, double * rhs, const int * perm) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];

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
      status = CPXmipopt (e.env, e.lp);
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
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  delete[] srhs;

  return solnstat;
}


int solve(Env & e, Problem & p, int * result, double * rhs, Thread * t) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];

  //for(int i = 0; i < p.objcnt; ++i)
  //  srhs[i] = rhs[perm[i]];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);

  // TODO Permutation applies here.
  for (int j_preimage = 0; j_preimage < p.objcnt; j_preimage++) {
    int j = t->perm(j_preimage);
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
      status = CPXmipopt (e.env, e.lp);
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
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  delete[] srhs;

  return solnstat;
}

template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Thread *t) {
//    std::list<int *> * my_feasibles, std::list<int *> * partner_feasibles,
  Env e;
  const bool sharing = (t->share_to != nullptr);
#ifdef DEBUG
  debug_mutex.lock();
  if (split) {
    std::cout << "Thread " << t->id << " has start ";
    std::cout << t->split_start << " and stop " << t->split_stop << std::endl;
  }
  if (sharing) {
    std::cout << "Thread " << t->id << " is sharing" << std::endl;
  }
  debug_mutex.unlock();
#endif
#ifdef FINETIMING
  double cplex_time = 0;
  double wait_time = 0;
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  double total_time = start.tv_sec + start.tv_nsec/1e9;
#endif
  {
    int status; // Only need status here.
    /* Initialize the CPLEX environment */
    e.env = CPXopenCPLEX (&status);
    if (status != 0) {
      std::cerr << "Could not open CPLEX environment." << std::endl;
    }

    /* Set to deterministic parallel mode */
    status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

    /* Set to only one thread */
    CPXsetintparam(e.env, CPXPARAM_Threads, cplex_threads);

    if (e.env == NULL) {
      std::cerr << "Could not open CPLEX environment." << std::endl;
    }

    status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
    if (status) {
      std::cerr << "Failure to turn off screen indicator." << std::endl;
    }
  }

  Problem p(pFilename, e);

  Solutions s(p.objcnt);

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
  if (split) {
    rhs[t->perm(p.objcnt)] = t->split_start;
  }

#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  double starttime = (start.tv_sec + start.tv_nsec/1e9);
#endif
  int solnstat = solve(e, p, result, rhs, t);
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
#ifdef DEBUG
  debug_mutex.lock();
  std::cout << "Thread " << t->id << " with constraints ";
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
  for(int i = 0; i < p.objcnt; ++i) {
    std::cout << result[i] << ",";
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  /* Need to add a result to the list here*/
  s.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
  // Note that if we are splitting, we aren't sharing.
  if (!split) {
    if (sharing && (solnstat != CPXMIP_INFEASIBLE)) {
      int *objectives = new int[p.objcnt];
      for (int i = 0; i < p.objcnt; ++i) {
        objectives[i] = result[i];
      }
    }
  }
  min = new int[p.objcnt];
  max = new int[p.objcnt];

#ifdef DEBUG
  debug_mutex.lock();
  std::cout << "Thread " << t->id << " using permutation ";
  for (int i = 0; i < p.objcnt; ++i) {
    std::cout << t->perm(i) <<  ", ";
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  if (sharing && (solnstat != CPXMIP_INFEASIBLE) && (p.objcnt > 1)) {
    int i = t->perm(1);
    if (t->share_to[i] != nullptr) {
      if (sense == MIN) {
        if (*t->share_to[i] < result[i]) {
          *t->share_to[i] = result[i];
        }
      } else {
        if (*t->share_to[i] > result[i]) {
          *t->share_to[i] = result[i];
        }
      }
    }
  }
  if (solnstat != CPXMIP_INFEASIBLE) {
    for (int j = 0; j < p.objcnt; j++) {
      min[j] = max[j] = result[j];
    }
  }


  /**
   * objective_counter marks which objectives are currently active. Note that
   * the first objective (highest priority) is 0, which we never make active.
   *
   * That is, if the current permutation is 1,2,3,4,5 and objective_counter is
   * 1, then we are trying to optimise objective 2. In other words, objectives
   * 3,4,5 all have their limits set to +- infinity, and these won't change
   * until objective_counter increases.
   *
   * If objective_counter is 3, then we are changing objectives 2,3 and 4.
   */

  /**
   * depth tracks which of the active objectives is the one being currently
   * modified. If objective_counter is 3, then the active objectives are 2,3
   * and 4. If depth is 1, then objective 2 is the one for which we are
   * changing the RHS.
   */

  /**
   * onwalk is true if the current active depth (as indicated by the depth
   * variable) was just increased.
   */
  for (int objective_counter = 1; objective_counter < p.objcnt; objective_counter++) {
    int objective = t->perm(objective_counter);
    int depth_level = 1; /* Track current "recursion" depth */
    int depth = t->perm(depth_level); /* Track depth objective */
    int onwalk = false; /* Are we on the move? */
    infcnt = 0; /* Infeasible count*/
    inflast = false; /* Last iteration infeasible?*/

    /* Set all constraints back to infinity*/
    for (int j_pre = 1; j_pre < p.objcnt; j_pre++) {
      int j = t->perm(j_pre);
      if ((!sharing) || (t->share_from[j] == nullptr)) {
        if (sense == MIN) {
          rhs[j] = CPX_INFBOUND;
        } else {
          rhs[j] = -CPX_INFBOUND;
        }
      } else {
        rhs[j] = *t->share_from[j];
      }
    }
    if (split) {
        rhs[t->perm(p.objcnt)] = t->split_start;
    }
    /* Set rhs of current depth */
    if (sense == MIN) {
      rhs[objective] = max[objective]-1;
    } else {
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
      std::cout << "Thread " << t->id;
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
        solnstat = solve(e, p, result, rhs, t);
#ifdef FINETIMING
        clock_gettime(CLOCK_MONOTONIC, &start);
        cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
        infeasible = ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD));
        /* Store result */
        s.insert(rhs, result, infeasible);
      }
#ifdef DEBUG
      debug_mutex.lock();
      std::cout << "Thread " << t->id << " with constraints ";
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
      if (sharing && !infeasible && (num_threads > 1) && (infcnt == 0)) {
        if (t->share_from[t->perm(0)] != nullptr) {
          if (sense == MIN) {
            if (result[t->perm(0)] >= *t->share_from[t->perm(0)]) {
              std::cout << "Thread " << t->id << " result found by partner, bailing." << std::endl;
            }
          } else {
            if (result[t->perm(0)] <= *t->share_from[t->perm(0)]) {
              std::cout << "Thread " << t->id << " result found by partner, bailing." << std::endl;
            }
          }
        }
      }
      debug_mutex.unlock();
#endif
      if (split) {
        if (!infeasible) {
          // check if we cross midpoint if we are thread 1
          if (sense == MIN) {
            if (result[p.objcnt-1] < t->split_stop) {
              infeasible = true;
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " reached split_stop";
              std::cout << " " << t->split_stop << std::endl;
              debug_mutex.unlock();
#endif
            }
          } else {
            if (result[p.objcnt-1] > t->split_stop) {
              infeasible = true;
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " reached split_stop";
              std::cout << " " << t->split_stop << std::endl;
              debug_mutex.unlock();
#endif
            }
          }
        }
      } else if (sharing) {
        // Not splitting.
        /* We want to keep the actual objective vector, and share it with our
        * partner as a relaxation. */
//        if (!infeasible) {
//          int *objectives = new int[p.objcnt];
//          for (int i = 0; i < p.objcnt; ++i) {
//            objectives[i] = result[i];
//          }
//          partner_feasibles->push_back(objectives);
//        }

        // Note that perm[1] is only shared with one other thread, that only
        // ever reads it, and that this thread always improves the value of
        // this shared limit, so we can use this faster update.
        if (!infeasible && (p.objcnt > 1) && (infcnt == 0)) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " updating bound on " << t->perm(1);
          std::cout << " to " << rhs[t->perm(1)] << std::endl;
          debug_mutex.unlock();
#endif
          if (t->share_to[t->perm(1)] != nullptr)
            *t->share_to[t->perm(1)] = rhs[t->perm(1)];
        //  rhs[perm[0]] = my_limit;
        }
      }

      // TODO Here, check _every_ value of i, not just i=0.
      if (sharing && !infeasible && (t->share_from[t->perm(0)] != nullptr)
          && (infcnt == 0)) {
        if (sense == MIN) {
          if (result[t->perm(0)] >= *t->share_from[t->perm(0)]) {
            // Pretend infeasible to backtrack properly
            infeasible = true;
          }
        } else {
          if (result[t->perm(0)] <= *t->share_from[t->perm(0)]) {
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

      // If we are sharing, and there are other threads, and this result
      // is infeasible, and it's a 2-objective problem
      if (sharing && infeasible) { //(lv->num_running_threads > 1) && (infcnt == 1)) { // Wait/share results of 2-objective problem
        int updated_objective = t->perm(infcnt+1);
#ifdef DEBUG
        debug_mutex.lock();
        std::cout << "Thread " << t->id <<  " done, ";
        std::cout << "infcnt is " << infcnt;
        std::cout << ", objective_counter is " << objective_counter;
        std::cout << ", updated_objective is " << updated_objective;
        if (sense == MIN) {
          std::cout << " and new value is " << max[updated_objective];
        } else {
          std::cout << " and new value is " << min[updated_objective];
        }
        std::cout << std::endl;
        //<< lv->num_running_threads << "still running.";
        debug_mutex.unlock();
#endif
        bool wait = false;
        Locking_Vars *lv = t->locks[updated_objective];
        if ( lv != nullptr ) {
          std::cout << "Found locking var" << std::endl;
          std::unique_lock<std::mutex> lk(lv->status_mutex);
          if ( lv->num_running_threads > 1) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id <<  " done, ";
            std::cout << " not last in." << std::endl;
            debug_mutex.unlock();
#endif
            lv->num_running_threads--; // This thread is no longer running.
            if (t->share_to[updated_objective] != nullptr) {
              if (sense == MIN) {
                if (*t->share_to[updated_objective] > max[updated_objective]) {
                  *t->share_to[updated_objective] = max[updated_objective];
                }
              } else {
                if (*t->share_to[updated_objective] < min[updated_objective]) {
                  *t->share_to[updated_objective] = min[updated_objective];
                }
              }
            }
            wait = true;
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            double start_wait = start.tv_sec + start.tv_nsec/1e9;
#endif
            lv->cv.wait(lk); // Wait for all other threads to update limits.
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            wait_time += (start.tv_sec + start.tv_nsec/1e9) - start_wait;
#endif
            if (t->share_to[updated_objective] != nullptr) {
              if (sense == MIN) {
                if (*t->share_to[updated_objective] > max[updated_objective]) {
                  max[updated_objective] = *t->share_to[updated_objective];
                }
              } else {
                if (*t->share_to[updated_objective] < min[updated_objective]) {
                  min[updated_objective] = *t->share_to[updated_objective];
                }
              }
            }
            lv->num_running_threads++; // This thread is starting up again.
            /* From each feasible result that our partner calculated, we can
            * "make up" an IP problem which would have that result as it's
            * solution. Since we know that the list from our partner, combined
            * with our own results, is the complete list of solutions, this is
            * doable. */
//            if (!my_feasibles->empty()) {
//              double *lp = new double[p.objcnt];
//              for( int i = 2; i < p.objcnt; ++i) {
//                lp[perm[i]] = rhs[perm[i]];
//              }
//              lp[perm[0]] = global_limits[perm[0]];
//              int *res = my_feasibles->front();
//              my_feasibles->pop_front();
//              if (sense == MIN) {
//                lp[perm[1]] = res[perm[1]] - 1;
//              } else {
//                lp[perm[1]] = res[perm[1]] + 1;
//              }
//#ifdef DEBUG
//              debug_mutex.lock();
//              std::cout << "Thread " << t->id << " adding relaxation ";
//              for (int j = 0; j < p.objcnt; ++j) {
//                if (lp[j] > 1e19)
//                  std::cout << "+∞,";
//                else if (lp[j] < -1e19)
//                  std::cout << "-∞,";
//                else
//                  std::cout << lp[j] << ",";
//              }
//              std::cout << " with result infeasible";
//              std::cout << std::endl;
//              debug_mutex.unlock();
//#endif
//              s.insert(lp, res, true);
//              int *last = res;
//              while (! my_feasibles->empty()) {
//                res = my_feasibles->front();
//                my_feasibles->pop_front();
//                if (sense == MIN) {
//                  lp[perm[1]] = res[perm[1]] - 1;
//                } else {
//                  lp[perm[1]] = res[perm[1]] + 1;
//                }
//#ifdef DEBUG
//                debug_mutex.lock();
//                std::cout << "Thread " << t->id << " adding relaxation ";
//                for (int j = 0; j < p.objcnt; ++j) {
//                  if (lp[j] > 1e19)
//                    std::cout << "+∞,";
//                  else if (lp[j] < -1e19)
//                    std::cout << "-∞,";
//                  else
//                    std::cout << lp[j] << ",";
//                }
//                std::cout << " with result ";
//                for (int j = 0; j < p.objcnt; ++j) {
//                  if (last[j] > 1e19)
//                    std::cout << "+∞,";
//                  else if (last[j] < -1e19)
//                    std::cout << "-∞,";
//                  else
//                    std::cout << last[j] << ",";
//                }
//                std::cout << std::endl;
//                debug_mutex.unlock();
//#endif
//                s.insert(lp, last, false);
//                delete[] last;
//                last = res;
//              }
//              delete[] last;
//              delete[] lp;
//            }
          } else { // (lv->num_running_threads == 1) Should be guaranteed?
            // This thread is last.
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id <<  " done, ";
            std::cout << " last in." << std::endl;
            std::cout << "max is " << max[updated_objective] << std::endl;
            if (t->share_to[updated_objective] != nullptr)
              std::cout << "share_to is " << *t->share_to[updated_objective] << std::endl;
            debug_mutex.unlock();
#endif
            if (t->share_to[updated_objective] != nullptr) {
              if (sense == MIN) {
                if (*t->share_to[updated_objective] < max[updated_objective]) {
                  *t->share_to[updated_objective] = max[updated_objective];
                } else {
                  max[updated_objective] = *t->share_to[updated_objective];
                }
              } else {
                if (*t->share_to[updated_objective] > min[updated_objective]) {
                  *t->share_to[updated_objective] = min[updated_objective];
                } else {
                  min[updated_objective] = *t->share_to[updated_objective];
                }
              }
            }
            lv->cv.notify_all();
            /* From each feasible result that our partner calculated, we can
            * "make up" an IP problem which would have that result as it's
            * solution. Since we know that the list from our partner, combined
            * with our own results, is the complete list of solutions, this is
            * doable. */
//            if (! my_feasibles->empty()) {
//              double *lp = new double[p.objcnt];
//              for( int i = 2; i < p.objcnt; ++i) {
//                lp[perm[i]] = rhs[perm[i]];
//              }
//              lp[perm[0]] = global_limits[perm[0]];
//              int *res = my_feasibles->front();
//              my_feasibles->pop_front();
//              if (sense == MIN) {
//                lp[perm[1]] = res[perm[1]] - 1;
//              } else {
//                lp[perm[1]] = res[perm[1]] + 1;
//              }
//#ifdef DEBUG
//              debug_mutex.lock();
//              std::cout << "Thread " << t->id << " adding relaxation ";
//              for (int j = 0; j < p.objcnt; ++j) {
//                if (lp[j] > 1e19)
//                  std::cout << "+∞,";
//                else if (lp[j] < -1e19)
//                  std::cout << "-∞,";
//                else
//                  std::cout << lp[j] << ",";
//              }
//              std::cout << " with result infeasible";
//              std::cout << std::endl;
//              debug_mutex.unlock();
//#endif
//              s.insert(lp, res, true);
//              int *last = res;
//              while (! my_feasibles->empty()) {
//                res = my_feasibles->front();
//                my_feasibles->pop_front();
//                if (sense == MIN) {
//                  lp[perm[1]] = res[perm[1]] - 1;
//                } else {
//                  lp[perm[1]] = res[perm[1]] + 1;
//                }
//#ifdef DEBUG
//                debug_mutex.lock();
//                std::cout << "Thread " << t->id << " adding relaxation ";
//                for (int j = 0; j < p.objcnt; ++j) {
//                  if (lp[j] > 1e19)
//                    std::cout << "+∞,";
//                  else if (lp[j] < -1e19)
//                    std::cout << "-∞,";
//                  else
//                    std::cout << lp[j] << ",";
//                }
//                std::cout << " with result ";
//                for (int j = 0; j < p.objcnt; ++j) {
//                  if (last[j] > 1e19)
//                    std::cout << "+∞,";
//                  else if (last[j] < -1e19)
//                    std::cout << "-∞,";
//                  else
//                    std::cout << last[j] << ",";
//                }
//                std::cout << std::endl;
//                debug_mutex.unlock();
//#endif
//                s.insert(lp, last, false);
//                delete[] last;
//                last = res;
//              }
//              delete[] last;
//              delete[] lp;
//            }

            /* We do modify partner_limit here, which normally we wouldn't do.
            * However, we need to know that partner_limit is reset before this
            * thread continues, and as there is no more waits, this is the
            * easiest way of achieving this. Note that we reset limits before
            * notifying the other thread, so the other thread will not start the
            * next run until limits are reset too. */

            // Note the typecasting to avoid the deleted copy constructor.
//            my_limit = (double)global_limits[perm[0]];
//            partner_limit = (double)global_limits[perm[1]];
          }
        }
//        if (wait) {
//          // Make sure that the second thread has finished its tasks.
//          // Note that we have to wait for this as otherwise we might notify
//          // the second thread before it is ready to be notified (causing 2nd
//          // thread to get stuck forever).
//          {
//            std::unique_lock<std::mutex> ready_lk(lv->ready_mutex);
//          }
//          // Tell second thread to run.
//          lv->ready_cv.notify_all();
//        } else {
//          std::unique_lock<std::mutex> ready_lk(lv->ready_mutex);
//          // Notify first thread to keep going
//          lv->cv.notify_all();
//#ifdef FINETIMING
//          clock_gettime(CLOCK_MONOTONIC, &start);
//          double start_wait = start.tv_sec + start.tv_nsec/1e9;
//#endif
//          // Wait for first thread to be ready
//          lv->ready_cv.wait(ready_lk);
//#ifdef FINETIMING
//          clock_gettime(CLOCK_MONOTONIC, &start);
//          wait_time += (start.tv_sec + start.tv_nsec/1e9) - start_wait;
//#endif
//        }
      }

      if ((p.objcnt > 2) && (infcnt == objective_counter) && (infcnt == p.objcnt - 2)) {
        if (t->share_to[t->perm(p.objcnt-1)] == nullptr)
          continue;
#ifdef DEBUG
        debug_mutex.lock();
        std::cout << "Thread " << t->id << " updating bound on " << t->perm(p.objcnt-1);
        std::cout << " to ";
        if (sense == MIN) {
          std::cout << max[t->perm(p.objcnt-1)];
        } else {
          std::cout << min[t->perm(p.objcnt-1)];
        }
        std::cout << std::endl;
        debug_mutex.unlock();
#endif
        if (sense == MIN) {
          *t->share_to[t->perm(p.objcnt-1)] = max[t->perm(p.objcnt-1)];
        } else {
          *t->share_to[t->perm(p.objcnt-1)] = min[t->perm(p.objcnt-1)];
        }
      }
      if (infcnt == objective_counter-1) {
        if ((p.objcnt > 2) && (objective_counter == p.objcnt - 1)) {
          if (t->share_to[t->perm(p.objcnt-1)] != nullptr) {
            if (sense == MIN) {
              *t->share_to[objective] = max[objective];
            } else {
              *t->share_to[objective] = min[objective];
            }
          }
        }
        /* Set all constraints back to infinity */
        for (int j = 0; j < p.objcnt; j++) {
//          if ((!sharing) || (t->share_from[j] == nullptr)) {
            if (sense == MIN)
              rhs[j] = CPX_INFBOUND;
            else
              rhs[j] = -CPX_INFBOUND;
//          } else {
//#ifdef DEBUG
//            debug_mutex.lock();
//            std::cout << "Thread " << t->id << " setting rhs[" << j;
//            std::cout << "] to " << *t->share_from[j] << std::endl;
//            debug_mutex.unlock();
//#endif
//            rhs[j] = *t->share_from[j];
//          }
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
        if (sense == MIN) {
          rhs[objective] = max[objective]-1;
          max[objective] = (int) -CPX_INFBOUND;
        } else {
          rhs[objective] = min[objective]+1;
          min[objective] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth_level = 1;
        depth = t->perm(depth_level);
        onwalk = false;
      } else if (inflast && infcnt != objective_counter) {
        if (sense == MIN) {
          rhs[depth] = CPX_INFBOUND;
        } else {
          rhs[depth] = -CPX_INFBOUND;
        }
        depth_level++;
        depth = t->perm(depth_level);
        if (sense == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = true;
      } else if (!onwalk && infcnt != 1) {
        if (sense == MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      } else if (onwalk && infcnt != 1)  {
        depth_level = 1;
        depth = t->perm(depth_level);
        if (sense == MIN) {
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
  for(int i = 0; i < p.objcnt; ++i) {
    Locking_Vars *lv = t->locks[i];
    if (lv != nullptr) {
      lv->num_running_threads--;
    }
  }
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  total_time = start.tv_sec + start.tv_nsec/1e9 - total_time;
  std::cout << "Thread " << t->id << " used " << cplex_time << "s in cplex";
  std::cout << ", waited for " << wait_time << "s";
  std::cout << " and " << total_time << "s overall." << std::endl;
#endif
  all.merge(s);
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
