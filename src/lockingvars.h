#ifndef LOCKINGVARS_H
#define LOCKINGVARS_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <iostream>
#include <ostream>
#include <vector>

#include "threadstate.h"


/**
 * Any pair of threads that are sharing solutions need to synchronise limits on
 * the other objectives. This struct contains the necessary variables to handle
 * such synchronisation.
 */
class Locking_Vars {
  public:
    /**
     * A mutex used to synchronise threads.
     */
    std::mutex status_mutex;
    /**
     * Used for synchronising the sharing of bounds.
     */
    std::condition_variable cv;
    /**
     * Used for synchronising when subproblems are complete.
     */
    std::condition_variable cv_subproblem_complete;


    void add_state(std::atomic<ThreadState> const *state);
    /**
     * Allows propogation of limits introduced from outside this cluster.
     *
     * Threads only look to one other thread for updated limits. As a result,
     * it can take multiple "rounds" for an update to propogate completely
     * through a cluster. This variable tracks this propogation.
     */
    std::atomic<bool> changed;
    /**
     * If one thread in the cluster finds a result, all threads need to know
     * this.
     */
    std::atomic<bool> found_any;


    /**
     * Basic constructor
     */
    Locking_Vars();
    /**
     * Check if all threads are done. That is, have they all either finished
     * this subproblem or reached a point where they can share bounds.
     */
    bool any_complete() const;
    /**
     * Check if all threads are complete. That is, have they found their solutions and exited.
     */
    bool all_done() const;
    /**
     * A list of thread statuses, so we can check the status of each thread that shares this lock.
     */
    std::vector<std::atomic<ThreadState> const *> states;
};

inline bool Locking_Vars::all_done() const {
  for(auto s: states) {
    if (*s == ThreadState::Running) {
      return false;
    }
  }
  return true;
}

inline bool Locking_Vars::any_complete() const {
  for(auto s: states) {
    if (*s == ThreadState::Complete) {
      return true;
    }
  }
  return false;
}


inline void Locking_Vars::add_state(std::atomic<ThreadState> const *state) {
  states.push_back(state);
}

/**
  * Basic constructor
  */
inline Locking_Vars::Locking_Vars() :
    changed(false), found_any(false) {
}


#endif /* LOCKINGVARS_H */
