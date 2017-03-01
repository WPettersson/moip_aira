#ifndef LOCKINGVARS_H
#define LOCKINGVARS_H

#include <atomic>
#include <condition_variable>
#include <mutex>

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
     * How many threads that share this Locking_Vars are currently running.
     * The last thread running is often used as a pseudo-master when
     * synchronising.
     */
    std::atomic<int> num_running_threads;
    /**
     * Used for synchronising as well.
     */
    std::condition_variable cv;
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
    Locking_Vars(int num_running_threads_);
    /**
     * Reset the number of running threads
     */
    void reset_num_running_threads();
    /**
     * The number of threads in this cluster.
     */
    int max_threads() const;
  private:
    /**
     * The number of threads in this cluster.
     */
    int max_running_threads;
};

/**
  * Reset the number of running threads
  */
inline void Locking_Vars::reset_num_running_threads() {
  num_running_threads = max_running_threads;
}

/**
  * The number of threads in this cluster.
  */
inline int Locking_Vars::max_threads() const {
  return max_running_threads;
}

/**
  * Basic constructor
  */
inline Locking_Vars::Locking_Vars(int num_running_threads_) :
    num_running_threads(num_running_threads_), changed(false),
    found_any(false), max_running_threads(num_running_threads_) {
}


#endif /* LOCKINGVARS_H */
