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
    std::mutex status_mutex;
    std::mutex ready_mutex;
    std::atomic<int> num_running_threads;
    std::condition_variable cv;
    std::condition_variable ready_cv;

    Locking_Vars(int num_running_threads_);
};

#endif /* LOCKINGVARS_H */
