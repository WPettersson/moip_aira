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
    std::condition_variable rv;
    std::atomic<bool> changed;


    Locking_Vars(int num_running_threads_);
    void reset_num_running_threads();
    int max_threads() const;
  private:
    int max_running_threads;
};

inline void Locking_Vars::reset_num_running_threads() {
  num_running_threads = max_running_threads;
}

inline int Locking_Vars::max_threads() const {
  return max_running_threads;
}

inline Locking_Vars::Locking_Vars(int num_running_threads_) :
    num_running_threads(num_running_threads_), changed(false),
    max_running_threads(num_running_threads_) {
}


#endif /* LOCKINGVARS_H */
