#include "thread.h"
#include "symgroup.h"
#include "lockingvars.h"

#ifdef DEBUG
#include <iostream>
#endif

Thread::Thread(int id_, int nObj, const int * perm__, int ** share_to_, int ** share_from_,
    int ** share_bounds_, int ** share_limit_, Locking_Vars ** locks_) :
  id(id_), split_start(0)//, split_stop(0)
{
  split_stop = 0;
  share_to = new int*[nObj];
#ifdef DEBUG
  std::cout << "Thread: " << id << std::endl;
#endif
  for(int i = 0; i < nObj; ++i) {
    share_to[i] = share_to_[i];
#ifdef DEBUG
    std::cout << "share_to[" << i << "] = " << share_to[i] << std::endl;
#endif
  }
  share_from = new int*[nObj];
  for(int i = 0; i < nObj; ++i) {
    share_from[i] = share_from_[i];
#ifdef DEBUG
    std::cout << "share_from[" << i << "] = " << share_from[i] << std::endl;
#endif
  }
  share_bounds = new int*[nObj];
  for(int i = 0; i < nObj; ++i) {
    share_bounds[i] = share_bounds_[i];
#ifdef DEBUG
    std::cout << "share_bounds[" << i << "] = " << share_bounds[i] << std::endl;
#endif
  }
  share_limit = new int*[nObj];
  for(int i = 0; i < nObj; ++i) {
    share_limit[i] = share_limit_[i];
#ifdef DEBUG
    std::cout << "share_limit[" << i << "] = " << share_limit[i] << std::endl;
#endif
  }
  locks = new Locking_Vars* [nObj];
  for(int i = 0; i < nObj; ++i) {
    locks[i] = locks_[i];
  }
  perm_ = new int[nObj];
  for(int i = 0; i < nObj; ++i) {
    perm_[i] = perm__[i];
  }
}

Thread::Thread(int id_, int nObj, double split_start_, double split_stop_) :
  id(id_), share_to(nullptr), share_from(nullptr), locks(nullptr), split_start(split_start_),
  split_stop(split_stop_)
{
  perm_ = new int[nObj];
  for(int i = 0; i < nObj; ++i) {
    perm_[i] = i; // If splitting, all threads use a common permutation
  }
}
