#include "thread.h"
#include "symgroup.h"
#include "lockingvars.h"

#ifdef DEBUG
#include <iostream>
#endif

Thread::Thread(int id_, int nObj_, const int * perm__, int ** share_to_, int ** share_from_,
    int ** share_bounds_, int ** share_limit_, Locking_Vars ** locks_, bool partnered_) :
  nObj(nObj_), id(id_), partnered(partnered_), split_start(0)//, split_stop(0)
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
#ifdef DEBUG
    std::cout << "locks[" << i << "] = " << locks[i] << std::endl;
#endif
  }
  perm_ = new int[nObj];
  for(int i = 0; i < nObj; ++i) {
    perm_[i] = perm__[i];
  }
}

void Thread::setNumObjectives(int newObj) {
  int * newPerm = new int[newObj];
  for(int i = 0; i < nObj; ++i) {
    newPerm[i] = perm_[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    newPerm[i] = i;
  }
  delete[] perm_;
  perm_ = newPerm;

  int ** new_share_to = new int*[newObj];
  for(int i = 0; i < nObj; ++i) {
    new_share_to[i] = share_to[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    new_share_to[i] = nullptr;
  }
  delete[] share_to;
  share_to = new_share_to;

  int ** new_share_from = new int*[newObj];
  for(int i = 0; i < nObj; ++i) {
    new_share_from[i] = share_from[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    new_share_from[i] = nullptr;
  }
  delete[] share_from;
  share_from = new_share_from;

  int ** new_share_bounds = new int*[newObj];
  for(int i = 0; i < nObj; ++i) {
    new_share_bounds[i] = share_bounds[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    new_share_bounds[i] = nullptr;
  }
  delete[] share_bounds;
  share_bounds = new_share_bounds;

  int ** new_share_limit = new int*[newObj];
  for(int i = 0; i < nObj; ++i) {
    new_share_limit[i] = share_limit[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    new_share_limit[i] = nullptr;
  }
  delete[] share_limit;
  share_limit = new_share_limit;

  Locking_Vars** new_locks = new Locking_Vars*[newObj];
  for(int i = 0; i < nObj; ++i) {
    new_locks[i] = locks[i];
  }
  for (int i = nObj; i < newObj; ++i) {
    new_locks[i] = nullptr;
  }
  delete[] locks;
  locks = new_locks;

  nObj = newObj;
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
