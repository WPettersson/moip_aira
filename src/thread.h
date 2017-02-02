#ifndef THREAD_H
#define THREAD_H

#include "lockingvars.h"

class Thread {
  public:
    int id;
    int perm(int) const;
    int ** share_to;
    int ** share_from;
    Locking_Vars ** locks;

    double split_start;
    double split_stop;


    Thread(int id_, int nObj, const int * perm_, int ** share_to_,
        int ** share_from_, Locking_Vars ** locks_);
    Thread(int id_, int nObj, double split_start_, double split_stop_);

  private:
    int * perm_;
};

inline int Thread::perm(int ind) const {
  return perm_[ind];
}

#endif /* THREAD_H */
