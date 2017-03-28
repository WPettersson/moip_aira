#ifndef THREAD_H
#define THREAD_H

#include "lockingvars.h"

class Thread {
  private:
    int nObj_;
  public:
    int id;
    int perm(int) const;
    int nObj() const;
    int ** share_to;
    int ** share_from;
    int ** share_bounds;
    int ** share_limit;
    Locking_Vars ** locks;

    bool partnered;

    double split_start;
    double split_stop;

    void setNumObjectives(int newObj);

    Thread(int id_, int nObj, const int * perm_, int ** share_to_,
        int ** share_from_, int ** share_bounds_, int ** share_limit_, 
        Locking_Vars ** locks_, bool partnered);
    Thread(int id_, int nObj, int totalObjCount, double split_start_, double split_stop_);

  private:
    int * perm_;
};

inline int Thread::nObj() const {
  return nObj_;
}

inline int Thread::perm(int ind) const {
  return perm_[ind];
}

#endif /* THREAD_H */
