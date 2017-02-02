#ifndef CLUSTER_H
#define CLUSTER_H

#include <list>

#include "thread.h"
#include "sense.h"
#include "lockingvars.h"

class Cluster {
  public:
    Cluster(int nThreads, int nObj, Sense sense, int nObjLeft, int * ordering, int ** share_to, int ** share_from, std::list<Thread*> & threads, Locking_Vars ** locks);

};

#endif /* CLUSTER_H */
