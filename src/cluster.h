#ifndef CLUSTER_H
#define CLUSTER_H

#include <list>

#include "thread.h"
#include "sense.h"
#include "lockingvars.h"

class Cluster {
  public:
    Cluster(int nThreads, int nObj, Sense sense, bool spread_threads, int
        nObjLeft, int * ordering, int ** share_to, int ** share_from, int **
        share_bounds, int ** share_limit, std::list<Thread*> & threads,
        Locking_Vars ** locks, bool dive);

};

#endif /* CLUSTER_H */
