#include <list>

#ifdef DEBUG
#include <iostream>
#endif

#ifndef CPX_INFBOUND
#include <climits>
#define INFBOUND INT_MAX
#else
#define INFBOUND CPX_INFBOUND
#endif

#include "cluster.h"
#include "sense.h"
#include "thread.h"
#include "lockingvars.h"
#include "symgroup.h" // S[n].size()

Cluster::Cluster(int nThreads, int nObj, Sense sense, bool spread_threads,
    int nObjLeft, int * ordering, int ** share_to_, int ** share_from_,
    int ** share_bounds_, int ** share_limit_, std::list<Thread*> & threads,
    Locking_Vars ** locks) {
  if (nThreads == 1) {
    // Build perm
    int * perm = new int[nObj];
    int i;
    for(i = 0; i < nObj - nObjLeft; ++i) {
      perm[nObj - i - 1] = nObj - 1 - ordering[i];
    }
    // Fill in rest, but don't use values already used
    int mapsto = nObj - 1;
    while (i < nObj) {
      bool found = false;
      for(int k = 0; k < i; ++k) {
        if (perm[nObj - k - 1] == mapsto) {
          found = true;
          break;
        }
      }
      if (!found) {
        perm[nObj - 1 - i] = mapsto;
        ++i;
      }
      --mapsto;
    }
#ifdef DEBUG
    std::cout << "Thread " << threads.size() << " using permutation ";
    std::cout << perm[0];
    for(int d = 1; d < nObj; ++d) {
      std::cout << "," << perm[d];
    }
    std::cout << std::endl;
#endif
    // Build thread
    Thread *t = new Thread(threads.size(), nObj, perm, share_to_, share_from_,
        share_bounds_, share_limit_, locks);
    threads.push_back(t);
//    threads.push_back( new Thread(threads.size(), nObj, perm, shared_limits));
  } else {
    int ** share_to = new int*[nObj];
    int ** share_from = new int*[nObj];
    int ** share_limit = new int*[nObj];
    int ** share_bounds = new int*[nObj];
    for(int j = 0; j < nObj; ++j) {
      share_to[j] = share_to_[j];
      share_bounds[j] = share_bounds_[j];
      share_from[j] = share_from_[j];
      share_limit[j] = share_limit_[j];
    }
    // Work out which objectives don't occur in ordering
    int ind = 0;
    int * objLeft = new int[nObj];
    for(int j = 0; j < nObj; ++j) {
      bool found = false;
      for(int k = 0; k < nObj - nObjLeft; ++k) {
        if (ordering[k] == j) {
          found = true;
          break;
        }
      }
      if (! found) {
        objLeft[ind] = j;
        ind++;
      }
    }
    int ** new_shares = new int*[nObj] {nullptr};
    int ** new_bounds = new int*[nObj] {nullptr};
    int ** new_limit = new int*[nObj] {nullptr};
    int num_sub_clusters = (ind < nThreads) ? ind : nThreads;
    for(int j = 0; j < num_sub_clusters; ++j) {
      int pos = nObj - objLeft[j] - 1;
      new_shares[pos] = new int;
      new_bounds[pos] = new int;
      new_limit[pos] = new int;
      if (sense == MIN) {
        *new_shares[pos] = INFBOUND;
        *new_bounds[pos] = -INFBOUND;
        *new_limit[pos] = INFBOUND;
      } else {
        *new_shares[pos] = -INFBOUND;
        *new_bounds[pos] = INFBOUND;
        *new_limit[pos] = -INFBOUND;
      }
    }

    if (spread_threads) {
      // Split up threads
      int perCluster = nThreads / nObjLeft;
      int withExtra = nThreads % nObjLeft;
      int i;
      for(i = 0; i < withExtra ; ++i) {
        ordering[ nObj - nObjLeft ] = objLeft[i];
        int pos = nObj - objLeft[i] - 1;
        locks[pos] = new Locking_Vars(perCluster + 1);
        int * old_share_to = share_to[pos];
        int * old_share_bounds = share_bounds[pos];
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = nObj - objLeft[j] - 1;
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
            share_bounds[obj] = new_bounds[obj];
          } else {
            share_limit[obj] = new_limit[obj];
            share_from[obj] = new_shares[obj];
          }
        }
        Cluster(perCluster + 1, nObj, sense, spread_threads, nObjLeft - 1,
            ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
        share_to[pos] = old_share_to;
        share_bounds[pos] = old_share_bounds;
      }
      if (perCluster > 0) {
        for( ; i < nObjLeft ; ++i) {
          ordering[ nObj - nObjLeft ] = objLeft[i];
          int pos = nObj - objLeft[i] - 1;
          int * old_share_to = share_to[pos];
          int * old_share_bounds = share_bounds[pos];
          locks[pos] = new Locking_Vars(perCluster);
          for(int j = 0; j < nObj; ++j) {
            int obj = nObj - objLeft[j] - 1;
            if (obj == pos) {
              share_to[obj] = new_shares[obj];
              share_bounds[obj] = new_bounds[obj];
            } else {
              share_from[obj] = new_shares[obj];
              share_limit[obj] = new_limit[obj];
            }
          }
          Cluster(perCluster, nObj, sense, spread_threads, nObjLeft - 1,
              ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
          share_to[pos] = old_share_to;
          share_bounds[pos] = old_share_bounds;
        }
      }
    } else { // grouping threads "near" each other
      int threads_remaining = nThreads;
      int i = 0; // Which objective we're currently assigning threads to.
      while (threads_remaining > 0) {
        int threads_to_use = (S[nObjLeft-1].size() > threads_remaining) ? threads_remaining : S[nObjLeft-1].size();
        ordering[ nObj - nObjLeft ] = objLeft[i];
        int pos = nObj - objLeft[i] - 1;
        locks[pos] = new Locking_Vars(threads_to_use);
        int * old_share_to = share_to[pos];
        int * old_share_bounds = share_bounds[pos];
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = nObj - objLeft[j] - 1;
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
            share_bounds[pos] = new_bounds[obj];
            share_limit[obj] = nullptr;
          } else {
            share_from[obj] = new_shares[obj];
            share_limit[obj] = new_limit[obj];
          }
        }
        Cluster(threads_to_use, nObj, sense, spread_threads, nObjLeft - 1,
            ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
        threads_remaining -= threads_to_use;
        share_to[pos] = old_share_to;
        share_bounds[pos] = old_share_bounds;
        // Need to reset share_limits ?
        ++i;
      }
    }
    delete[] objLeft;
    delete[] new_shares;
    delete[] new_limit;
    delete[] share_to;
    delete[] share_from;
    delete[] share_limit;
    delete[] share_bounds;
  }
}
