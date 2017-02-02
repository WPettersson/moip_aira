#include <list>

#ifdef DEBUG
#include <iostream>
#endif

#ifndef CPX_INFBOUND
#include <climits>
#define CPX_INFBOUND INT_MAX
#endif

#include "cluster.h"
#include "sense.h"
#include "thread.h"
#include "lockingvars.h"
#include "symgroup.h" // S[n].size()

Cluster::Cluster(int nThreads, int nObj, Sense sense, bool spread_threads, int nObjLeft, int * ordering, int ** share_to, int ** share_from, std::list<Thread*> & threads, Locking_Vars ** locks) {
  if (nThreads == 1) {
#ifdef DEBUG
    std::cout << "Thread " << threads.size() << " has ordering ";
    std::cout << ordering[0];
    for(int d = 1; d < nObj - nObjLeft; ++d) {
      std::cout << "," << ordering[d];
    }
    std::cout << std::endl;
#endif
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
    Thread *t = new Thread(threads.size(), nObj, perm, share_to, share_from, locks);
    threads.push_back(t);
//    threads.push_back( new Thread(threads.size(), nObj, perm, shared_limits));
  } else {
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
    int num_sub_clusters = (ind < nThreads) ? ind : nThreads;
    for(int j = 0; j < num_sub_clusters; ++j) {
      int pos = nObj - objLeft[j] - 1;
      new_shares[pos] = new int;
      if (sense == MIN) {
        *new_shares[pos] = CPX_INFBOUND;
      } else {
        *new_shares[pos] = -CPX_INFBOUND;
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
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = objLeft[j];
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
          } else {
            share_from[obj] = new_shares[obj];
          }
        }
        Cluster(perCluster + 1, nObj, sense, spread_threads, nObjLeft - 1, ordering, share_to, share_from,
            threads, locks);
        share_to[pos] = nullptr;
      }
      if (perCluster > 0) {
        for( ; i < nObjLeft ; ++i) {
          ordering[ nObj - nObjLeft ] = objLeft[i];
          int pos = nObj - objLeft[i] - 1;
          locks[pos] = new Locking_Vars(perCluster);
          for(int j = 0; j < nObj; ++j) {
            int obj = objLeft[j];
            if (obj == pos) {
              share_to[obj] = new_shares[obj];
            } else {
              share_from[obj] = new_shares[obj];
            }
          }
          Cluster(perCluster, nObj, sense, spread_threads, nObjLeft - 1, ordering, share_to, share_from,
              threads, locks);
          share_to[pos] = nullptr;
        }
      }
    } else { // grouping threads "near" each other
      int threads_remaining = nThreads;
      int i = 0; // Which objective we're currently assigning threads to.
      while (threads_remaining > 0) {
        int threads_to_use = (S[nObjLeft-1].size() > threads_remaining) ? threads_remaining : S[nObjLeft-1].size();
        ordering[ nObj - nObjLeft ] = objLeft[i];
        std::cout << "Assigning " << threads_to_use << " threads." << std::endl;
        std::cout << "working on " << (nObj - nObjLeft) << std::endl;
        std::cout << "ordering now is " << ordering[nObj - nObjLeft] << std::endl;
        int pos = nObj - objLeft[i] - 1;
        locks[pos] = new Locking_Vars(threads_to_use);
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = objLeft[j];
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
          } else {
            share_from[obj] = new_shares[obj];
          }
        }
        Cluster(threads_to_use, nObj, sense, spread_threads, nObjLeft - 1, ordering, share_to, share_from,
            threads, locks);
        threads_remaining -= threads_to_use;
        share_to[pos] = nullptr;
        ++i;
      }
    }
    delete[] objLeft;
    delete[] new_shares;
  }
}
