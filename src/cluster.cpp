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
    for(i = 0; i < nObj; ++i) {
      perm[i] = ordering[i];
    }
#ifdef DEBUG
    std::cout << "Thread " << threads.size() << " using permutation ";
    std::cout << perm[0];
    for(int d = 1; d < nObj; ++d) {
      std::cout << "," << perm[d];
    }
    std::cout << std::endl;
#endif
    // If there is only one objective left to add, then at the previous stage
    // we must've had 2 objectives left, and 2 threads to use (else the Thread
    // object would be constructed there), so this thread would have a partner.
    bool has_partner = (nObjLeft == 1);
    // Build thread
    Thread *t = new Thread(threads.size(), nObj, perm, share_to_, share_from_,
        share_bounds_, share_limit_, locks, has_partner);
    threads.push_back(t);
//    threads.push_back( new Thread(threads.size(), nObj, perm, shared_limits));
  } else {
    int *my_ordering = new int[nObj];
#ifdef DEBUG
    std::cout << "Working with permutation ";
    std::cout << ordering[0];
    for(int d = 1; d < nObj; ++d) {
      std::cout << "," << ordering[d];
    }
    std::cout << ", " << nThreads << " threads left,";
    std::cout << " and " << nObjLeft << " objectives left.";
    std::cout << std::endl;
#endif
    for(int i = 0; i < nObj; ++i) {
      my_ordering[i] = ordering[i];
    }
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
    // Work out how many different objectives will occur
    // in position "nObjLeft - 1"
    int ind = nObjLeft;
    int ** new_shares = new int*[nObj] {nullptr};
    int ** new_bounds = new int*[nObj] {nullptr};
    int ** new_limit = new int*[nObj] {nullptr};
    int num_sub_clusters = (ind < nThreads) ? ind : nThreads;
    for(int j = 0; j < num_sub_clusters; ++j) {
      int pos = my_ordering[nObjLeft - j - 1];
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
        int pos = my_ordering[nObjLeft-1];
        locks[pos] = new Locking_Vars(perCluster + 1);
        int * old_share_to = share_to[pos];
        int * old_share_bounds = share_bounds[pos];
        int * old_share_limit = share_limit[pos];
        int ** old_from = new int*[nObj];
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          old_from[obj] = share_from[obj];
        }
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
            share_bounds[obj] = new_bounds[obj];
            share_limit[obj] = new_limit[obj];
          } else {
            share_from[obj] = new_shares[obj];
          }
        }
        Cluster(perCluster + 1, nObj, sense, spread_threads, nObjLeft - 1,
            my_ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
        share_to[pos] = old_share_to;
        share_bounds[pos] = old_share_bounds;
        share_limit[pos] = old_share_limit;
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          share_from[obj] = old_from[obj];
        }
        // Rotate ordering
        int first = my_ordering[0];
        for (int c = 0; c < nObjLeft-1; ++c) {
          my_ordering[c] = my_ordering[c+1];
        }
        my_ordering[nObjLeft - 1] = first;
        locks[pos] = nullptr;
      }
      if (perCluster > 0) {
        for( ; i < nObjLeft ; ++i) {
          int pos = my_ordering[nObjLeft-1];
          int * old_share_to = share_to[pos];
          int * old_share_bounds = share_bounds[pos];
          int * old_share_limit = share_limit[pos];
          locks[pos] = new Locking_Vars(perCluster);
          int ** old_from = new int*[nObj];
          for(int j = 0; j < nObjLeft; ++j) {
            int obj = my_ordering[j];
            old_from[obj] = share_from[obj];
          }
          for(int j = 0; j < nObjLeft; ++j) {
            int obj = my_ordering[j];
            if (obj == pos) {
              share_to[obj] = new_shares[obj];
              share_bounds[obj] = new_bounds[obj];
              share_limit[obj] = new_limit[obj];
            } else {
              share_from[obj] = new_shares[obj];
            }
          }
          Cluster(perCluster, nObj, sense, spread_threads, nObjLeft - 1,
              my_ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
          share_to[pos] = old_share_to;
          share_bounds[pos] = old_share_bounds;
          share_limit[pos] = old_share_limit;
          for(int j = 0; j < nObjLeft; ++j) {
            int obj = my_ordering[j];
            share_from[obj] = old_from[obj];
          }
          // Rotate ordering
          int first = my_ordering[0];
          for (int c = 0; c < nObjLeft-1; ++c) {
            my_ordering[c] = my_ordering[c+1];
          }
          my_ordering[nObjLeft - 1] = first;
          locks[pos] = nullptr;
        }
      }
    } else { // grouping threads "near" each other
      int threads_remaining = nThreads;
      while (threads_remaining > 0) {
        int threads_to_use = (S[nObjLeft-1].size() > threads_remaining) ? threads_remaining : S[nObjLeft-1].size();
        int pos = my_ordering[nObjLeft-1];
        locks[pos] = new Locking_Vars(threads_to_use);
        int * old_share_to = share_to[pos];
        int * old_share_bounds = share_bounds[pos];
        int * old_share_limit = share_limit[pos];
        int ** old_from = new int*[nObj];
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          old_from[obj] = share_from[obj];
        }
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          if (obj == pos) {
            share_to[obj] = new_shares[obj];
            share_bounds[pos] = new_bounds[obj];
            share_limit[obj] = new_limit[obj];
          } else {
            share_from[obj] = new_shares[obj];
          }
        }
        Cluster(threads_to_use, nObj, sense, spread_threads, nObjLeft - 1,
            my_ordering, share_to, share_from, share_bounds, share_limit, threads, locks);
        threads_remaining -= threads_to_use;
        share_to[pos] = old_share_to;
        share_bounds[pos] = old_share_bounds;
        share_limit[pos] = old_share_limit;
        for(int j = 0; j < nObjLeft; ++j) {
          int obj = my_ordering[j];
          share_from[obj] = old_from[obj];
        }
        // Rotate ordering
        int first = my_ordering[0];
        for (int c = 0; c < nObjLeft-1; ++c) {
          my_ordering[c] = my_ordering[c+1];
        }
        my_ordering[nObjLeft - 1] = first;
        locks[pos] = nullptr;
      }
    }
    delete[] my_ordering;
    delete[] new_shares;
    delete[] new_limit;
    delete[] share_to;
    delete[] share_from;
    delete[] share_limit;
    delete[] share_bounds;
  }
}
