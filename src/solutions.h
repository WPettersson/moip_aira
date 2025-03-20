#ifndef SOLUTIONS_H
#define SOLUTIONS_H

#include <list>
#include <mutex>
#include "result.h"
#include "sense.h"

class Solutions {

  public:
    Solutions(int numObjectives);
    ~Solutions();
    const Result * find(const double *ip, const Sense sense) const;
    void insert(const double *lp, const int *result, const bool infeasible);
    void merge(Solutions& other);
    void sort();

    // Iterator functionality
    std::list<Result*>::const_iterator begin() const;
    std::list<Result*>::const_iterator end() const;

    /*
     * Sort the list of results, and remove any duplicates.
     */
    void sort_unique();

  private:
    int objective_count;
    // Note that a pointer to a Result object should only ever exist in a
    // single list. We allocate objects inside insert(), and only free them in
    // the Solutions destructor.
    std::list<Result*> store_;
    std::mutex mutex;
};

inline Solutions::Solutions(int numObjectives) : objective_count(numObjectives)
{ }


inline void Solutions::merge(Solutions& other) {
  std::unique_lock<std::mutex> lk(mutex);
  store_.splice(store_.begin(), other.store_);
}

inline std::list<Result *>::const_iterator Solutions::begin() const {
  return store_.begin();
}

inline std::list<Result *>::const_iterator Solutions::end() const {
  return store_.end();
}

inline void Solutions::sort_unique() {
  store_.sort( [](const Result *a, const Result *b) {return (*a) < (*b);} );
  store_.unique( [](const Result *a, const Result *b) {return (*a) ==( *b);} );
}

#endif /* SOLUTIONS_H */
