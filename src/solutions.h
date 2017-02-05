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

    // Sorting
    static void sort( bool (*cmp)(const Result * a, const Result * b) );

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
  // Note that only one list is allowed to own a pointer at any one time.
  // We aren't used unique_ptr to avoid overheads, we are trusting that our
  // coding is correct.
  while (other.store_.size() > 0) {
    store_.push_back(other.store_.front());
    other.store_.pop_front();
  }
}

inline std::list<Result *>::const_iterator Solutions::begin() const {
  return store_.begin();
}

inline std::list<Result *>::const_iterator Solutions::end() const {
  return store_.end();
}

inline void Solutions::sort() {
  store_.sort( [](const Result *a, const Result *b) {return (*a) < (*b);} );
}

#endif /* SOLUTIONS_H */
