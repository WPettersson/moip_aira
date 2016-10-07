#ifndef SYMGROUP_H
#define SYMGROUP_H

class SymGroup {
  private:
    const int *data_;
    int N_;
  public:
    SymGroup(int N, const int *data);
    const int* operator[](int index) const;
    int size() const;
};

inline SymGroup::SymGroup(int N, const int *data) : data_(data), N_(N) {
}

inline const int* SymGroup::operator[](int index) const {
  return (data_ + N_ * index);
}

#include "symgroup_extern.h"

#endif /* SYMGROUP_H */
