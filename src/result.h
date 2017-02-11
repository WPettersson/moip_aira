#ifndef RESULT_H
#define RESULT_H

#include <iostream>


/* List node. */
class Result {
  public:
    int *result;
    int objective_count ;
    bool infeasible;
    double *ip;
    friend bool operator<(const Result& lhs, const Result& rhs);
    friend bool operator!=(const Result& lhs, const Result& rhs);
    friend bool operator==(const Result& lhs, const Result& rhs);

};

std::ostream& operator<<(std::ostream& out, const Result& r);


inline bool operator!=(const Result& lhs, const Result& rhs) {
  return !(lhs == rhs);
}

#endif /* RESULT_H */
