#ifndef RESULT_H
#define RESULT_H

/* List node. */
class Result {
  public:
    int *result;
    int objective_count ;
    bool infeasible;
    double *ip;
    bool operator<(const Result *other) const;
};

#endif /* RESULT_H */
