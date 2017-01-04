#ifndef PROBLEM_H
#define PROBLEM_H

#include "sense.h"

enum filetype_t { UNKNOWN, LP, MOP };

class Problem {
  public:
    int objcnt; // Number of objectives
    double* rhs;
    int** objind; // Objective indices
    double** objcoef; // Objective coefficients
    Sense objsen; // Objective sense. Note that all objectives must have the same
                // sense (i.e., either all objectives are to be minimised, or
                // all objectives are to be maximised).
    int* conind;
    char* consense;

    /**
     * cplex_threads is a problem parameter, because each thread needs to know
     * the number of cplex threads to use.
     */
    int cplex_threads;


    const char* filename_;
    filetype_t filetype;

    const char* filename();

    Problem(const char* filename, int cplex_threads);
    ~Problem();

};

inline const char * Problem::filename() {
  return filename_;
}

inline Problem::~Problem() {
  // If objcnt == 0, then no problem has been assigned and no memory allocated
  if (objcnt == 0)
    return;
  for(int j = 0; j < objcnt; ++j) {
    delete[] objind[j];
    delete[] objcoef[j];
  }
  delete[] objind;
  delete[] objcoef;
  delete[] rhs;
  delete[] conind;
  delete[] consense;
}
#endif /* PROBLEM_H */
