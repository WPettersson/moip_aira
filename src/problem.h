#ifndef PROBLEM_H
#define PROBLEM_H

#include "sense.h"

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

    ~Problem();

};


inline Problem::~Problem() {
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
