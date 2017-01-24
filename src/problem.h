#ifndef PROBLEM_H
#define PROBLEM_H

#include "sense.h"
#include "env.h"

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

    double mip_tolerance;

    filetype_t filetype;

    const char* filename();

    Problem(const char* filename, Env& env);
    ~Problem();

  private:
    int read_lp_problem(Env& e);
    int read_mop_problem(Env& e);
    const char* filename_;

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
