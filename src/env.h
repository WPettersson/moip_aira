#ifndef ENV_H
#define ENV_H

#include <ilcplex/cplex.h>

class Env {
  public:
   CPXENVptr env;
   CPXLPptr lp;
};
#endif /* ENV_H */
