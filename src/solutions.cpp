
#include <cstring>

#include "result.h"
#include "solutions.h"

const Result * Solutions::find(const double *ip, const Sense sense) const {
  bool t1,t2,t3;
  for (auto& res: store_) {

    t1 = true;
    t2 = false;
    t3 = true;
    if (sense == MIN) {
      for (int i = 0; i < objective_count; i++){
        /* t1: All values of candidate must be >= than ip */
        if (res->ip[i] < ip[i]) {
          t1 = false;
          break;
        }
        /* t2: At least one inequality must be strict */
        if (res->ip[i] > ip[i]) {
          t2 = true;
        }
        /* t3: All values of candidate res->lt must be <= than ip */
        if (!res->infeasible) {
          if (res->result[i] > ip[i]) {
            t3 = false;
            break;
          }
        }
      }
    }
    else {
      for (int i = 0; i < objective_count; i++){
        /* t1: All values of candidate must be <= than ip */
        if (res->ip[i] > ip[i]) {
          t1 = false;
          break;
        }
        /* t3: All values of candidate result must be >= than ip */
        if (!res->infeasible) {
          if (res->result[i] < ip[i]) {
            t3 = false;
            break;
          }
        }
      }
    }
    /* If all conditions are met copy problem & solution and return */
    if (t1 && t3) {
      return res;
    }
  }
  return nullptr;
}
void Solutions::insert(const double *lp, const int *result,
    const bool infeasible) {
  Result * r = new Result;
  r->ip = new double[objective_count];
  memcpy(r->ip, lp, objective_count * sizeof(double));
  if (infeasible) {
    r->infeasible = true;
    r->result = nullptr;
  }
  else {
    r->infeasible = false;
    r->result = new int[objective_count];
    memcpy(r->result, result, objective_count * sizeof(int));
  }
  store_.push_back(r);
}

Solutions::~Solutions() {
  for(auto r: store_) {
    delete[] r->ip;
    if (r->result)
      delete[] r->result;
  }
}
