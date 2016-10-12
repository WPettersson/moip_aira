#include "result.h"

/* Compare two solutions of type Result
 * True if a ≥ b, False otherwise.
 * This might seem backwards, but the sort is only used for final display of
 * results, so this sorts the results in descending order.
 * Note that infeasible > feasible, but we don't show infeasible results so it
 * doesn't really matter. */
bool Result::operator<(const Result * other) const {
  if (this->infeasible) {
      return true;
  }
  if (other->infeasible) {
    return false;
  }
  for (int i = 0; i < objective_count; i++) {
    if (this->result[i] < other->result[i]) {
      return false;
    }
    else if (this->result[i] > other->result[i]) {
      return true;
    }
  }
  return false;
}
