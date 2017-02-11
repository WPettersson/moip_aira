#include "result.h"

/* Compare two solutions of type Result
 * True if a ≥ b, False otherwise.
 * This might seem backwards, but the sort is only used for final display of
 * results, so this sorts the results in descending order.
 * Note that infeasible > feasible, but we don't show infeasible results so it
 * doesn't really matter. */
bool operator<(const Result& lhs, const Result& rhs) {
  if (lhs.infeasible) {
      return true;
  }
  if (rhs.infeasible) {
    return false;
  }
  if (lhs.objective_count < rhs.objective_count)
    return false;
  if (lhs.objective_count > rhs.objective_count)
    return true;
  for (int i = 0; i < lhs.objective_count; i++) {
    if (lhs.result[i] < rhs.result[i]) {
      return false;
    }
    else if (lhs.result[i] > rhs.result[i]) {
      return true;
    }
  }
  return false;
}

bool operator==(const Result& lhs, const Result& rhs) {
  if (lhs.infeasible != rhs.infeasible) {
    return false;
  }
  if (rhs.infeasible) {
    return true;
  }
  if (lhs.objective_count != rhs.objective_count)
    return false;
  for (int i = 0; i < lhs.objective_count; i++) {
    if (lhs.result[i] != rhs.result[i]) {
      return false;
    }
  }
  return true;
}

std::ostream& operator<<(std::ostream& out, const Result& r) {
  if (r.infeasible) {
    out << "Infeasible";
    return out;
  }
  out << " Feasible : ";
  for (int i = 0; i < r.objective_count; ++i) {
    out << r.result[i];
    if (i < r.objective_count-1)
      out <<",";
  }
  return out;
}
