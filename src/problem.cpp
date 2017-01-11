#include "problem.h"
#include <cstring>

Problem::Problem(const char * filename, int cplex_threads):
      objcnt(0), cplex_threads(cplex_threads),
      mip_tolerance(1e-4), filename_(filename)
{
  filetype = UNKNOWN;
  int len = strlen(filename);
  if ((len > 3) && ('.' == filename[len-3]) &&
      ('l' == filename[len-2] && 'p' == filename[len-1])) {
    filetype = LP;
  } else if ((len > 4) && ('.' == filename[len-4]) &&
      ('m' == filename[len-3] && 'o' == filename[len-2] &&
      'p' == filename[len-1])) {
    filetype = MOP;
  }
}

