#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "problem.h"
#include "env.h"
#include "errors.h"

Problem::Problem(const char * filename, Env& env):
      objcnt(0), mip_tolerance(1e-4), filename_(filename)
{
  filetype = UNKNOWN;
  int len = strlen(filename);
  if ((len > 3) && ('.' == filename[len-3]) &&
      ('l' == filename[len-2] && 'p' == filename[len-1])) {
    filetype = LP;
    read_lp_problem(env);
  } else if ((len > 4) && ('.' == filename[len-4]) &&
      ('m' == filename[len-3] && 'o' == filename[len-2] &&
      'p' == filename[len-1])) {
    filetype = MOP;
    read_mop_problem(env);
  }
}

int Problem::read_lp_problem(Env& e) {
  int status;
  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, filename());

  if (e.lp == NULL) {
    std::cerr << "Failed to create LP." << std::endl;
    return -ERR_CPLEX;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, filename(), NULL);
  if (status) {
    std::cerr << "Failed to read and copy the problem data." << std::endl;
    return -ERR_CPLEX;
  }

  /* Get last rhs and determine the number of objectives.*/
  int cur_numcols = CPXgetnumcols(e.env, e.lp);
  int cur_numrows = CPXgetnumrows(e.env, e.lp);
  int cur_numnz = CPXgetnumnz(e.env, e.lp);

  rhs = new double[cur_numrows];

  /* Get RHS of last objective - this tells us the objective count */
  status = CPXgetrhs (e.env, e.lp, rhs, cur_numrows-1, cur_numrows-1);

  if (status) {
    std::cerr << "Failed to get RHS." << std::endl;
    return -ERR_CPLEX;
  }

  objcnt = static_cast<int>(rhs[0]);

  /* Create a pair of multidimensional arrays to store the objective
  * coefficients and their indices indices first */
  objind = new int*[objcnt];
  for(int j = 0; j < objcnt; j++) {
    objind[j] = new int[cur_numcols];
    for(int i = 0; i < cur_numcols; i++){
      objind[j][i] = i;
    }
  }
  /* Now coefficients */
  objcoef = new double*[objcnt];
  for(int j = 0; j < objcnt; j++) {
    objcoef[j] = new double[cur_numcols];
    for (int i = 0; i < cur_numcols; ++i)
      objcoef[j][i] = 0;
  }

  /* Parse out the objectives working backwards from the last constraint */
  int * rmatbeg = new int[cur_numrows];
  int * rmatind = new int[cur_numnz];
  double * rmatval = new double[cur_numnz];
  int rmatspace = cur_numnz;

  int nzcnt, surplus;
  status = CPXgetrows (e.env, e.lp, &nzcnt, rmatbeg, rmatind, rmatval,
                      rmatspace, &surplus, cur_numrows-objcnt, cur_numrows-1);

  if (status) {
    std::cerr << "Couldn't get rows." << std::endl;
    return -ERR_CPLEX;
  }

  for (int j = 0; j < objcnt; j++) {
    int to;
    int from = rmatbeg[j];
    if (j == objcnt-1) {
      to = nzcnt-1;
    }
    else {
      to = rmatbeg[(j+1)] - 1;
    }
    for (int k = from; k <= to; k++) {
      objcoef[j][rmatind[k]] = rmatval[k];
    }
  }
  delete[] rmatbeg;
  delete[] rmatind;
  delete[] rmatval;

  /* Setup problem for solving */

  /* Resize rhs to fit only objective function constraints */
  delete[] rhs;
  rhs = new double[objcnt];

  /* Get objective sense */
  int cpx_sense = CPXgetobjsen(e.env, e.lp);
  objsen = (cpx_sense == CPX_MIN ? MIN : MAX);
  /* Set objective constraint sense and RHS */
  consense = new char[objcnt];
  for (int j = 0; j < objcnt; j++) {
    if (objsen == MIN) {
      consense[j] = 'L'; /* Set sense to <= */
      rhs[j] = CPX_INFBOUND;
    }
    else {
      consense[j] = 'G'; /* Set sense to >= */
      rhs[j] = -CPX_INFBOUND;
    }
  }

  /* Specify index of objective constraints */
  conind = new int[objcnt];
  for (int j = 0, k = objcnt; j < objcnt; j++, k--) {
    conind[j] = cur_numrows-k;
  }

  /* Set sense of objective constraints */
  status = CPXchgsense (e.env, e.lp, objcnt, conind, consense);
  if (status) {
    std::cerr << "Failed to change constraint sense" << std::endl;
    return -ERR_CPLEX;
  }

  /* Set rhs of objective constraints */
  status = CPXchgrhs (e.env, e.lp, objcnt, conind, rhs);
  if (status) {
    std::cerr << "Failed to change constraint rhs" << std::endl;
    return -ERR_CPLEX;
  }
  return 0;
}


int Problem::read_mop_problem(Env& e) {
  int status;
  /* Create the problem, using the filename as the problem name */
  e.lp = CPXcreateprob(e.env, &status, filename());

  if (e.lp == NULL) {
    std::cerr << "Failed to create problem." << std::endl;
    return -ERR_CPLEX;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(e.env, e.lp, filename(), NULL);
  if (status) {
    std::cerr << "Failed to read and copy the problem data." << std::endl;
    return -ERR_CPLEX;
  }


  size_t cur_numcols = CPXgetnumcols(e.env, e.lp);
  size_t cur_numrows = CPXgetnumrows(e.env, e.lp);


  char** colNames = new char*[cur_numcols];
  size_t sz = cur_numcols * 1024;
  char* store = new char[sz];
  int surplus;
  status = CPXgetcolname(e.env, e.lp, colNames, store, sz, &surplus, 0, cur_numcols-1);
  if (status) {
    std::cerr << "Error retrieving column names." << std::endl;
    return -ERR_CPLEX;
  }


  // Iterate through file until we get to ROWS
  // count number of objectives
  // store names
  // break when we get COLUMNS
  // read a b c
  // if b == objective_name[i] objcoeff[i][ name_to_index[a]] = int(c)
  std::fstream mop_file(filename(), std::ios::in);
  std::string line;
  // Find "ROWS" line
  while (std::getline( mop_file, line)) {
    if ("ROWS" == line) {
      break;
    }
  }
  // Read objective names
  std::vector<std::string> objNames;
  while (std::getline( mop_file, line)) {
    // Objective functions have N as the second character
    if ('N' != line[1]) {
      break;
    }
    std::stringstream ss(line);
    std::string n, name;
    ss >> n; // Skip "N"
    ss >> name;
    objNames.push_back(name);
  }
  objcnt = static_cast<int>(objNames.size());
  int new_nzcnt = objcnt * cur_numcols;
  double * newRowMatVal = new double[new_nzcnt] {0};
  int * newRowMatInd = new int[new_nzcnt];

  // We will store each coefficient, in order, so set matrix up
  for(int i = 0; i < objcnt; ++i) {
    for(size_t j = 0; j < cur_numcols; ++j) {
      newRowMatInd[i*cur_numcols + j] = j;
    }
  }
  int * newRowMatBeg = new int[objNames.size()];
  for(int i = 0; i < objcnt; ++i) {
    newRowMatBeg[i] = i*cur_numcols;
  }

  /* Create a pair of multidimensional arrays to store the objective
  * coefficients and their indices indices first */
  objind = new int*[objcnt];
  for(int j = 0; j < objcnt; j++) {
    objind[j] = new int[cur_numcols];
    for(size_t i = 0; i < cur_numcols; i++){
      objind[j][i] = i;
    }
  }
  /* Now coefficients */
  objcoef = new double*[objcnt];
  for(int j = 0; j < objcnt; j++) {
    objcoef[j] = new double[cur_numcols];
    for (size_t i = 0; i < cur_numcols; ++i) {
      objcoef[j][i] = 0;
    }
  }
  // Find columns
  while (std::getline( mop_file, line)) {
    if ("COLUMNS" == line) {
      break;
    }
  }

  /* Parse out the objectives working backwards from the last constraint */
  while (std::getline( mop_file, line)) {
    std::stringstream ss(line);
    std::string col, obj;
    signed int val;
    ss >> col;
    ss >> obj;
    ss >> val;
    size_t colInd;
    int objInd;
    for(colInd = 0; colInd < cur_numcols; ++colInd) {
      if ( col.compare(colNames[colInd]) == 0) {
        break;
      }
    }
    // Skip if we don't recognise because why not?
    if (colInd == cur_numcols)
      continue;

    for(objInd = 0; objInd < objcnt; ++objInd) {
      if (objNames[objInd] == obj) {
        break;
      }
    }
    if (objInd == objcnt) {
      // Just an inequality, which we've already read.
      continue;
    }
    objcoef[objInd][colInd] = val;
    newRowMatVal[objInd * cur_numcols + colInd] = val;
  }

  // Now we need to set up the RHS.
  rhs = new double[objcnt];

  if (status) {
    std::cerr << "Failed to get RHS." << std::endl;
    return -ERR_CPLEX;
  }
  /* Get objective sense */
  int cpx_sense = CPXgetobjsen(e.env, e.lp);
  objsen = (cpx_sense == CPX_MIN ? MIN : MAX);
  /* Set objective constraint sense and RHS */
  consense = new char[objcnt];
  for (int j = 0; j < objcnt; j++) {
    if (objsen == MIN) {
      consense[j] = 'L'; /* Set sense to <= */
      rhs[j] = CPX_INFBOUND;
    }
    else {
      consense[j] = 'G'; /* Set sense to >= */
      rhs[j] = -CPX_INFBOUND;
    }
  }

  delete[] colNames;
  colNames = new char*[objNames.size()];
  for(size_t i = 0; i < objNames.size(); ++i) {
    colNames[i] = new char[objNames[i].size()+1];
    for(size_t j = 0; j < objNames[i].size(); j++) {
      colNames[i][j] = objNames[i][j];
    }
    colNames[i][objNames.size()] = '\0';
  }
  CPXaddrows(e.env, e.lp, 0 /*ccnt*/, objNames.size(), new_nzcnt, rhs, consense, newRowMatBeg, newRowMatInd, newRowMatVal, colNames, NULL /*rowname*/);

  conind = new int[objcnt];
  /* Specify index of objective constraints */
  for (int j = 0; j < objcnt; j++) {
    conind[j] = cur_numrows+j;
  }

  /* Set sense of objective constraints */
  status = CPXchgsense (e.env, e.lp, objcnt, conind, consense);
  if (status) {
    std::cerr << "Failed to change constraint sense" << std::endl;
    return -ERR_CPLEX;
  }

  delete[] newRowMatVal;
  delete[] newRowMatInd;
  delete[] newRowMatBeg;
  return 0;
}
