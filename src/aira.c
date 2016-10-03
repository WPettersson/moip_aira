#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include <libgen.h>
#include <time.h>
#include <glib.h>

/* List node. */
typedef struct {
  double *ip;
  int *result;
  int infeasible;
} Problem;

/* List of possible relaxations */
GSList* list = NULL;

/* CPLEX environment and problem */
CPXENVptr env = NULL;
CPXLPptr  lp = NULL;

/* Objective count */
int objcnt; 

/* Solve CLMOIP and return status */
int solve(int objcnt, int **objind, double **objcoef, int *conind, double *rhs, int *result, int *ipcount);

/* Add new solution/problem to list of possible relaxations */
void add_to_list(double *ip, int *result, int infeasible);

/* Search list for relaxations.
 * ip: Current problem
 * ipr: Place holder for the relaxation to the current problem
 */
int getrelaxation(double *ip, double *ipr, int *result, int *infeasible, int objsen);

/* Compare two solutions (from struct Problem) 
 * returns negative value if a < b; zero if a = b; positive value if a > b.
 * */
gint icmp(gconstpointer a, gconstpointer b);

/* Clean up Problem nodes in the list */
void  free_list(gpointer data);

int main (int argc, char *argv[])
{

  int status = 0; /* Operation status */
  int solnstat; /* Solution status (from CPLEX) */
  int i,j,k; /* counters */

  /* LP File Processing  */
  int from, to; /* First row, last row */
  int cur_numrows, cur_numcols, cur_numnz; /* # rows in current problem, # cols, # nonzero values*/
  int rmatspace, cur_objsen; /* Length of rmatind & rmatval, objective sense of current objective */
  int nzcnt, surplus; /* Nonzero row value count, array length check */
  int *rmatbeg, *rmatind; /* Row indices, column indices of rmatval */ 
  double *rhs, *ipr, *rmatval; /* Constraint RHS, relaxation ip, nozero row values */ 
  int **objind; /* Objective indices */
  double **objcoef; /* Objective coefficients */
  int *conind; /* Constraint indices */
  char *consense; /* Constraint senses */

  /* Algorithm */
  int infcnt; /* Infeasibility count */
  int *infeasible; /* Infeasibility flag for relaxation */
  int inflast, relaxed; /* Last iteration was feasible?, did we find a relaxation? */ 
  int *result, *max, *min; /* Current result, max result, min result */

  /* LP and log file */
  FILE *outfp;
  char *lpfn, *outfn, *tmpstr;

  /* Utility singly linked lists */
  GSList *iterator = NULL;
  GSList *printed = NULL; 
  GSList *tmplst = NULL;
  
  /* Timing */
  clock_t starttime, endtime;
  time_t startelapsed, endelapsed; 
  double cpu_time_used, elapsedtime; 
  int *ipcount;
  int solcount;

  if (argc < 2) {
    printf("\nUsage: aira <filename.lp>\n\n");
    exit(1);
  }
  
  /* Last arg is the lp file */
  lpfn = argv[argc-1];

  /* Output file: <filename>.out */ 
  tmpstr = (char *) malloc(strlen(lpfn) * sizeof(char));
  outfn = (char *) malloc(strlen(lpfn) * sizeof(char));
  tmpstr = strrchr(lpfn, '.');
  strncpy(outfn, lpfn, tmpstr-lpfn);
  outfn[tmpstr-lpfn] = '\0';
  strcat(outfn, ".out");
  outfp = fopen(outfn, "w+");

  /* Initialize the CPLEX environment */
  env = CPXopenCPLEX (&status);
  
  /* Set to deterministic parallel mode */
  status=CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);  

  if (env == NULL) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
    goto TERMINATE;
  }

 status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
 if (status) {
   fprintf (stderr, 
       "Failure to turn off screen indicator, error %d.\n", status);
   goto TERMINATE;
 }

  /* Create the problem, using the filename as the problem name */
  lp = CPXcreateprob(env, &status, lpfn);

  if (lp == NULL) {
    fprintf (stderr, "Failed to create LP.\n");
    goto TERMINATE;
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXreadcopyprob(env, lp, lpfn, NULL);
  if (status) {
    fprintf (stderr, "Failed to read and copy the problem data.\n");
    goto TERMINATE;
  }
  
  /* Get last rhs and determine the number of objectives.*/
  cur_numcols = CPXgetnumcols(env, lp);
  cur_numrows = CPXgetnumrows(env, lp);
  cur_numnz = CPXgetnumnz(env, lp);

  rhs = (double *) malloc(cur_numrows * sizeof(double));
  if (rhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Get RHS of last objective - this tells us the objective count */
  status = CPXgetrhs (env, lp, rhs, cur_numrows-1, cur_numrows-1);
  
  if (status) { 
    fprintf (stderr, "Failed to get RHS.\n"); 
    goto TERMINATE;
  }

  objcnt = (int) rhs[0];
  
  /* Create a pair of multidimensional arrays to store the objective 
   * coefficients and their indices indices first */
  objind =  malloc(objcnt * sizeof(int *));
  if (objind == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for(j=0; j<objcnt; j++) {
    objind[j] = (int *) malloc(cur_numcols * sizeof(int));
    if (objind[j] == NULL) {
      fprintf(stderr, "Error! memory is not available\n");
      goto TERMINATE;
    }
    for(i=0; i<cur_numcols;i++){
      objind[j][i] = i;
    }
  }
  /* Now coefficients */
  objcoef = malloc(objcnt * sizeof(double *));
  if (objcoef == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for(j=0; j<objcnt; j++) {
    objcoef[j] = (double *) malloc(cur_numcols * sizeof(double));
    if (objcoef[j] == NULL) {
      fprintf(stderr, "Error! memory is not available\n");
      goto TERMINATE;
    }
    memset(objcoef[j], 0, cur_numcols * sizeof(double));
  }

  /* Parse out the objectives working backwards from the last constraint */
  rmatbeg = (int *) malloc(cur_numrows * sizeof(int));
  rmatind = (int *) malloc(cur_numnz * sizeof(int));
  rmatval = (double *) malloc(cur_numnz * sizeof(double));
  if (rmatbeg == NULL || rmatind == NULL || rmatval == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  rmatspace = cur_numnz;

  status = CPXgetrows (env, lp, &nzcnt, rmatbeg, rmatind, rmatval,
                      rmatspace, &surplus, cur_numrows-objcnt, cur_numrows-1);

  if (status) {
    fprintf (stderr, "Couldn't get rows.\n");
    goto TERMINATE;
  }  

  for (j=0; j<objcnt; j++) {
    from = rmatbeg[j];
    if (j==objcnt-1) {
      to = nzcnt-1;
    }
    else {
      to = rmatbeg[(j+1)] - 1;
    }
    for (k=from;k<=to;k++){
      objcoef[j][rmatind[k]] = rmatval[k]; 
    }
  }

  /* Setup problem for solving */

  /* Resize rhs to fit only objective function constraints */
  free(rhs);
  rhs = (double *) malloc(objcnt * sizeof(double));
  if (rhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Get objective sense */
  cur_objsen = CPXgetobjsen(env, lp);

  /* Set objective constraint sense and RHS */
  consense = (char *) malloc(objcnt * sizeof(char));
  if (consense == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  for (j=0; j<objcnt; j++) {
    if (cur_objsen == CPX_MIN) {
      consense[j] = 'L'; /* Set sense to <= */
      rhs[j] = CPX_INFBOUND; 
    }
    else {
      consense[j] = 'G'; /* Set sense to >= */
      rhs[j] = -CPX_INFBOUND; 
    }
  }
 
  /* Specify index of objective constraints */
  conind = (int *) malloc(objcnt * sizeof(int)); 
  if (conind == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  k = objcnt;
  for (j=0; j<objcnt; j++) {
    conind[j] = cur_numrows-k;
    k--;
  }
  
  /* Set sense of objective constraints */
  status = CPXchgsense (env, lp, objcnt, conind, consense);
  if (status) {
     fprintf (stderr, "Failed to change constraint sense\n");
     goto TERMINATE;
  }

  /* Set rhs of objective constraints */
  status = CPXchgrhs (env, lp, objcnt, conind, rhs);
  if (status) {
     fprintf (stderr, "Failed to change constraint rhs\n");
     goto TERMINATE;
  }

  /* Allocate an array to store the result */
  result = (int *) malloc(objcnt * sizeof(int)); 
  if (result == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  
  /* Allocate an array to track the maximum values */
  max = (int *) malloc(objcnt * sizeof(int)); 
  if (max == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  /* Allocate an array to track the minimum values */
  min = (int *) malloc(objcnt * sizeof(int)); 
  if (min == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  ipcount = (int *) malloc(sizeof(int)); 
  if (ipcount == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }
  *ipcount = 0;

  fprintf(outfp, "\nUsing improved algorithm\n");
  ipr = (double *) malloc(objcnt * sizeof(double));
  if (ipr == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    goto TERMINATE;
  }

  /* Start the timer */
  starttime = clock();
  startelapsed = time(NULL);
  
  solnstat = solve(objcnt, objind, objcoef, conind, rhs, result, ipcount);
  memcpy(max, result, objcnt * sizeof(int));
  memcpy(min, result, objcnt * sizeof(int));

  infeasible = (int *) malloc(sizeof(int));
  *infeasible = FALSE;
  
  /* Need to add a result to the list here*/
  add_to_list(rhs, result, *infeasible);

  for (i=1; i<objcnt; i++) {
    int depth = 1; /* Track current "recursion" depth */
    int onwalk = FALSE; /* Are we on the move? */
    infcnt = 0; /* Infeasible count*/
    inflast = FALSE; /* Last iteration infeasible?*/
    
    /* Set all contraints back to infinity*/
    for (j=0; j<objcnt; j++) {
      if (cur_objsen == CPX_MIN) {
        rhs[j] = CPX_INFBOUND;
      }
      else {
        rhs[j] = -CPX_INFBOUND;
      }
    }

    /* Set rhs of current depth */
    if (cur_objsen == CPX_MIN) {
      rhs[i] = max[i]-1;
    }
    else {
      rhs[i] = min[i]+1;
    }
    max[i] = (int) -CPX_INFBOUND;
    min[i] = (int) CPX_INFBOUND;

    while (infcnt < i) {
      relaxed = FALSE;
    
      /* Look for possible relaxations to the current problem*/
      solnstat = getrelaxation(rhs, ipr, result, infeasible, cur_objsen);
      relaxed = solnstat ? TRUE : FALSE;
      
      if(!relaxed) {
        /* Solve in the absence of a relaxation*/
        solnstat = solve(objcnt, objind, objcoef, conind, rhs, result, ipcount);
      }
      
      if (solnstat == CPXMIP_INFEASIBLE || *infeasible) { 
        infcnt++;
        inflast = TRUE;
      }
      else {
        infcnt = 0;
        inflast = FALSE;
        /* Update maxima */
        for (j=1; j<objcnt; j++) {
          if (result[j] > max[j]) {
            max[j] = result[j];
          }
        }
        /* Update minima */
        for (j=1; j<objcnt; j++) {
          if (result[j] < min[j]) {
            min[j] = result[j];
          }
        }
      }
      /* Store result */
      if (!relaxed) {
        add_to_list(rhs, result, *infeasible);
      }

      if (infcnt == i-1) {
        /* Set all contraints back to infinity */
        for (j=0; j<objcnt; j++) {
          if (cur_objsen == CPX_MIN) {
            rhs[j] = CPX_INFBOUND;
          }
          else {
            rhs[j] = -CPX_INFBOUND;
          }
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set 
         * current level to min objective function value  +1 */
        if (cur_objsen == CPX_MIN) {
          rhs[i] = max[i]-1;
          max[i] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[i] = min[i]+1;
          min[i] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth = 1;
        onwalk = FALSE;
      }
      else if (inflast && infcnt != i) {
        rhs[depth] = CPX_INFBOUND;
        depth++;
        if (cur_objsen == CPX_MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = TRUE;
      }
      else if (!onwalk && infcnt != 1) {
        if (cur_objsen == CPX_MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      }
      else if (onwalk && infcnt != 1)  {
        depth = 1;
        if (cur_objsen == CPX_MIN) {
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        }
        else {
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = FALSE;
      }
    }
  } 

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  endelapsed = time (NULL);
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  elapsedtime=(double) endelapsed - startelapsed;

  list = g_slist_sort(list, (GCompareFunc)icmp);
  list = g_slist_reverse(list);
  solcount = 0;

  for (iterator=list; iterator; iterator=iterator->next) {
    if (printed != NULL) {
      tmplst = g_slist_find_custom (printed, (gconstpointer)iterator->data, (GCompareFunc)icmp);
    }
    if (tmplst == NULL) {
       for (i=0; i<objcnt; i++) {
          fprintf(outfp,"%d \t", ((Problem*)iterator->data)->result[i]);
        }
        fprintf(outfp,"\n");
      printed = g_slist_append(printed,((Problem*)iterator->data));
      solcount++;
    }
  }

  fprintf(outfp,"\n---\n%8.4f CPU seconds\n", cpu_time_used);
  fprintf(outfp,"%8.4f elapsed seconds\n", elapsedtime);
  fprintf(outfp,"%8d IPs solved\n", *ipcount);
  fprintf(outfp,"%8d Solutions found\n", solcount);

  
  TERMINATE:
   /* Free up memory as necessary. */
   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }
   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }
   if (rhs != NULL) free(rhs);
   if (rmatbeg != NULL) free(rmatbeg);
   if (rmatind != NULL) free(rmatind);
   if (rmatval != NULL) free(rmatval);
   if (objind !=NULL) {
    for(j=0; j<objcnt; j++) {
      if (objind[j] != NULL) free(objind[j]);
    }
   }
   if (objind != NULL) free(objind);
   if (objcoef !=NULL) {
    for(j=0; j<objcnt; j++) {
      if (objcoef[j] != NULL) free(objcoef[j]);
    }
   }
   if (objcoef != NULL) free(objcoef);
   if (conind != NULL) free(conind);
   if (consense != NULL) free(consense);
   if (outfp != NULL) fclose(outfp);
   if (printed != NULL) g_slist_free(printed);
   if (list != NULL) g_slist_free_full(list, free_list);
     
   return (status);
}

/* Solve CLMOIP and return solution status */
int solve(int objcnt, int **objind, double **objcoef, int *conind, double *rhs,
    int *result, int *ipcount) {
  
  int j;
  int cur_numcols, status, solnstat;
  double objval, *srhs;
  srhs = (double *) malloc(objcnt * sizeof(double));
  if (srhs == NULL) {
    fprintf(stderr, "Error! memory is not available\n");
    return -1;
  }

  memcpy(srhs, rhs, objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(env, lp);
  
  for (j=0;j<objcnt;j++) {
    status = CPXchgobj(env, lp, cur_numcols, objind[j], objcoef[j]);
    if (status) {
       fprintf (stderr, "Failed to set objective.\n");
    }

    status = CPXchgrhs (env, lp, objcnt, conind, srhs);
    if (status) {
      fprintf (stderr, "Failed to change constraint srhs\n");
    }
  
    /* solve for current objective*/
    status = CPXmipopt (env, lp);
    if (status) {
       fprintf (stderr, "Failed to optimize LP.\n");
    }

    (*ipcount)++;
    
    solnstat = CPXgetstat (env, lp);
    if (solnstat == CPXMIP_INFEASIBLE) {
       break;
    }
    status = CPXgetobjval (env, lp, &objval);
    if ( status ) {
       fprintf (stderr, "Failed to obtain objective value.\n");
    }
    
    result[j] = round(objval);
    srhs[j] = round(objval);
  }

  return solnstat;
}

/* Add new solution/problem to list of possible relaxations */
void add_to_list(double *ip, int *result, int infeasible) {
  Problem* p = (Problem*) malloc(sizeof(Problem)); 
  p->ip = malloc(objcnt * sizeof(double));
  memcpy(p->ip, ip, objcnt * sizeof(double));
  if (infeasible) {
    p->infeasible = TRUE;
    p->result = NULL;
  }
  else {
    p->infeasible = FALSE;
    p->result = malloc(objcnt * sizeof(int));
    memcpy(p->result, result, objcnt * sizeof(int));
  }
  list = g_slist_append(list, p);
}

/* Search list for relaxations.
 * ip: Current problem
 * ipr: Place holder for the relaxation to the current problem
 */
int getrelaxation(double *ip, double *ipr, int *result, int *infeasible, int objsen) {

  GSList *iterator;
  int t1,t2,t3,i;
  for (iterator=list; iterator; iterator=iterator->next) {

    t1 = TRUE;
    t2 = FALSE;
    t3 = TRUE;

    if(objsen == CPX_MIN) {
      for (i=0; i<objcnt; i++){
        /* t1: All values of candidate must be >= than ip */
        if (((Problem*)iterator->data)->ip[i] < ip[i]) {
          t1 = FALSE;
          break;
        }
        /* t2: At least one inequality must be strict */
        if (((Problem*)iterator->data)->ip[i] > ip[i]) {
          t2 = TRUE;
        }
        /* t3: All values of candidate result must be <= than ip */
        if (!((Problem*)iterator->data)->infeasible) {
          if (((Problem*)iterator->data)->result[i] > ip[i]) {
            t3 = FALSE;
            break;
          }
        }
      }
    }
    else {
      for (i=0; i<objcnt; i++){
        /* t1: All values of candidate must be <= than ip */
        if (((Problem*)iterator->data)->ip[i] > ip[i]) {
          t1 = FALSE;
          break;
        }
        /* t2: At least one inequality must be strict */
        if (((Problem*)iterator->data)->ip[i] < ip[i]) {
          t2 = TRUE;
        }
        /* t3: All values of candidate result must be >= than ip */
        if (!((Problem*)iterator->data)->infeasible) {
          if (((Problem*)iterator->data)->result[i] < ip[i]) {
            t3 = FALSE;
            break;
          }
        }
      }
    }
    /* If all conditions are met copy problem & solution and return */
    if (t1 && t2 && t3) {
      memcpy(ipr, ((Problem*)iterator->data)->ip, objcnt * sizeof(double));
      memcpy(infeasible, &((Problem*)iterator->data)->infeasible, sizeof(int));
      if (!*infeasible) {
        memcpy(result, ((Problem*)iterator->data)->result, objcnt * sizeof(int));
      }
      return 1;
    }
  }
  return 0;
}

/* Compare two solutions (from struct Problem) 
 * returns negative value if a < b; zero if a = b; positive value if a > b. */
gint icmp(gconstpointer a, gconstpointer b) {
  int i;
  for (i=0;i<objcnt;i++) {
    if (((Problem*)a)->result[i] < ((Problem*)b)->result[i]) { 
      return (gint)-1;
    }
    else if (((Problem*)a)->result[i] > ((Problem*)b)->result[i]){ 
      return (gint)1;
    }
  }
  return (gint)0;
}

/* Clean up Problem nodes in the list */
void free_list(gpointer data) {
  if (((Problem*)data)->ip != NULL) {
    free(((Problem*)data)->ip);
  }
  if (((Problem*)data)->result != NULL) {
    free(((Problem*)data)->result);
  }
} 

