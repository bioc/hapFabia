#include <stdlib.h>
#include <stdio.h>  
#include <limits.h> 
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

int split_sparse_matrixB(int narg, const char *arg1, const char *arg2, const char *arg3, const char *arg4, const char *arg5);
 
int vcftoFABIAB(const char *arg1, const char *arg2, const char *arg3, const char *arg4);

void split_sparse_matrix(SEXP nargS,SEXP arg1S,SEXP arg2S,SEXP arg3S,SEXP arg4S,SEXP arg5S) {

  
  int narg = (int)(INTEGER(nargS)[0]);
  const char *arg1=CHAR(STRING_ELT(arg1S,0));
  const char *arg2=CHAR(STRING_ELT(arg2S,0));
  const char *arg3=CHAR(STRING_ELT(arg3S,0));
  const char *arg4=CHAR(STRING_ELT(arg4S,0));
  const char *arg5=CHAR(STRING_ELT(arg5S,0));
  
  split_sparse_matrixB(narg,arg1,arg2,arg3,arg4,arg5);

  return;

}

void vcftoFABIA(SEXP arg1S,SEXP arg2S,SEXP arg3S,SEXP arg4S) {

  const char *arg1=CHAR(STRING_ELT(arg1S,0));
  const char *arg2=CHAR(STRING_ELT(arg2S,0));
  const char *arg3=CHAR(STRING_ELT(arg3S,0));
  const char *arg4=CHAR(STRING_ELT(arg4S,0));
  
  vcftoFABIAB(arg1,arg2,arg3,arg4);

  return;
}




 R_CallMethodDef callMethods[]  = {
       {"split_sparse_matrix", (DL_FUNC) &split_sparse_matrix, 6},
       {"vcftoFABIA", (DL_FUNC) &vcftoFABIA, 4},
       {NULL, NULL, 0}
     };



 void  R_init_myLib(DllInfo *info)
     {
	 R_registerRoutines(info, NULL, callMethods, NULL, NULL);
     }


int main() {

  return(1);
}


