#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int split_sparse_matrixB(int narg, const char *arg1, const char *arg2, const char *arg3, const char *arg4, const char *arg5);
 
void Rprintf(const char * RSformat, ...) {

  printf(RSformat);

  return;

}

int main(int argc, char* argv[]) {

  char agr1[500],agr2[500],agr3[500],agr4[500],agr5[500]; 
  int narg;

    narg=argc;
    agr1[0]=0;
    agr2[0]=0;
    agr3[0]=0;
    agr4[0]=0;
    agr5[0]=0;
    strcat(agr1,argv[1]);
    strcat(agr2,argv[2]);
    if (argc>3) {
      strcat(agr3,argv[3]);
    }
    if (argc>4) {
      strcat(agr4,argv[4]);
    }
     if (argc>5) {
      strcat(agr5,argv[5]);
    }


  split_sparse_matrixB(narg,agr1,agr2,agr3,agr4,agr5);

  return 0;
}


