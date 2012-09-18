#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 
int vcftoFABIAB(int narg, const char *arg1, const char *arg2, const char *arg3);

void Rprintf(const char * RSformat, ...) {

  printf(RSformat);

  return;

}

int main(int argc, char* argv[]) {

  char agr1[500],agr2[500],agr3[500];
  int narg;

  narg=argc;
  agr1[0]=0;
  agr2[0]=0;
  agr3[0]=0;
  strcat(agr1,argv[1]);
  strcat(agr2,argv[2]);
  if (argc>3) {
    strcat(agr3,argv[3]);
  }

  vcftoFABIAB(narg,agr1,agr2,agr3);

  return 0;
}



