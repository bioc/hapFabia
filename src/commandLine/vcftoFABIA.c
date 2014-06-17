#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 
int vcftoFABIAB(const char *arg1, const char *arg2, const char *arg3, const char *arg4);

void Rprintf(const char * RSformat, ...) {

  printf("%s",RSformat);

  return;

}

int main(int argc, char* argv[]) {

  char agr1[500],agr2[500],agr3[500],agr4[500];
  int narg,i;

  narg=argc;
  agr1[0]=0;
  strcat(agr1,argv[1]);
  agr2[0]=0;
  strcat(agr2,argv[2]);
  agr3[0]=0;
  strcat(agr3,"NA");
  agr4[0]=0;
  strcat(agr4,"NA");
  if (narg>2) {
    i=narg;
    for (i = 3; i < argc; i++)
    {
        char const *option =  argv[i];
        if (option[0] == '-')
        {
            switch (option[1])
            {
                case 's':
		  i++;
		  agr3[0]=0;
		  strcat(agr3,argv[i]);
		  break;
                case 'o':
		  i++;
		  agr4[0]=0;
		  strcat(agr4,argv[i]);
		  break;
                default:
                    printf("flag %s not recognised\n", option);
                    break;
            }
        }
        else
        {   
	  printf("flag %s not recognised\n", option);
	}



    }

  }

  vcftoFABIAB(agr1,agr2,agr3,agr4);

  return 0;
}



