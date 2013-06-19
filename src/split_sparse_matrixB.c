#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void Rprintf(const char *, ...);


char *sput_i(int integer, char *string)
{
  if (integer / 10 != 0) {
    string = sput_i(integer / 10, string);
  }
  *string++ = (char)('0' + integer % 10);
  return string;
}

char *sput_ip1(int integer, char *string)
{
  int digit;

  digit = (integer % 10 + 1) % 10;
  if (integer / 10 != 0) {
    string = (digit == 0 ? sput_ip1 : sput_i)(integer / 10, string);
    *string++ = (char)('0' + digit);
  } else {
    if (digit == 0) {
      *string++ = '1';
    }
    *string++ = (char)('0' + digit);
  }
  return string;
}


void itoaS(int integer, char *string)
{
  if (0 > integer) {
    ++integer;
    *string++ = '-';
    *sput_ip1(-integer, string) = '\0';
  } else {
    *sput_i(integer, string) = '\0';
  }
}



int split_sparse_matrixB(int narg, char *agr1,  char *agr2,  char *agr3 , char *agr4,  char *agr5 ) {

  FILE *pFile,*pFile1;
    
  int i,i1,j,j1,individuals,snps,size,shift,ig,parts,
    annotation,inI,start,end,diff,rest;

    //snps =  3530000;
    
    //individuals = 1092;

    //individuals = 1094;
    //snps = 1844361;


    size = atoi(agr3);
    shift = atoi(agr4);
    annotation = atoi(agr5);


    if (shift>size) {
      Rprintf("shift %d is larger than size %d, therefore shift = %d!\n",shift,size,size);
      shift=size;
    }
 
    char sst[500],str[500];

    char snpName[500], major[2000],minor[2000],filter[200],
      info[2000],format[200],qual[200];
    int chrom,pos,change;
    double freq;





    sst[0]=0;
    strcat(sst,agr1);
    strcat(sst,agr2);
    pFile = fopen(sst,"r");

    if (pFile==NULL) {
      Rprintf("File >%s< not found! Stop.\n", sst);
      return(-1);
    }

    fscanf(pFile,"%d\n",&individuals);  
    fscanf(pFile,"%d\n",&snps);  
 
    parts = (snps-size)/shift+2; 
    rest = snps-(parts-1)*shift;

    char **fileN = calloc(parts, sizeof(char *));
    fileN[0] =  calloc((long) 200*parts, sizeof(char));
    for(i=0; i < parts; i++)
    {
	fileN[i] =  fileN[0]+i*200;
    }

    int *startP = calloc(parts , sizeof(int));
    int *endP = calloc(parts , sizeof(int));

    for(i = 0; i < parts; i ++)
      {

	sst[0]=0;
	strcat(sst,agr1);
	strcat(sst,"_");
	str[0]=0;
	itoaS((i*shift),str);
	strcat(sst,str);
	strcat(sst,"_");
	if (i<(parts-1)) {
	  str[0]=0;
	  itoaS((i*shift+size),str);
	  strcat(sst,str);
	} else {
	  str[0]=0;
	  itoaS((i*shift+rest),str);
	  strcat(sst,str);
	}
	strcat(sst,agr2);
	pFile1 = fopen(sst,"w");
	if (pFile1==NULL) {
	  Rprintf("File >%s< cannot be opened! Stop.\n", sst);
	  return(-1);
	}
	fprintf(pFile1,"%d\n",individuals);  
	if (i<(parts-1)) {
	  fprintf(pFile1,"%d\n",size);  
	} else {
	  fprintf(pFile1,"%d\n",rest);  
	} 
	fclose (pFile1);

	strcat(fileN[i],sst);
	
	
      }


    double *Lval = calloc((long) snps , sizeof(double));
    
    int *Lind = calloc((long) snps , sizeof(int));



    for(i = 0; i < individuals; i ++)
    {

      Rprintf("Write Individual: %d\r",(i+1));

	fscanf(pFile,"%d\n",&ig);
	for(j = 0; j < ig; j ++)
	    fscanf(pFile,"%d ",&Lind[j]);
	fscanf(pFile,"\n");
	for(j = 0; j < ig; j ++)
	    fscanf(pFile,"%lf ",&Lval[j]);
	fscanf(pFile,"\n");

	for(i1 = 0; i1 < parts; i1 ++)
	{
	  startP[i1]=-1;
	}
	i1=0;
	j1=0;
	for(j = 0; j < ig; j ++)
	  {
	    inI= Lind[j];
	    while (inI>=(i1*shift)) {
	      if ((inI<(i1*shift+size))&&(startP[i1]==-1)) {
		  startP[i1]=j;
		}
	      i1++;
	    }
	    if (j==(ig-1)) {
	      for(; j1 < parts; j1++) {
		endP[j1]=j;
	      }
	    } else {
	      while (inI>=(j1*shift+size)) {
		endP[j1]=j;
		j1++;
	      }
	    }
	  }

	endP[parts-1]++;

	for(i1 = 0; i1 < parts; i1 ++)
	{
	    
	  pFile1 = fopen(fileN[i1],"a");
	  if (pFile1==NULL) {
	    Rprintf("File >%s< cannot be opened! Stop.\n", fileN[i1]);
	    return(-1);
	  }
	  

	  if (startP[i1]==-1)
	  {
	    fprintf(pFile1,"0\n");  
	    fprintf(pFile1,"\n");  
	    fprintf(pFile1,"\n");  
	  } else {
	    
	    start=startP[i1];
	    end=endP[i1];
	    
	    diff=end-start;	
	    fprintf(pFile1,"%d\n",diff);  
	    
	    for (j1 = start; j1<end; j1++) {
	      j=Lind[j1]-i1*shift;
	      if (i1<(parts-1)) {
		if ((j>=0)&&(j<size))
		  {
		    fprintf(pFile1,"%d ",j); 
		  } else {
		  Rprintf("index of %d larger than end of %d. Stop!\n",Lind[j1],(i1*shift+size));
		  return(-1);
		}
	      } else {
		if ((j>=0)&&(j<rest))
		  {
		    fprintf(pFile1,"%d ",j); 
		  } else {
		  Rprintf("index of %d larger than end of %d. Stop!\n",Lind[j1],(i1*shift+rest));

		  return(-1);
		}

	      } 
	    }
	    fprintf(pFile1,"\n");  
	    for (j1 = start; j1<end; j1++) {
	      fprintf(pFile1,"%.3f ",Lval[j1]); 
	    }
	    fprintf(pFile1,"\n");  
	  }
	  
	  
	  fclose (pFile1);
	  
	}
	
    }
    fclose (pFile);
  



    Rprintf("\n");    



    if (annotation>0) {

    Rprintf("Write Annotation.\n");    
 
    sst[0]=0;
    strcat(sst,agr1);
    strcat(sst,"_annot.txt");
    pFile = fopen(sst,"r");
    if (pFile==NULL) {
      Rprintf("File >%s< not found! Stop.\n", sst);
      return(-1);
    }
 

    fscanf(pFile,"Individuals: %d\n",&individuals);  
    fscanf(pFile,"SNVs       : %d\n",&snps);  
 


    char **fileN1 = calloc(parts, sizeof(char *));
    fileN1[0] =  calloc((long) 200*parts, sizeof(char));
    for(i=0; i < parts; i++)
    {
	fileN1[i] =  fileN1[0]+i*200;
    }


   for(i = 0; i < parts; i ++)
      {

	sst[0]=0;
	strcat(sst,agr1);
	strcat(sst,"_");
	str[0]=0;
	itoaS((i*shift),str);
	strcat(sst,str);
	strcat(sst,"_");
	if (i<(parts-1)) {
	  str[0]=0;
	  itoaS((i*shift+size),str);
	  strcat(sst,str);
	} else {
	  str[0]=0;
	  itoaS((i*shift+rest),str);
	  strcat(sst,str);
	}
	strcat(sst,"_annot.txt");
	pFile1 = fopen(sst,"w");
	if (pFile1==NULL) {
	  Rprintf("File >%s< cannot be opened! Stop.\n", sst);
	  return(-1);
	}
	fprintf(pFile1,"%d\n",individuals);  
	if (i<(parts-1)) {
	  fprintf(pFile1,"%d\n",size);  
	} else {
	  fprintf(pFile1,"%d\n",rest);  
	} 
	fclose (pFile1);

	strcat(fileN1[i],sst);
	
	
      }


   i1=0;
   j1=1;
   for(j = 0; j < snps; j ++)
   {
     fscanf(pFile,"%d %d %s %s %s %s %s %s %s %lf %d\n",&chrom,&pos, snpName, major, minor,qual,filter,info,format,&freq,&change);
     
     

     if (j>=i1*shift+size) {
       i1++;
     }
     if (j>=j1*shift) {
       if (j1<parts) j1++;
     }


     for (i=i1; i<j1; i++) { 
       pFile1 = fopen(fileN1[i],"a");
       if (pFile1==NULL) {
	 Rprintf("File >%s< cannot be opened! Stop.\n", fileN1[i]);
	 return(-1);
       }
       fprintf(pFile1,"%d %d %s %s %s %s %s %s %s %.8f %d\n",chrom,pos, snpName, major, minor,qual,filter,info,format,freq,change);
       fclose (pFile1);
     }
     }   
   fclose(pFile);
  

   
    }
     
    Rprintf("\n");    
    Rprintf("\n");    

    return(0);
    
}

