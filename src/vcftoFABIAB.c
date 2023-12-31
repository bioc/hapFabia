#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *strtok_rS(char *s, const char *delim, char **save_ptr)
{
    char *token;

    if (s == NULL)
        s = *save_ptr;

    /* Scan leading delimiters. */
    s += strspn (s, delim);
    if (*s == '\0')
        return NULL;

    /* Find the end of the token. */
    token = s;
    s = strpbrk (token, delim);
    if (s == NULL)
        /* This token finishes the string. */
        *save_ptr = strchr (token, '\0');
    else
    {
        /* Terminate the token and make *SAVE_PTR point past it. */
        *s = '\0';
        *save_ptr = s + 1;
    }
    return token;
}


ssize_t getdelimS(char **linebuf, size_t *linebufsz, int delimiter, FILE *file)
{
   static const int GROWBY = 4096;      /* grow strings by */
   int ch;
   size_t idx = 0;
   if (file == NULL || linebuf == NULL || linebufsz == NULL) {
      return -1;
   }
   if (*linebuf == NULL || *linebufsz < 2) {
      *linebuf = malloc(GROWBY);
      if (!*linebuf) {
         return -1;
      }
      *linebufsz += GROWBY;
   }
   while (1) {
      ch = fgetc(file);
      if (ch == EOF)
         break;
     
         /* grow the line buffer as necessary */ 
         while (idx > *linebufsz - 2) {
         *linebuf = realloc(*linebuf, *linebufsz += GROWBY);
         if (!*linebuf) {
            return -1;
         }
      }
      (*linebuf)[idx++] = (char) ch;
      if ((char) ch == delimiter)
         break;
   }
   if (idx != 0)
      (*linebuf)[idx] = 0;
   else if (ch == EOF)
      return -1;
   return idx;
}


ssize_t getlineS(char **linebuf, size_t *n, FILE *file) 
{
   return (getdelimS(linebuf, n, '\n', file));
}



void Rprintf(const char *, ...);

int vcftoFABIAB(const char *agr1, const char *agr2, const char *agr3, const char *agr4) {

  FILE *pFile, *pFile1;
  
  int i=0,j=0,header_lines=0,individuals=0,ig=0,haplo=0,dsFound=0,Csnps=0;
  
  long snps=0;
  int lineC=0;
  int countL=0,samples=0,posG=0,GTpos=0,DSpos=0,bo=0,ret=0;
  char *p,*p1,*saveptr1=NULL,*saveptr2=NULL;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  
  float eps,s;
  
  
  unsigned int hpp;
  
  
  unsigned short gh1,gh2;
  float dsI;
  
  
  char sst[500],outS[500],IsnpName[5000], Imajor[2000000],Iminor[2000000],Ifilter[2000],Iinfo[20000],Iformat[2000],Iqual[2000];
  
  unsigned short **g1= NULL;
  unsigned short **g2= NULL;
  float **ds= NULL;
  unsigned short *chrom= NULL;
  unsigned int *pos= NULL;
  float *freq= NULL;
  float *Lval= NULL;
  unsigned int *Lind= NULL;
  unsigned short *change= NULL;
  char **major= NULL;
  char **minor= NULL;
  char **snpName= NULL;
  char **filter= NULL;
  char **qual= NULL;
  char **info= NULL;
  char **format= NULL;
      


  eps=0.6;
 
  sst[0]=0;
  strcat(sst,agr2);
  strcat(sst,agr1);
  strcat(sst,".vcf");

  pFile = fopen (sst,"r");

  if (pFile==NULL) {
    Rprintf("File >%s< not found! Stop.\n", sst);
    return(-1);
  }

  if (strcmp(agr3,"NA")==0) {
    lineC=0;
    while ((getlineS(&line, &len, pFile)) != -1)   
      ++lineC;
    rewind(pFile);
    Csnps=1;
  } else {
    snps = atoi(agr3); 
  }
  
  
  if (strcmp(agr4,"NA")==0) {
    outS[0]=0;
    strcat(outS,agr1);
  } else {
    outS[0]=0;
    strcat(outS,agr4);
  }
 
  countL=0;
  j=0;
  while ((read = getlineS(&line, &len, pFile)) != -1) {
    j++;
    if ( line[0]=='#' ) {
      if ( line[1]=='C' ) {
	
	if (Csnps==1) {
	  snps = lineC-j;
	} 
	
	header_lines = j;

	Rprintf ("SNVs: %d\n", snps);
	
	snps+=10;
	
	j=-1;

	samples=0;
	p = strtok_rS (line,"\t",&saveptr1);
	
	sst[0]=0;
	strcat(sst,outS);
	strcat(sst,"_individuals.txt");
	pFile1 = fopen (sst,"w");
	
	while (p != NULL)
	  {
	    samples++;
	    if (samples>9) {
	      fprintf(pFile1,"%d %s\n", (samples-9),p);
	    }
	    p = strtok_rS (NULL, "\t",&saveptr1);
	  }
	fclose (pFile1);
	individuals = samples-9;
	Rprintf ("Individuals: %d\n", individuals);
      }
      
    } else {
      
      //##################################################################
      if (j==0) { 
      g1 = calloc(individuals, sizeof(unsigned short *));
      g1[0] =  calloc((long) individuals*snps, sizeof(unsigned short));
      for(i=0; i < individuals; i++)
	{
	  g1[i] =  g1[0]+i*snps;
	}
      g2 = calloc(individuals, sizeof(unsigned short *));
      g2[0] =  calloc((long) individuals*snps, sizeof(unsigned short));
      for(i=0; i < individuals; i++)
	{
	  g2[i] =  g2[0]+i*snps;
	}
      
      ds = calloc(individuals, sizeof(float *));
      ds[0] =  calloc((long) individuals*snps, sizeof(float));
      for(i=0; i < individuals; i++)
	{
	  ds[i] =  ds[0]+i*snps;
	}
      
      
     
      chrom = calloc((long) snps , sizeof(unsigned short));
      pos = calloc((long) snps , sizeof(unsigned int));
      freq = calloc((long) snps , sizeof(float));
      Lval = calloc((long) snps , sizeof(float));
      Lind = calloc((long) snps , sizeof(unsigned int));
      change = calloc((long) snps , sizeof(unsigned short));
      
      // char *Imajor = calloc(2000000, sizeof(char));
      // char *Iminor = calloc(2000000, sizeof(char));

      major = calloc(snps, sizeof(char *));
      major[0] =  calloc((long) 200*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  major[i] =  major[0]+i*200;
	}
      
     minor = calloc(snps, sizeof(char *));
      minor[0] =  calloc((long) 200*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  minor[i] =  minor[0]+i*200;
	}
      
      snpName = calloc(snps, sizeof(char *));
      snpName[0] =  calloc((long) 50*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  snpName[i] =  snpName[0]+i*50;
	}
      
      
      filter = calloc(snps, sizeof(char *));
      filter[0] =  calloc((long) 20*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  filter[i] =  filter[0]+i*20;
	}
      
      
      qual = calloc(snps, sizeof(char *));
      qual[0] =  calloc((long) 20*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  qual[i] =  qual[0]+i*20;
	}
      
      info = calloc(snps, sizeof(char *));
      info[0] =  calloc((long) 200*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  info[i] =  info[0]+i*200;
	}
      
      
      format = calloc(snps, sizeof(char *));
      format[0] =  calloc((long) 20*snps, sizeof(char));
      for(i=0; i < snps; i++)
	{
	  format[i] =  format[0]+i*20;
	}
      


      }
      //##################################################################

	if (j%1000==0) {
	  Rprintf("Read SNV: %d\r",j);
	}


	  i=0;

	  samples=0;
	  p = strtok_rS (line,"\t",&saveptr1);
  

	  while (p != NULL)
	    {

	      switch(samples) {
	      case 0:  
		ret=sscanf(p,"%hu", &gh1);
		if (ret<1) {
		  chrom[j] = 1;
		} else {
		chrom[j] = gh1;
		}
		break;

	      
	      case 1:  
		ret=sscanf(p,"%u", &hpp);
		if (ret<1) {
		  pos[j] = 0;
		} else {
		  pos[j] = hpp;
		}

		break;

	      case 2:  
		ret=sscanf(p,"%s", IsnpName);
		if (ret<1) {
		  strncat(snpName[j],"",49);
		} else {
		  strncat(snpName[j],IsnpName,49);
		}
		
		break;

	      case 3:  
		ret=sscanf(p,"%s", Imajor);
		if (ret<1) {
		  strncat(major[j],"",190);
		} else {
		  strncat(major[j],Imajor,190);
		}

		break;

	      case 4:  
		ret=sscanf(p,"%s",  Iminor);
		if (ret<1) {
		  strncat(minor[j],"",190);
		} else {
		  strncat(minor[j],Iminor,190);
		}

		break;

	      case 5:  
		ret=sscanf(p,"%s", Iqual);
		if (ret<1) {
		  strncat(qual[j],"",19);
		} else {
		  strncat(qual[j],Iqual,19);
		}

		break;

	      case 6:  
		ret=sscanf(p,"%s", Ifilter);
		if (ret<1) {
		  strncat(filter[j],"",19);
		} else {
		  strncat(filter[j],Ifilter,19);
		}

		break;

	      case 7:  
		ret=sscanf(p,"%s", Iinfo);
		if (ret<1) {
		  strncat(info[j],"",190);
		} else {
		  strncat(info[j],Iinfo,190);
		}

		break;



	    case 8:
		ret=sscanf(p,"%s",Iformat);
		if (ret<1) {
		  strncat(format[j],"",19);
		} else {
		  strncat(format[j],Iformat,19);
		}
		p1 =  strtok_rS (p,":",&saveptr2);
		GTpos=0;
		DSpos=0;
		posG=0;
		while (p1 != NULL)
		  {
		    posG++;
		    if (strcmp (p1,"GT") == 0) {
		      GTpos=posG;
		    }
		    if (strcmp (p1,"DS") == 0) {
		      DSpos=posG;
		    }
		    
		    p1 = strtok_rS (NULL, ":",&saveptr2);
		  }


		if (GTpos==0) {
		  Rprintf("Error in SNV %d and Line %d: no GT found!\n",(countL- header_lines), countL);
		  return(-1);
		}
		break;


	    default:	      
		p1 =  strtok_rS (p,":",&saveptr2);
		gh1=0;
		gh2=0;
		dsI=0.0;
		posG=0;

		while (p1 != NULL)
		  {
		    posG++;
		    if (posG==GTpos) {
		      //phased
		      bo=sscanf(p1,"%hu|%hu", &gh1, &gh2);
		      //unphased
		      if (bo<2) {
			bo=sscanf(p1,"%hu/%hu", &gh1, &gh2);
		      }
		      //haplotypes or pure genotype data
		      if (bo<2) {
			bo=sscanf(p1,"%hu", &gh1);
			gh2=0;
			haplo=1;    
		      }
		    } 

		    if (posG==DSpos) {
		      bo=sscanf(p1,"%f", &dsI);
		      if (bo>0) {
			dsFound=1;
		      }
		    }
		    
		    p1 = strtok_rS (NULL, ":",&saveptr2);
		  }

		g1[i][j] = gh1;
		g2[i][j] = gh2;
		ds[i][j] = dsI;

		i++;

	      break;

	    }


	      samples++;

	      p = strtok_rS (NULL, "\t",&saveptr1);

	    }

    }
    countL++;
  }
  fclose (pFile);

  snps= j+1;
 
  Rprintf ("\n Reconfirmed SNVs: %d\n", snps);



    // check for minor allele and change
    for(j = 0; j < snps; j ++)
      {

	for(ig=0,i = 0; i < individuals; i ++)
	  {
	    if (g1[i][j]>0)
	      {
		ig++;
	      }
	    if (haplo<1){
	    if (g2[i][j]>0)
	      {
		ig++;
	      }
	    }
	  }

	if (haplo==1){
	  freq[j] = (float) (1.0*ig)/(1.0*individuals);
	  change[j]=0;
	  if (ig > ( (int) individuals/2) )
	    {
	      change[j]=1;
	      for(i = 0; i < individuals; i ++)
		{
		  g1[i][j]= 1 - g1[i][j];
		  g2[i][j]= 0;
		  ds[i][j] = 1.0 - ds[i][j];
		}
	    }
	} else {
	  freq[j] = (float) (1.0*ig)/(2.0*individuals);
	  change[j]=0;
	  if (ig>individuals)
	    {
	      change[j]=1;
	      for(i = 0; i < individuals; i ++)
		{
		  g1[i][j]= 1 - g1[i][j];
		  g2[i][j]= 1 - g2[i][j];
		  ds[i][j] = 2.0 - ds[i][j];
		}
	    }
	}


      }


    //check for high probable snps
 /* for(j = 0; j < snps; j ++) */
 /*       { */

 /* 	for(i = 0; i < individuals; i ++) { */

 /* 	  if (p1[i][j]<0.8) { */
 /* 	      p1[i][j]=0.0; */
 /* 	  } */
 /* 	  if (p2[i][j]<0.8) { */
 /* 	      p2[i][j]=0.0; */
 /* 	  } */

 /* 	} */

 /*      } */









    sst[0]=0;
    strcat(sst,outS);
    strcat(sst,"_annot.txt");
    pFile = fopen (sst,"w");
    if (pFile==NULL) {
      Rprintf("File >%s< cannot be opened! Stop.\n", sst);
      return(-1);
    }
    fprintf(pFile,"Individuals: %d\n",individuals);  
    fprintf(pFile,"SNVs       : %d\n",(int) snps);  
    for(j = 0; j < snps; j ++)
    {
      fprintf(pFile,"%d %d %s %s %s %s %s %s %s %.8f %d\n",chrom[j],pos[j], snpName[j], major[j], minor[j],qual[j],filter[j],info[j],format[j],freq[j],change[j]);
    }
    fclose (pFile);

    Rprintf("\n");    

    sst[0]=0;
    strcat(sst,outS);
    strcat(sst,"_matH.txt");
    pFile = fopen (sst,"w");
    if (pFile==NULL) {
      Rprintf("File >%s< cannot be opened! Stop.\n", sst);
      return(-1);
    }
    if (haplo==1) {
      fprintf(pFile,"%d\n",individuals);  
    } else {
      fprintf(pFile,"%d\n",2*individuals);  
    }
    fprintf(pFile,"%d\n",(int) snps);  
    for(i = 0; i < individuals; i ++)
    {
      Rprintf("Write Individual Haplotype: %d\r",(i+1));
	ig=0;
	for(j = 0; j < snps; j ++)
	{
	    if (g1[i][j]>0){
		Lind[ig] = j;
		Lval[ig] = 1.0;
		ig++;
	    }
	    
	}
	fprintf(pFile,"%d\n",ig); 
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%d ",Lind[j]);
	fprintf(pFile,"\n");
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%.1f ",Lval[j]);
	fprintf(pFile,"\n");

        if (haplo<1) {
	ig=0;
	for(j = 0; j < snps; j ++)
	{
	    if (g2[i][j]>0){
		Lind[ig] = j;
		Lval[ig] = 1.0;
		ig++;
	    }
	    
	}
	fprintf(pFile,"%d\n",ig); 
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%d ",Lind[j]);
	fprintf(pFile,"\n");
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%.1f ",Lval[j]);
	fprintf(pFile,"\n");
	}
    }
    fclose (pFile);

    Rprintf("\n");    

    sst[0]=0;
    strcat(sst,outS);
    strcat(sst,"_matG.txt");
    pFile = fopen (sst,"w");
    if (pFile==NULL) {
      Rprintf("File >%s< cannot be opened! Stop.\n", sst);
      return(-1);
    }
    fprintf(pFile,"%d\n",individuals);  
    fprintf(pFile,"%d\n",(int) snps);  
    for(i = 0; i < individuals; i ++)
    {
      Rprintf("Write Individual Genotype: %d\r",(i+1));
	ig=0;
	for(j = 0; j < snps; j ++)
	{
	  if (haplo<1) {
	    s = g1[i][j]+g2[i][j];
	  } else {
	    s = g1[i][j];
	  }
	    if (s>eps){
		Lind[ig] = j;
		Lval[ig] = s;
		ig++;
	    }
	    
	}
	fprintf(pFile,"%d\n",ig); 
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%d ",Lind[j]);
	fprintf(pFile,"\n");
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%.1f ",Lval[j]);
	fprintf(pFile,"\n");

    }
    fclose (pFile);

    Rprintf("\n");    

    sst[0]=0;
    strcat(sst,outS);
    strcat(sst,"_matD.txt");
    pFile = fopen (sst,"w");
    if (pFile==NULL) {
      Rprintf("File >%s< cannot be opened! Stop.\n", sst);
      return(-1);
    }
    fprintf(pFile,"%d\n",individuals);  
    fprintf(pFile,"%d\n",(int) snps);  
    if (dsFound==1) {
    for(i = 0; i < individuals; i ++)
    {
      Rprintf("Write Individual Dosage: %d\r",(i+1));
	ig=0;
	for(j = 0; j < snps; j ++)
	{
	  s = ds[i][j];
	    if (s>eps){
		Lind[ig] = j;
		Lval[ig] = s;
		ig++;
	    }
	    
	}
	fprintf(pFile,"%d\n",ig); 
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%d ",Lind[j]);
	fprintf(pFile,"\n");
	for(j = 0; j < ig; j ++)
	    fprintf(pFile,"%.3f ",Lval[j]);
	fprintf(pFile,"\n");

    }
    } 
    fclose (pFile);


    Rprintf("\n");    
    Rprintf("\n");    



    return(0);
}

