#include "ign.h"

void Output(int iflag, int iter, double t, double deltat, N_Vector y0, double *scal)
{
  int i, ierr, Nvars, Nspec ;
  double sumY;

  Nvars = NV_LENGTH_S(y0) ;
#ifdef ALLSPEC
  Nspec = Nvars-1;
#else
  Nspec = Nvars  ;
#endif

  if ( tempNmsfr == NULL ) 
    tempNmsfr = (double *) malloc( (Nspec+1) * sizeof(double) ) ;
  if ( rhsvals   == NULL ) 
    rhsvals   = (double *) malloc( (Nspec+1) * sizeof(double) ) ;

  if ( iflag == 0 ) 
  {
    char *snameshdr = (char *) malloc( Nspec * LENGTHOFSPECNAME * sizeof(char) ) ; 
    ierr = TC_getSnames( Nspec, snameshdr ) ;
    fprintf(myfile,"# 1: iteration, #2: time, #3: deltat, #4: temperature,\n") ;
    for ( i = 0; i<Nspec; i++ )
      fprintf(myfile,"#%d: Mass fraction of %s \n",
              i+5,&snameshdr[i*LENGTHOFSPECNAME]);
    fprintf(myfile1,"#1: time, #2: temperature,\n") ;
    for ( i = 0; i<Nspec; i++ )
      fprintf(myfile1,"#%d: Mass fraction of %s \n",
              i+3,&snameshdr[i*LENGTHOFSPECNAME]);
    fprintf(myfile2,"#1: time, \n") ;
    for ( i = 0; i<Nspec; i++ )
      fprintf(myfile2,"#%d: Molar concentration of %s [kmol/m3]\n",
             i+2,&snameshdr[i*LENGTHOFSPECNAME]);
    fprintf(myfile3,"# 1: time, #2: specific enthalpy [J/kg]") ;
    free(snameshdr) ;
    return ;
  }

#ifdef ALLSPEC
  /* Integrate T+(Nspec) */
  for (i = 0; i < Nvars; i++) scal[i] = NV_Ith_S(y0,i) ;
#else
  /* Integrate T+(Nspec-1); make sure Y_Nspec = 1-sum(Y_i) */
  scal[0] = NV_Ith_S(y0,0) ;
  sumY = 0.0 ;
  for ( i = 1; i < Nvars; i++)
  {
    scal[i] = NV_Ith_S(y0,i) ;
    sumY += NV_Ith_S(y0,i) ;
  }
  scal[Nspec] = 1.0-sumY ;
#endif

  /* Output */
  printf(" Iter = %-7d,  t[s] = %14.6e, dt[s] = %14.6e, T[K] = %14.6e\n",
         iter, t, deltat, scal[0] );

  fprintf(myfile,"%-10d  %20.12e  %20.12e", iter, t, deltat );
  for ( i = 0 ; i<Nspec+1 ; i++)   fprintf(myfile,"  %20.12e", scal[i] );
  fprintf(myfile,"\n");

  fprintf(myfile1,"%20.12e", t);
  for ( i = 0 ; i<Nspec+1 ; i++)   fprintf(myfile1,"  %20.12e", scal[i] );
  fprintf(myfile1,"\n");
  
  TC_getMs2Cc ( scal, Nspec+1, &tempNmsfr[1] ) ;
  fprintf(myfile2,"%20.12e  %20.12e", t, scal[0]);
  for ( i = 1 ; i<Nspec+1 ; i++)   fprintf(myfile2,"  %20.12e", tempNmsfr[i] );
  fprintf(myfile2,"\n");

  TC_getMs2HmixMs ( scal, Nspec+1, &tempNmsfr[0] ) ;
  fprintf(myfile3,"%20.12e  %20.12e", t, tempNmsfr[0]);

  return ;

}
