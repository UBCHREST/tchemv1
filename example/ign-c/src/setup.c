#include "ign.h"

#define FNAMEMAX 100

void setup(int *NiterMax, int *oFreq, double *Tini, double *Temp_id, 
           double *deltat, double *deltatMax, double *tEnd, double *deltaTemp, 
           double *pressure, double *pfac, unsigned int *withTab, int *getIgnDel, int *getSens, 
           double *CVrelt, double *CVsmall, int *CVmaxord, int *CVmaxnumsteps, 
           char *SpecName, double *SpecMsFr, int *specinno,
           int *getDeltatSeq, double **deltatSeq, 
           char *mechfile, char *thermofile)
{
  int i ;
  char filDeltatSeq[FNAMEMAX] ;

  /* 
      _       __             _ _       
   __| | ___ / _| __ _ _   _| | |_ ___ 
  / _` |/ _ \ |_ / _` | | | | | __/ __|
 | (_| |  __/  _| (_| | |_| | | |_\__ \
  \__,_|\___|_|  \__,_|\__,_|_|\__|___/

  */
  *NiterMax  = 100000    ; /* Maximum no. of iterations */
  *oFreq     = 10        ; /* Output frequency          */
  *Tini      = 1000.0    ; /* Initial temperature  [K]  */
  *Temp_id   = 1500.0    ; /* Threshold temperature for ignition delay [K]  */
  *deltat    = 1.e-10    ; /* Initial time step size [s]  */
  *deltatMax = 1.e-4     ; /* Maximum time step size [s]  */
  *tEnd      = 2.0e0     ; /* End integration time   [s]  */
  *deltaTemp = 1.0       ; /* Maximum temperature change per time step [K]  */
  *pressure  = 1.01325e5 ; /* Pressure [Pa] */
  *pfac      = 1.0       ; /* Pressure scale factor [ ]  */

  *getDeltatSeq = 0      ; /* Do not use a sequence of dt  */
  *deltatSeq    = NULL   ;

  *CVrelt        = 1.e-12 ;
  *CVsmall       = 1.e-20 ;
  *CVmaxord      = 5      ;
  *CVmaxnumsteps = 10000  ;

  strcpy(mechfile  ,"chem.inp" ) ; 
  strcpy(thermofile,"therm.dat") ;

  *withTab = 0 ;  /* no tabulation           */
  *getIgnDel = 0 ;/* no ignition delay stop  */
  *getSens   = 0 ;/* no sensitivity analysis */ 

  /* read setup file */
  FILE *fsetup = fopen("input.dat","r") ;
  char key[20] ;
  *specinno = 0 ;
  while ( fscanf(fsetup, "%s", key) != EOF ) 
  {
    if ( strncmp(key,"END",3) == 0 ) break ;
    if ( strncmp(key,"NiterMax"    , 8) == 0 ) fscanf(fsetup, "%d" , NiterMax   ) ;
    if ( strncmp(key,"oFreq"       , 5) == 0 ) fscanf(fsetup, "%d" , oFreq      ) ;
    if ( strncmp(key,"withTab"     , 7) == 0 ) fscanf(fsetup, "%d" , withTab    ) ;
    if ( strncmp(key,"getIgnDel"   , 9) == 0 ) fscanf(fsetup, "%d" , getIgnDel  ) ;
    if ( strncmp(key,"getSens"     , 7) == 0 ) fscanf(fsetup, "%d" , getSens    ) ;
    if ( strncmp(key,"Tini"        , 4) == 0 ) fscanf(fsetup, "%le", Tini       ) ;
    if ( strncmp(key,"Temp_id"     , 7) == 0 ) fscanf(fsetup, "%le", Temp_id    ) ;
    if ( strncmp(key,"deltatMax"   , 9) == 0 ) fscanf(fsetup, "%le", deltatMax  ) ;
    if ( strncmp(key,"deltat"      , 6) == 0 ) fscanf(fsetup, "%le", deltat     ) ;
    if ( strncmp(key,"deltaTemp"   , 9) == 0 ) fscanf(fsetup, "%le", deltaTemp  ) ;
    if ( strncmp(key,"tEnd"        , 4) == 0 ) fscanf(fsetup, "%le", tEnd       ) ;
    if ( strncmp(key,"pfac"        , 4) == 0 ) fscanf(fsetup, "%le", pfac       ) ;
    if ( strncmp(key,"mech"        , 4) == 0 ) fscanf(fsetup, "%s" , mechfile   ) ;
    if ( strncmp(key,"thermo"      , 6) == 0 ) fscanf(fsetup, "%s" , thermofile ) ;

    if ( strncmp(key,"getDeltatSeq",12) == 0 ) fscanf(fsetup, "%d" , getDeltatSeq ) ;
    if ( strncmp(key,"filDeltatSeq",12) == 0 ) fscanf(fsetup, "%s" , filDeltatSeq ) ;

    if ( strncmp(key,"spec"      ,4) == 0 ) 
    {  
      fscanf(fsetup, "%s%le" ,&SpecName[(*specinno)*LENGTHOFSPECNAME],&SpecMsFr[*specinno]) ; 
      (*specinno)++ ;
    }    
  }
  fclose(fsetup) ;

  /* Convert species names to upper case */
  for ( i = 0 ; i<(*specinno)*LENGTHOFSPECNAME ; i++)
    SpecName[i] = toupper(SpecName[i]);

  printf("------------------------------------------\n") ;
  printf("Run parameters : \n") ;
  printf("    NiterMax  = %d\n"     ,*NiterMax  ) ;
  printf("    oFreq     = %d\n"     ,*oFreq     ) ;
  printf("    Tini      = %20.12e\n",*Tini      ) ;
  printf("    Temp_id   = %20.12e\n",*Temp_id   ) ;
  printf("    deltat    = %20.12e\n",*deltat    ) ;
  printf("    deltatMax = %20.12e\n",*deltatMax ) ;
  printf("    deltaTemp = %20.12e\n",*deltaTemp ) ;
  printf("    tEnd      = %20.12e\n",*tEnd      ) ;
  printf("    pfac      = %20.12e\n",*pfac      ) ;
  printf("    Kinetic model : %s\n",mechfile   ) ;
  printf("    Thermo props  : %s\n",thermofile ) ;
  printf("    Create tables : %d\n",*withTab    ) ;
  printf("    Do ign del    : %d\n",*getIgnDel  ) ;
  printf("    Do sensitivity: %d\n",*getSens    ) ;
  printf("    Species mole fractions : \n"    ) ;
  for ( i = 0 ; i<*specinno ; i++)
    if ( fabs(SpecMsFr[i]) > 1.e-15 )
      printf("    X_%s = %20.12e\n",&SpecName[i*LENGTHOFSPECNAME],SpecMsFr[i]);
  printf("------------------------------------------\n") ;

  /* Read sequence of time steps from file if necessary */
  if ( (*getDeltatSeq) == 1 )
  {
    FILE *fseq = NULL ;
    fseq = fopen(filDeltatSeq,"r");
    if ( fseq != NULL )
    {
      fscanf(fseq,"%d",getDeltatSeq);
      (*deltatSeq) = (double *) malloc((*getDeltatSeq)*sizeof(double)) ;
      for ( i=0; i<(*getDeltatSeq); i++ ) fscanf(fseq,"%le",&((*deltatSeq)[i])) ;
      fclose(fseq) ;

      printf("There are %d entries in the list of deltat : \n",(*getDeltatSeq));
      for ( i=0; i<(*getDeltatSeq); i++ ) printf("%d: %24.16e\n",i,(*deltatSeq)[i]) ;
      printf("------------------------------------------\n") ;

    }
    else
    {
      printf("Error in setup() : Could not open %s -> Abort\n",filDeltatSeq) ;
      exit(-1);
    }
  } /* done if (*getDeltatSeq) == 1 */

  return ;

}
