#include "ign.h"

int main()
{

  /* Integer parameters */
  int Nspec, Nvars, NiterMax, oFreq    ;
  int CVmaxord, CVmaxnumsteps ;

  /* Double parameters */
  double Tini, Temp_id, deltat, deltatMax, tEnd ;
  double CVrelt, CVsmall ;
  double deltaTemp,pressure,pfac ;

  /* Character array parameters */
  int sizecMax = 100, specinno ;
  char *mechfile, *thermofile ;

  char   *SpecName ; 
  double *SpecMsFr ;

  unsigned int withTab ;  /* tabulation flag    */
  int getIgnDel ;         /* Ignition delay run */

  /* Scalar array */
  double *scal ;

  /* list of time steps */
  int    getDeltatSeq ;
  double *deltatSeq    ;

  /* work parameters */
  int i, kspec ;
  double sumY ;
  double time_id1, time_id2 ;

  /* CVODE-related parameters */
  int cvflag, iter, ierr ;
  N_Vector y0, absT ;
  realtype relT, tstart, tret ;
  double *udata ;
  double t ;


  mechfile   = (char *) malloc( sizecMax * sizeof(char) ) ; 
  thermofile = (char *) malloc( sizecMax * sizeof(char) );
  memset(mechfile  , 0, sizecMax) ;
  memset(thermofile, 0, sizecMax) ;

  SpecName = (char   *) malloc( MAXSPECIN * LENGTHOFSPECNAME * sizeof(char  ) ) ; 
  SpecMsFr = (double *) malloc( MAXSPECIN *                    sizeof(double) ) ; 
  memset(SpecName  , 0, sizecMax) ;
  for ( i = 0 ; i < MAXSPECIN; i++ ) SpecMsFr[i] = 0.0 ;

  setup(&NiterMax, &oFreq, &Tini, &Temp_id, &deltat, &deltatMax, &tEnd, &deltaTemp, 
        &pressure, &pfac, &withTab, &getIgnDel, 
        &CVrelt, &CVsmall, &CVmaxord, &CVmaxnumsteps, 
        SpecName, SpecMsFr, &specinno, 
        &getDeltatSeq, &deltatSeq, 
        mechfile, thermofile) ;

  /*
                           _               
                  ___  ___| |_ _   _ _ __  
                 / __|/ _ \ __| | | | '_ \ 
                 \__ \  __/ |_| |_| | |_) |
                 |___/\___|\__|\__,_| .__/ 
		 |_|     

  */
  /* Initialize TC library */
  TC_initChem( mechfile, thermofile, (int) withTab, 1.0) ; 
  free(mechfile  ) ;
  free(thermofile) ;

  /* Set pressure */
  pressure *= pfac ;
  TC_setThermoPres(pressure) ;

  Nspec = TC_getNspec() ;
  printf("Nspec   = %d\n",Nspec) ;
  scal = ( double * ) malloc( (Nspec+1)*sizeof(double) ) ;
 
  /* Set initial conditions */
  printf("Set initial conditions \n") ;
  scal[0] = Tini ;
  for ( i = 1 ; i<Nspec+1 ; i++) scal[i] = 0.0;

  for ( i = 0 ; i < specinno ; i++)
  {
    int ispec = TC_getSpos( &SpecName[i*LENGTHOFSPECNAME], strlen(&SpecName[i*LENGTHOFSPECNAME]) ) ;
    if (ispec < 0 )
    {
      printf("Error : Could not find species %s -> Exit !\n",&SpecName[i*LENGTHOFSPECNAME]);
      exit(1) ;
    }
    else
      printf("Index of species %s is %d\n",&SpecName[i*LENGTHOFSPECNAME],ispec) ;
    scal[ispec+1] = SpecMsFr[i] ;
  }
  free(SpecName  ) ;
  free(SpecMsFr  ) ;

  /* convert mole to mass fractions */
  double *msfr = (double *) malloc(Nspec * sizeof(double)) ;
  TC_getMl2Ms( &(scal[1]), Nspec, msfr ) ;
  for ( i = 1 ; i < Nspec+1 ; i++) scal[i] = msfr[i-1];
  free(msfr) ;

  for ( i = 1 ; i < Nspec+1 ; i++)
    if ( fabs(scal[i]) > 1.e-15 )
      printf("Mass fraction of species %d is %20.12e \n",i-1,scal[i]);

#ifdef ALLSPEC
  Nvars  = Nspec+1 ; /* T+(Nspec) */
#else
  Nvars  = Nspec   ; /* T+(Nspec-1) */
#endif
  tempNmsfr = (double *) malloc( (Nspec+1) * sizeof(double) ) ;
  rhsvals   = (double *) malloc( (Nspec+1) * sizeof(double) ) ;

  time_id1 = 0.0;
  time_id2 = 0.0 ;

  /* 
             CVODE setup 
  */
  tstart=0.0;   /* start time */

  /* initial condition array */
  y0 = N_VNew_Serial( Nvars ) ;

  /* cvode tolerances */
  relT = CVrelt ;
  absT = N_VNew_Serial( Nvars ) ;
  
  /* set initial conditions and absolute tolerances */
  for( i = 0; i < Nvars ; i++)
  {
    double yscal = MAX(scal[i],0.0) ;
    NV_Ith_S(y0  ,i) = yscal ;
    NV_Ith_S(absT,i) = MAX( relT*yscal, CVsmall ) ;
  }

  /* declare work space (mostly for jacobian) */
  udata = (double *) malloc ( (Nvars*Nvars+1)* sizeof(double))  ;
  udata[Nvars*Nvars] = Temp_id ;

  /* Create cvode solver */
  void *cvode_mem = NULL;
  initCVODE (&cvode_mem, tstart, tEnd, y0, Nvars, udata, relT, absT, 
             CVmaxord, CVmaxnumsteps ) ;

  printf("Starting time advancement : \n") ;
  t    = tstart ;
  iter = 0      ;
#ifndef NO_OUTPUT
  myfile  = fopen("ignsol.hdr","w" ) ;
  myfile1 = fopen("ys.hdr"    ,"w" ) ;
  myfile2 = fopen("cs.hdr"    ,"w" ) ;
  myfile3 = fopen("h.hdr"     ,"w" ) ;
  Output(0, iter, t, deltat, y0, scal) ;
  fclose(myfile);
  fclose(myfile1);
  fclose(myfile2);
  fclose(myfile3);

  myfile  = fopen("ignsol.dat","w" ) ;
  myfile1 = fopen("ys.out"    ,"w" ) ;
  myfile2 = fopen("cs.out"    ,"w" ) ;
  myfile3 = fopen("h.out"     ,"w" ) ;
  Output(1, iter, t, deltat, y0, scal) ;
#endif
  
  iter = NiterMax ;
  t    = tEnd     ;

#ifdef CHEMMODRUN
  /* retrieve the pre-exponential factor of reaction 10 */
  int reacid = 10 ;
  int posid  = 0 ;
  double preexp, preexpmod ;
  ierr = TC_getArhenFor(reacid,posid,&preexp) ;

  double *ysave    = (double *) malloc(Nvars * sizeof(double));
  double *absTsave = (double *) malloc(Nvars * sizeof(double));
  int iscal ;
  for ( iscal = 0; iscal < Nvars; iscal++ )
  {
    ysave[iscal]    = NV_Ith_S(y0,iscal  ) ;
    absTsave[iscal] = NV_Ith_S(absT,iscal) ;
  }
  for ( i=0; i<10; i++)
  {
    for ( iscal = 0; iscal < Nvars; iscal++ )
    {
      NV_Ith_S(y0,iscal)    = ysave[iscal] ;
      NV_Ith_S(absT,iscal)  = absTsave[iscal] ;
    }
    deltat    = 10.0 ; /* seconds */
    deltatMax = 10.0 ; /* seconds */
    deltaTemp = 100.0 ; /* Kelvin */
    preexpmod = ((double) i+1)*10.0*preexp ;
    ierr = TC_chgArhenFor(reacid, posid, preexpmod) ;
    time_id1 = doIgnReinit( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
			    oFreq, &iter, scal, CVsmall, y0, relT, absT, 
			    myfile, myfile1, myfile2, myfile3, cvode_mem) ;
    printf("Ignition %d done : %20.12e\n",i,time_id1) ;
  }
#else
  if ( getDeltatSeq > 0 )
    time_id1 = doIgnDtSeq( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
		      oFreq, &iter, scal, 
                      getDeltatSeq, deltatSeq, 
                      CVsmall, y0, relT, absT, 
                      myfile, myfile1, myfile2, myfile3, cvode_mem) ;
  else
    time_id1 = doIgn( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
                      oFreq, &iter, scal, 
                      CVsmall, y0, relT, absT, 
                      myfile, myfile1, myfile2, myfile3, cvode_mem) ;
#endif


#ifndef NO_OUTPUT
  Output(1, iter, t, deltat, y0, scal) ;
  fclose(myfile) ;
  fclose(myfile1) ;
  fclose(myfile2) ;
  fclose(myfile3) ;  
#endif

  FILE *myfileID=fopen( "tid.dat", "w" );
  fprintf(myfileID,"%20.12e  %20.12e\n",time_id1,time_id2) ;
  printf("%20.12e  %20.12e\n",time_id1,time_id2) ;
  fclose(myfileID) ;

  /* Clean-up cvode */
  N_VDestroy_Serial( y0   ) ;
  N_VDestroy_Serial( absT ) ;
  CVodeFree( &cvode_mem ) ;

  TC_reset() ;

  free(tempNmsfr ) ;
  free(rhsvals   ) ;
  free(udata     ) ;
  free(scal      ) ;

  if ( deltatSeq != NULL ) free(deltatSeq) ;

  return ( 0 ) ;

}



