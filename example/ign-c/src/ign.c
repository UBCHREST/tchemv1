#include "ign.h"

void getmeps(double *rtol, double *atol) 
{

  double sml = 1.0 ;
  while ( 1.0+sml != 1.0 ) sml *= 0.5 ;
  *rtol = sqrt(2.0*sml) ;
  *atol = *rtol ;
  return ;

}
int main()
{

  /* Integer parameters */
  int Nspec, Nvars, Nreac, NiterMax, oFreq    ;
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
  int getSens   ;         /* Sensitivity run    */

  /* Scalar array */
  double *scal ;

  /* list of time steps */
  int    getDeltatSeq ;
  double *deltatSeq    ;

  /* work parameters */
  int i, kspec ;
  double sumY ;
  double time_id1,time_id2 ;

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
        &pressure, &pfac, &withTab, &getIgnDel, &getSens, 
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
  Nreac = TC_getNreac() ;
  printf("Nspec/Nreac = %d/%d\n",Nspec,Nreac) ;
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
  myfile  = fopen("msfrini.dat","w" ) ;
  for ( i = 0 ; i < Nspec ; i++) fprintf(myfile,"%20.12e\n",msfr[i]);  
  fclose(myfile);
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
  time_id1 = 0.0  ;
  time_id2 = 0.0  ;

  if ( getDeltatSeq > 0 )
    time_id1 = doIgnDtSeq( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
		      oFreq, &iter, scal, 
                      getDeltatSeq, deltatSeq, 
                      CVsmall, y0, relT, absT, 
                      myfile, myfile1, myfile2, myfile3, cvode_mem) ;
  else
    if ( getSens > 0 )
    {

      double rtol, atol;
      int iscal, ir ;
      getmeps(&rtol,&atol) ;

      double *ysave    = (double *) malloc(Nvars * sizeof(double));
      double *absTsave = (double *) malloc(Nvars * sizeof(double));
      double arrh, arrhmod, perturb ;
      for ( iscal = 0; iscal < Nvars; iscal++ )
      {
	ysave[iscal]    = NV_Ith_S(y0,iscal  ) ;
	absTsave[iscal] = NV_Ith_S(absT,iscal) ;
      }

      for ( ir = 0; ir < Nreac; ir++) 
      {
        /* restore IC */
	for ( iscal = 0; iscal < Nvars; iscal++ )
	{
	  NV_Ith_S(y0,iscal)    = ysave[iscal] ;
	  NV_Ith_S(absT,iscal)  = absTsave[iscal] ;
	}
        /* get Arrhenius factor and perturb "-" */
	ierr = TC_getArhenFor(ir,getSens-1,&arrh) ;
	perturb = fabs(rtol*arrh) ;
	if ( perturb == 0.0 ) perturb = atol ;
	arrhmod = arrh-perturb;
	ierr = TC_chgArhenFor(ir,getSens-1, arrhmod) ;
	iter = NiterMax ;
	t    = tEnd     ;
	time_id1 = doIgn( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
			oFreq, &iter, scal, CVsmall, y0, relT, absT, 
			myfile, myfile1, myfile2, myfile3, cvode_mem) ;

        /* restore IC */
	for ( iscal = 0; iscal < Nvars; iscal++ )
	{
	  NV_Ith_S(y0,iscal)    = ysave[iscal] ;
	  NV_Ith_S(absT,iscal)  = absTsave[iscal] ;
	}

        /* get perturb "+" */
 	arrhmod = arrh+perturb;
	ierr = TC_chgArhenFor(ir, getSens-1, arrhmod) ;
	iter = NiterMax ;
	t    = tEnd     ;
	time_id2 = doIgn( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
			oFreq, &iter, scal, CVsmall, y0, relT, absT, 
			myfile, myfile1, myfile2, myfile3, cvode_mem) ;
	double sval = (time_id2-time_id1)/(2.0*perturb); 
	printf("Reaction %d: %20.12e %20.12e, %20.12e %20.12e\n", ir, time_id1,time_id2, sval, sval*arrh);
	ierr = TC_chgArhenFor(ir, getSens-1, arrh) ;
      } /* done loop over all reactions */
    } /* done if sensitivity analysis */
    else
      time_id1 = doIgn( getIgnDel, tstart, &t, &deltat, deltatMax, deltaTemp, 
			oFreq, &iter, scal, 
			CVsmall, y0, relT, absT, 
			myfile, myfile1, myfile2, myfile3, cvode_mem) ;


#ifndef NO_OUTPUT
  Output(1, iter, t, deltat, y0, scal) ;
  fclose(myfile) ;
  fclose(myfile1) ;
  fclose(myfile2) ;
  fclose(myfile3) ;  
#endif

  FILE *myfileID=fopen( "tid.dat", "w" );
  fprintf(myfileID,"%20.12e\n",time_id1) ;
  printf("%20.12e\n",time_id1) ;
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



