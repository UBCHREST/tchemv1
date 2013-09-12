#include "ign.h"

double doIgn(int getIgnDel, double tsta, double *tEnd, double *deltat, double deltatMax, double deltaTemp, 
             int oFreq, int *NiterMax, double *scal, double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
             FILE *myfile, FILE *myfile1, FILE *myfile2, FILE *myfile3, 
             void *cvode_mem)
{

  double t, tret, sumY; 
  int    i, iter, cvflag, Nspec, Nvars;
  
  int foundid2 = 0 ;
  double time_id1 = 0.0    , time_id2 = 0.0     ;
  double Temp_m2  = scal[0], Temp_m1  = scal[0] ;
  double time_m2  = 0.0,     time_m1  = 0.0     ;
  double der2 ;

  Nvars = NV_LENGTH_S(y0) ;
#ifdef ALLSPEC
  Nspec = Nvars-1;
#else
  Nspec = Nvars  ;
#endif

  t    = tsta ;
  iter = 0    ;
  while ( ( t < (*tEnd) ) && (iter < (*NiterMax) ) )
  {
    /* call cvode */
    t    += (*deltat) ;
    iter += 1         ;
    cvflag = CVode(cvode_mem, t, y0, &tret, CV_NORMAL);

    if (cvflag == CV_ROOT_RETURN) 
    {
      int rootsfound[1] ;
      int flagroot = CVodeGetRootInfo(cvode_mem, rootsfound);
      Check_CVflag(&flagroot, "CVodeGetRootInfo", 1) ;
      printf("Found root : %20.12e\n",tret);
      time_id1 = tret ;
      if ( getIgnDel == 1 ) { *tEnd = tret ; (*NiterMax) = iter ; return ( time_id1 ) ; }
    }

    /* reset time based on end time from cvode */
    t = tret ;

    /*
       calculate new time step:
       - reduce/increase time step if \Delta T >< 0.25K
     */
    (*deltat) = MIN(
		    deltaTemp/MAX(fabs(NV_Ith_S(y0,0)-scal[0]),TCSMALL)*(*deltat),
		    2.0*(*deltat)
		);

     /* reset time step based on set limits */
     /* if ( iter < 100 ) */
     /*   (*deltat) = MIN((*deltat),1.e-4) ; */
     /* else */
       (*deltat) = MIN((*deltat),deltatMax) ;


    scal[0] = NV_Ith_S(y0,0) ;

#ifndef NO_OUTPUT
    /* output solution */
    if ( iter % oFreq == 0 )
      Output(1, iter, t, *deltat, y0, scal) ;
#endif

    /* re-calculate absolute tolerances based on recent solution */
    for( i = 0; i < Nvars ; i++)
    {
      double yscal = MAX( NV_Ith_S(y0,i), 0.0 ) ;
      NV_Ith_S(absT,i) = MAX( relT*yscal, CVsmall ) ;
    }

    /* Ignition delay time - second derivative */
    if ( ( foundid2 == 0 ) && ( NV_Ith_S(y0,0) > 1300.0 ) )
    {
      der2 = (NV_Ith_S(y0,0)-Temp_m1)/(t-time_m1)-(Temp_m1-Temp_m2)/(time_m1-time_m2) ;
      if ( der2 <= 0.0 )
      {
    	foundid2 = 1 ;
    	time_id2 = time_m1 ;
	printf("Found time_id2 : %20.12e\n",time_id2) ;
      }
    }

    /* Cycle temperatures and times */
    Temp_m2 = Temp_m1        ;
    Temp_m1 = NV_Ith_S(y0,0) ;
    time_m2 = time_m1        ;
    time_m1 = t              ;

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

  (*NiterMax) = iter ; 
  *tEnd       = tret ;

  return ( time_id1 ) ;

} /* end of doIgn */
