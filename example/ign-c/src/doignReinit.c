#include "ign.h"

double doIgnReinit(int getIgnDel, double tsta, double *tEnd, double *deltat, double deltatMax, double deltaTemp, 
             int oFreq, int *NiterMax, double *scal, double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
             FILE *myfile, FILE *myfile1, FILE *myfile2, FILE *myfile3, 
             void *cvode_mem)
{

  double t, tret, sumY; 
  int    i, iter, cvflag, Nspec, Nvars;
  double time_id1 = 0.0 ;

  Nvars = NV_LENGTH_S(y0) ;
#ifdef ALLSPEC
  Nspec = Nvars-1;
#else
  Nspec = Nvars  ;
#endif

  t    = tsta ;
  iter = 0    ;
  time_id1 = -100.0 ;

  /* Re-initialize cvode */
  cvflag = CVodeReInit(cvode_mem, tsta, y0 );
  Check_CVflag(&cvflag, "CVodeReInit", 1) ;
  cvflag = CVodeSVtolerances(cvode_mem, relT, absT);
  Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;
  cvflag = CVodeSetStopTime(cvode_mem, *tEnd);
  Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;

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
      /* printf("Found root : %20.12e\n",tret); */
      time_id1 = tret ;
      return ( time_id1 ) ; 
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
    (*deltat) = MIN((*deltat),deltatMax) ;
    scal[0] = NV_Ith_S(y0,0) ;

    /* re-calculate absolute tolerances based on recent solution */
    for( i = 0; i < Nvars ; i++)
    {
      double yscal = MAX( NV_Ith_S(y0,i), 0.0 ) ;
      NV_Ith_S(absT,i) = MAX( relT*yscal, CVsmall ) ;
    }

  }

  return ( time_id1 ) ;

} /* end of doIgn */
