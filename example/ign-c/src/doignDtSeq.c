#include "ign.h"

double doIgnDtSeq(int getIgnDel, double tsta, double *tEnd, double *deltat, double deltatMax, 
                  double deltaTemp, int oFreq, int *NiterMax, double *scal, 
                  int getDeltatSeq, double *deltatSeq,
                  double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
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
  while ( iter < getDeltatSeq ) 
  {

/*#ifdef DONOTDOIT*/
    /* Re-initialize cvode */
    cvflag = CVodeReInit(cvode_mem, t, y0 );
    Check_CVflag(&cvflag, "CVodeReInit", 1) ;

    cvflag = CVodeSVtolerances(cvode_mem, relT, absT);
    Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;

    cvflag = CVodeSetStopTime(cvode_mem, *tEnd);
    Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;
/*#endif*/

    /* call cvode */
    (*deltat) = deltatSeq[iter] ;
    t    += deltatSeq[iter] ;
    iter += 1               ;
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

  } /* end of while loop */

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
    sumY += scal[i] ;
  }
  scal[Nspec] = 1.0-sumY ;
#endif

  (*NiterMax) = iter ; 
  (*tEnd)     = tret ;

  return ( time_id1 ) ;

} /* end of doIgn */
