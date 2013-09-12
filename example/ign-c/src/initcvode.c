#include "ign.h"

void initCVODE (void **cvode_mem, double tstart, double tEnd, N_Vector y0, int Nvars,
                double *udata, realtype relT, N_Vector absT, 
                int CVmaxord, int CVmaxnumsteps )
{
  int cvflag ;

  /* Create cvode solver */
  (*cvode_mem) = CVodeCreate(CV_BDF, CV_NEWTON);
  Check_CVflag(*cvode_mem, "CVodeCreate", 0) ;

  /* Allocate memory */
  cvflag = CVodeInit(*cvode_mem, &chemrhswrapper, tstart, y0);
  Check_CVflag(&cvflag, "CVodeInit", 1) ;
  cvflag = CVodeSVtolerances(*cvode_mem, relT, absT);
  Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;

  /*  Set dense solver */
  cvflag = CVDense(*cvode_mem, Nvars );
  Check_CVflag(&cvflag, "CVDense", 1) ;
  
  /* Set work array */
  cvflag = CVodeSetUserData(*cvode_mem,  (void *) udata);
  Check_CVflag(&cvflag, "CVodeSetUserData", 1) ;

#ifdef USEJAC
  /* Set dense Jacobian */
  cvflag = CVDlsSetDenseJacFn(*cvode_mem, &chemjacwrapper);
  Check_CVflag(&cvflag, "CVDlsSetDenseJacFn", 1) ;
#endif

  /* Set maximum order for the integratiom method */
  cvflag = CVodeSetMaxOrd(*cvode_mem, CVmaxord);
  Check_CVflag(&cvflag, "CVodeSetMaxOrd", 1) ;
  
  /* Set maximum number of steps */
  cvflag = CVodeSetMaxNumSteps(*cvode_mem, CVmaxnumsteps);
  Check_CVflag(&cvflag, "CVodeSetMaxNumSteps", 1) ;

  /* Set stop value */
  cvflag = CVodeSetStopTime(*cvode_mem, tEnd);
  Check_CVflag(&cvflag, "CVodeSetStopTime", 1) ;

#ifdef GETTIG
  /* Call CVodeRootInit to specify the root function g with 1 components */
  cvflag = CVodeRootInit(*cvode_mem, 1, &g);
  Check_CVflag(&cvflag, "CVodeRootInit", 1) ;
#endif

  return ;
}
