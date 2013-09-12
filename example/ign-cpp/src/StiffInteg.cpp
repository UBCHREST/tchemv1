#include "StiffInteg.h"

// Static members need to be pre-declared
int StiffInteg::next_index_ = 0;
StiffInteg::OMap_t *StiffInteg::omap_ = NULL;

int StiffInteg::compute(double tend, double *deltat, double deltatMax, 
			int NiterMax, int oFreq)
{
   
  if ( !isInit_ ) this->init() ;

  /* Maximum allowed temperature jump */
  double deltaTemp = 1.0 ;

  /* Initialize time and temperature history */
  bool foundid2 = false ;
  double der2num = 0.0 ;
  double Temp_m2  = scal_[0], Temp_m1 = scal_[0] ; 
  double time_m2  = 0.0,      time_m1 = 0.0 ;
  double time_id1 = 0.0,      time_id2 = 0.0 ;


  /* initial condition array */
  N_Vector y0 ;
  y0 = N_VNew_Serial( Nvars_ ) ;

  /* cvode tolerances */
  realtype relT ;
  N_Vector absT ;
  relT = CVrelt_ ;
  absT = N_VNew_Serial( Nvars_ ) ;
  
  /* set initial conditions and absolute tolerances */
  for(int i = 0; i < Nvars_ ; i++)
  {
    double yscal = MAX(scal_[i],0.0) ;
    NV_Ith_S(y0  ,i) = yscal ;
    NV_Ith_S(absT,i) = MAX( relT*yscal, CVsmall_ ) ;
  }

  /* declare integration limits */
  realtype tstart=0.0, tret ;

  /* declare user space for jacobian */
  double *udata ;
  udata = new double[Nvars_*Nvars_+2] ;
  udata[0] = (double) ( this->my_index_ ) ;
  udata[Nvars_*Nvars_+1] = Temp_id ;

  /* cvode flag */
  int cvflag ;

  /* Create cvode solver */
  void *cvode_mem = NULL;
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  this->Check_CVflag(cvode_mem, "StiffInteg::compute : CVodeCreate", 0) ;
  
  /* Allocate memory */
  cvflag = CVodeInit(cvode_mem, &chemrhswrapper, tstart, y0);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeInit", 1) ;
  cvflag = CVodeSVtolerances(cvode_mem, relT, absT);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeSVtolerances", 1) ;

  /* Set dense solver */
  cvflag = CVDense(cvode_mem, Nvars_ );
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVDense", 1) ;
  
  /* Set work array */
  cvflag = CVodeSetUserData(cvode_mem, (void *) udata);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeSetUserData", 1) ;

#ifdef USEJAC
  /* Set dense Jacobian */
  cvflag = CVDlsSetDenseJacFn(cvode_mem, &chemjacwrapper);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVDlsSetDenseJacFn", 1) ;
#endif

  /* Set maximum order for the integratiom method */
  cvflag = CVodeSetMaxOrd(cvode_mem, CVmaxord_);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeSetMaxOrd", 1) ;
  
  /* Set maximum number of steps */
  cvflag = CVodeSetMaxNumSteps(cvode_mem, CVmaxnumsteps_);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeSetMaxNumSteps", 1) ;

  /* Set stop value */
  cvflag = CVodeSetStopTime(cvode_mem, tend);
  this->Check_CVflag(&cvflag, "StiffInteg::compute : CVodeSetStopTime", 1) ;

  /* Call CVodeRootInit to specify the root function g with 1 components */
  cvflag = CVodeRootInit(cvode_mem, 1, &findignwrap);
  this->Check_CVflag(&cvflag, "CVodeRootInit", 1) ;

  double t = tstart ;
  int    iter = 0 ;

#ifndef NO_OUTPUT
  /* output initial condition */
  std::ofstream myfile ( "ignsol.dat" );
  myfile.precision(12) ;
  myfile<<std::setw(10)<<iter<<std::scientific<<" "<<t<<" "<<*deltat<<" " ;
  for (int i = 0 ; i<Nspec_+1 ; i++) myfile <<scal_[i]<<" " ;
  myfile<<std::endl ;
  
  std::ofstream myfile1 ( "ys.out" );
  myfile1.precision(12) ;
  myfile1<<std::scientific<<t<<" " ;
  for (int i = 0 ; i<Nspec_+1 ; i++) myfile1<<scal_[i]<<" " ;
  myfile1<<std::endl ;
  
  std::ofstream myfile2 ( "cs.out" );
  myfile2.precision(12) ;
  TC_getMs2Cc ( scal_, Nspec_+1, &tempNmsfr[1] ) ;
  myfile2<<std::scientific<<t<<" "<<scal_[0]<<" " ;
  for (int i = 1 ; i<Nspec_+1 ; i++) myfile2 << tempNmsfr[i]<<" " ;
  myfile2<<std::endl ;
  
  std::ofstream myfile3 ( "h.out" );
  myfile3.precision(12) ;
  TC_getMs2HmixMs ( scal_, Nspec_+1, &tempNmsfr[0] ) ;
  myfile3<<std::scientific<<t<<" "<<tempNmsfr[0] <<std::endl ;
#endif
  
  std::cout<<"Starting time advancement : "<<std::endl ;

  while ( ( t < tend ) && (iter < NiterMax ) )
  {
    /* call cvode */
    t    += (*deltat) ;
    iter += 1         ;
    cvflag = CVode(cvode_mem, t, y0, &tret, CV_NORMAL);

    if (cvflag == CV_ROOT_RETURN) {
      int rootsfound[1] ;
      int flagroot = CVodeGetRootInfo(cvode_mem, rootsfound);
      this->Check_CVflag(&flagroot, "CVodeGetRootInfo", 1) ;
      printf("Found time_id1 : %20.12e\n",tret);
      time_id1 = tret ;
    }

    /* reset time based on end time from cvode */
    t = tret ;

    /* 
       calculate new time step:
       - reduce/increase time step if \Delta T >< 0.25K 
     */
    *deltat = MIN(
		  dTmax_/MAX(fabs(NV_Ith_S(y0,0)-scal_[0]),TCSMALL)*(*deltat),
		  2.0*(*deltat)
		  );

    /* reset time step based on set limits */
    if ( iter < 100 ) 
      *deltat = MIN(*deltat,1.e-4) ;
    else
      *deltat = MIN(*deltat,deltatMax) ;

#ifdef ALLSPEC
    /* Integrate T+(Nspec) */
    for (int i = 0; i < Nvars_; i++)
      scal_[i] = NV_Ith_S(y0,i) ;
#else
    /* Integrate T+(Nspec-1); make sure Y_Nspec = 1-sum(Y_i) */
    scal_[0] = NV_Ith_S(y0,0) ;
    double sumY = 0.0 ;
    for (int i = 1; i < Nvars_; i++)
    {
      scal_[i] = NV_Ith_S(y0,i) ;
      sumY += NV_Ith_S(y0,i) ;
    }
    scal_[Nspec_] = 1.0-sumY ;
#endif

#ifndef NO_OUTPUT
    /* output solution */
    if ( iter % oFreq == 0 )
    {
      std::cout<<" Iter = "<<std::setw(7)<<iter<<"  t[s] = "
	       <<std::fixed<<std::scientific<<std::setprecision(6)
	       <<t<<"  dt[s] = "<<(*deltat)<<"   T[K] = "<<scal_[0]<<std::endl ;

      myfile<<std::setw(10)<<iter<<" "<<std::scientific<<t<<" "<<*deltat<<" " ;
      for (int i = 0 ; i < Nspec_+1 ; i++) myfile << scal_[i]<<" " ;
      myfile<<std::endl ;

      myfile1<<std::scientific<<t<<" " ;
      for (int i = 0 ; i<Nspec_+1 ; i++) myfile1 << scal_[i]<<" " ;
      myfile1<<std::endl ;
  
      TC_getMs2Cc ( scal_, Nspec_+1, &tempNmsfr[1] ) ;
      myfile2<<std::scientific<<t<<" "<<scal_[0]<<" " ;
      for (int i = 1 ; i<Nspec_+1 ; i++) myfile2 << tempNmsfr[i]<<" " ;
      myfile2<<std::endl ;
  
      TC_getMs2HmixMs ( scal_, Nspec_+1, &tempNmsfr[0] ) ;
      myfile3<<std::scientific<<t<<" "<<tempNmsfr[0] <<std::endl ;
  
    }
#endif

    /* calculate absolute tolerances based on recent solution */
    for(int i = 0; i < Nvars_ ; i++)
    {
      double yscal = MAX(scal_[i],0.0) ;
      NV_Ith_S(absT,i) = MAX( relT*yscal, CVsmall_ ) ;
    }

    /* Ignition delay time - second derivative */
    if ( !foundid2 && ( scal_[0] > 1300.0 ) )
    {
      der2num = (scal_[0]-Temp_m1)/(t-time_m1)-(Temp_m1-Temp_m2)/(time_m1-time_m2) ;
      if ( der2num<=0.0 )
      {
	foundid2 = true ;
	time_id2 = time_m1 ;
	printf("Found time_id2 : %20.12e\n",time_id2);
      }
    }

    /* Cycle temperatures and times */
    Temp_m2 = Temp_m1 ;
    Temp_m1 = scal_[0] ;
    time_m2 = time_m1 ;
    time_m1 = t ;

  }


#ifndef NO_OUTPUT
  /* Output last soution */
  std::cout<<" Iter = "<<std::setw(7)<<iter<<"  t[s] = "
	   <<std::scientific<<std::setprecision(6)
	   <<t<<"  dt[s] = "<<(*deltat)<<"   T[K] = "<<scal_[0]<<std::endl ;

  myfile<<std::setw(10)<<iter<<" "<<std::scientific<<t<<" "<<*deltat<<" " ;
  for (int i = 0 ; i<Nspec_+1 ; i++) myfile << scal_[i]<<" " ;
  myfile<<std::endl ;

  myfile1<<std::scientific<<t<<" " ;
  for (int i = 0 ; i<Nspec_+1 ; i++) myfile1 << scal_[i]<<" " ;
  myfile1<<std::endl ;
  
  TC_getMs2Cc ( scal_, Nspec_+1, &tempNmsfr[1] ) ;
  myfile2<<std::scientific<<t<<" "<<scal_[0]<<" " ;
  for (int i = 1 ; i<Nspec_+1 ; i++) myfile2 << tempNmsfr[i]<<" " ;
  myfile2<<std::endl ;
  
  TC_getMs2HmixMs ( scal_, Nspec_+1, &tempNmsfr[0] ) ;
  myfile3<<std::scientific<<t<<" "<<tempNmsfr[0] <<std::endl ;
  
  /* close files */
  myfile. close() ;
  myfile1.close() ;
  myfile2.close() ;
  myfile3.close() ;
#endif

  std::ofstream myfileID ( "tid.dat" );
  myfileID.precision(12) ;
  myfileID<<time_id1<<"  "<<time_id2<<std::endl<<std::flush ;
  myfileID.close() ;
  std::cout<<time_id1<<"  "<<time_id2<<std::endl<<std::flush ;

  /* Destroy temporary arrays */
  N_VDestroy_Serial( y0   ) ;
  N_VDestroy_Serial( absT ) ;

  /* destroy solver */
  CVodeFree( &cvode_mem ) ;

  /* Cleanup */
  delete [] udata ; udata = NULL ;
  return ( 0 ) ;

}

int StiffInteg::init()
{

  /* cvode parameters */
  CVrelt_        = 1.e-10 ;
  CVsmall_       = 1.e-20 ;
  CVmaxord_      = 5      ;
  CVmaxnumsteps_ = 10000  ;

  /* temperature and species mass fractions */
  tempNmsfr = new double[Nvars_+1] ;
  rhsvals   = new double[Nvars_+1] ;

  /* default temperature threshold for ignition delay */
  Temp_id = 1500.0 ;

  useJacAnl = true ;
  isInit_   = true ;

  return ( 0 ) ;

}

void StiffInteg::chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{

  /* temperature */
  tempNmsfr[0] = NV_Ith_S(y,0) ;

#ifdef ALLSPEC
  /* here Nvars_ = 1+(Nspec_) */
  for (int iscal = 1; iscal < Nvars_; iscal++)
    tempNmsfr[iscal] = NV_Ith_S(y,iscal) ;
#else
  /* here Nvars_ = 1+(Nspec_-1) */
  double sumY = 0.0 ;
  for (int iscal = 1; iscal < Nvars_; iscal++)
  {
    tempNmsfr[iscal] = NV_Ith_S(y,iscal) ;
    sumY += tempNmsfr[iscal] ;
  }
  tempNmsfr[Nvars_] = 1.0-sumY ;
#endif

  /* get src term from tchem library */
  TC_getSrc ( tempNmsfr, Nspec_+1, rhsvals ) ;

  for (int iscal = 0; iscal < Nvars_ ; iscal++)
    NV_Ith_S(ydot,iscal) = rhsvals[iscal] ;

  return ;

}

#ifdef USEJAC
void StiffInteg::chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp)
{

  /* temperature */
  tempNmsfr[0] = NV_Ith_S(y,0) ;

#ifdef ALLSPEC
  /* Nvars_ = T+(Nspec_) */
  for (int i = 1; i < Nvars_; i++)
    tempNmsfr[i] = NV_Ith_S(y,i) ;
  TC_getJacTYN ( tempNmsfr, Nspec_, jactmp, (unsigned int) useJacAnl ) ;
#else
  /* species, Nvars_ = T+(Nspec_-1) */
  double sumY = 0.0 ;
  for (int i = 1; i < Nvars_; i++)
  {
    tempNmsfr[i] = NV_Ith_S(y,i) ;
    sumY += tempNmsfr[i] ;
  }
  tempNmsfr[Nvars_] = 1.0-sumY ;

  TC_getJacTYNm1 ( tempNmsfr, Nspec_, jactmp, (unsigned int) useJacAnl ) ;
#endif

  for (int j = 0 ; j < Nvars_ ; j++)
  {
    realtype *colJ = DENSE_COL(Jac,j) ;
    for ( int i = 0, iJac = j*Nvars_ ; i < Nvars_ ; i++,iJac++ )
      colJ[i] = jactmp[iJac] ;
  }

  return ;
  
}

#endif

int StiffInteg::findign(realtype t, N_Vector y, realtype *gout, void *udata)
{

  realtype temperature;
  int N = NV_LENGTH_S(y) ;
  double *udataloc = (double *)udata ;
  temperature = NV_Ith_S(y,0) ;

  gout[0] = temperature-udataloc[N*N] ;

  return(0);

}


int StiffInteg::Check_CVflag(void *flagvalue, char *funcname, int opt)
{

  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if ( ( opt == 0 ) && ( flagvalue == NULL ) ) 
  {
    std::cout<<"CVODE_ERROR: "<<std::string(funcname)
	     <<" failed - returned NULL pointer"<<std::endl<<std::flush ;
    std::terminate() ;
  }

  /* Check if flag < 0 */
  else if ( opt == 1 ) 
  {
    errflag = (int *) flagvalue;
    if ( *errflag < 0 ) 
    {
      std::cout<<"CVODE_ERROR: "<<std::string(funcname)
	       <<" failed with flag = "<<*((int *) flagvalue)<<std::endl<<std::flush ;
      std::terminate() ;
    }
  }

  return ( 0 ) ;

}
