#ifndef StiffIntegHSeen
#define StiffIntegHSeen

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include "math.h"

#include "TC_interface.h"

/*  CVODE headers  */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define MAX(A,B) ( ((A) > (B)) ? (A) : (B) )
#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )

class StiffInteg
{
    
 public:

  StiffInteg( double *scal, int Nspec, double dTmax )  {
	
    isInit_    = false ;

    scal_   = scal  ;
    Nspec_  = Nspec ;
    dTmax_  = dTmax ;

#ifdef ALLSPEC
    Nvars_  = Nspec+1 ; /* T+(Nspec)   */
#else
    Nvars_  = Nspec   ; /* T+(Nspec-1) */
#endif

    /* add class to the list of objects of the same type */
    if (! StiffInteg::omap_) StiffInteg::omap_ = new StiffInteg::OMap_t;
    
    this->my_index_ = StiffInteg::next_index_++;
    (*StiffInteg::omap_)[my_index_] = this ;
    
    /* local arrays */
    tempNmsfr = 0 ;
    rhsvals   = 0 ;
    
    std::cout<<"StiffInteg : contructor : index : "<<my_index_
	     <<std::endl<<std::flush ;
    
  }

  ~StiffInteg() {

    /* Reset pointer */
    scal_    = 0 ; 
    if (tempNmsfr != 0) delete []tempNmsfr ;
    if (rhsvals   != 0) delete []rhsvals   ;

    /* Remove this class from the list of objects */
    OMap_t::iterator me = StiffInteg::omap_->find(this->my_index_);
    if (StiffInteg::omap_ && (me != StiffInteg::omap_->end())) 
      StiffInteg::omap_->erase(me);
    
    std::cout<<"StiffInteg : destructor : index : "<<my_index_
	     <<std::endl<<std::flush ;
    
  }

  int setRelT(double relT) {
    if ( !isInit_ ) init() ;
    CVrelt_ = relT ;
    return ( 0 ) ;
  } ;
  
  int setSmall(double small) {
    if ( !isInit_ ) init() ;
    CVsmall_ = small ;
    return ( 0 ) ;
  } ;

  int setTempThresh( double Tid ) {
    if ( !isInit_ ) init() ;
    Temp_id = Tid;
    return ( 0 ) ;
  }
  
  int compute(double tEnd,double *deltat, double deltatMax, int NiterMax, int oFreq) ;
  
 private:

  bool isInit_ ;
  int init() ;
  
  double *scal_ ;

  double CVrelt_, CVsmall_ ;
  int CVmaxord_, CVmaxnumsteps_ ;
  bool useJacAnl ;
  
  int Nvars_, Nspec_ ;
  double *tempNmsfr, *rhsvals, dTmax_, Temp_id ;

  /* Related to map of StiffInteg objects */
  typedef std::map<int, StiffInteg*> OMap_t ;
  static  int next_index_ ;
  static  OMap_t *omap_ ;
  int     my_index_ ;
  
  static int chemrhswrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data)
  {
    int indxobj = ((int*) f_data)[0] ;
    
    OMap_t::iterator it = omap_->find((int) indxobj);
    
    if (it == omap_->end()) 
    {
      std::cout<<"chemrhswrapper(): the callback object is not a valid entry in the map"
	       <<std::endl<<std::flush ;
      std::terminate() ;
    }
    
    /* Perform callback to the member function of the proper StiffInteg class */
    (it->second)->chemrhs(t,y,ydot,f_data) ;
    
    return ( 0 );
  }
  void chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) ;

#ifdef USEJAC
  static int chemjacwrapper(int N, realtype t, 
                            N_Vector y, N_Vector fy, DlsMat J, void *udata,
			    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3     ) 
  {
    int indxobj = ((int*) udata)[0] ;

    OMap_t::iterator it = omap_->find((int) indxobj);
    
    if (it == omap_->end()) 
    {
      std::cout<<"chemjacwrapper(): the callback object is not a valid entry in the map"
	       <<std::endl<<std::flush ;
      std::terminate() ;
    }
	
    /* Perform callback to the member function of the proper StiffInteg class */
    double *jactmp = &((double*) udata)[1] ;
    (it->second)->chemjac(t,y,J,jactmp) ;

    return ( 0 );
  }
  void chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp) ; 
#endif

  static int findignwrap(realtype t, N_Vector y, realtype *gout, void *udata)
  {

     int indxobj = ((int*) udata)[0] ;

     OMap_t::iterator it = omap_->find((int) indxobj);
    
      if (it == omap_->end()) 
      {
        std::cout<<"chemjacwrapper(): the callback object is not a valid entry in the map"
  	       <<std::endl<<std::flush ;
        std::terminate() ;
      }
	
      /* Perform callback to the member function of the proper StiffInteg class */
      double *utmp = &((double*) udata)[1] ;
      (it->second)->findign(t,y,gout,utmp) ;

      return ( 0 );
  }

  int findign(realtype t, N_Vector y, realtype *gout, void *udata) ;

  /* Error messages from cvode */
  int Check_CVflag(void *flagvalue, char *funcname, int opt) ;

} ;

#endif
