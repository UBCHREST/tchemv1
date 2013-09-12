#ifndef IGNHSEEN
#define IGNHSEEN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

/*  CVODE headers  */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense                */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM        */
#include <sundials/sundials_types.h> /* definition of type realtype          */

#include "TC_interface.h"
#include "TC_params.h"

#define MAX(A,B) ( ((A) > (B)) ? (A) : (B) )
#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )

#define MAXSPECIN 200


static double *tempNmsfr=NULL, *rhsvals=NULL ;

/* Output files */

FILE *myfile, *myfile1, *myfile2, *myfile3 ;

void initCVODE (void **cvode_mem, double tstart, double tEnd, N_Vector y0, int Nvars,
                double *udata, realtype relT, N_Vector absT, 
                int CVmaxord, int CVmaxnumsteps ) ;

double doIgn(int getIgnDel, double tsta, double *tend, double *deltat, double deltatMax, double deltaTemp, 
             int oFreq, int *NiterMax, double *scal, double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
             FILE *myfile, FILE *myfile1, FILE *myfile2, FILE *myfile3, 
             void *cvode_mem) ;

double doIgnReinit(int getIgnDel, double tsta, double *tEnd, double *deltat, double deltatMax, double deltaTemp, 
             int oFreq, int *NiterMax, double *scal, double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
             FILE *myfile, FILE *myfile1, FILE *myfile2, FILE *myfile3, 
		   void *cvode_mem)  ;

double doIgnDtSeq(int getIgnDel, double tsta, double *tEnd, double *deltat, double deltatMax, 
                  double deltaTemp, int oFreq, int *NiterMax, double *scal, 
                  int getDeltatSeq, double *deltatSeq,
                  double CVsmall, N_Vector y0, realtype relT, N_Vector absT, 
                  FILE *myfile, FILE *myfile1, FILE *myfile2, FILE *myfile3, 
                  void *cvode_mem) ;

void Output(int iflag, int iter, double t, double deltat, N_Vector y0, double *scal) ;

void setup(int *NiterMax, int *oFreq, double *Tini, double *Temp_id, 
           double *deltat, double *deltatMax, double *tEnd, double *deltaTemp, 
           double *pressure, double *pfac, unsigned int *withTab, int *getIgnDel, 
           double *CVrelt, double *CVsmall, int *CVmaxord, int *CVmaxnumsteps, 
           char *SpecName, double *SpecMsFr, int *specinno, 
           int *getDeltatSeq, double **deltatSeq, 
           char *mechfile, char *thermofile) ;

void chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) ;
int chemrhswrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data) ;
#ifdef USEJAC
void chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp) ; 
int chemjacwrapper(int N, realtype t, 
                   N_Vector y, N_Vector fy, DlsMat J, void *udata,
		   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) ;
#endif
int g(realtype t, N_Vector y, realtype *gout, void *user_data) ;
int Check_CVflag(void *flagvalue, char *funcname, int opt) ;

#endif
