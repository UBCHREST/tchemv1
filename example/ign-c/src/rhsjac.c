#include "ign.h"

void chemrhs(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{

  int iscal, Ns, Nv ;
  double sumY ;

  Nv = NV_LENGTH_S(y) ;

  if ( tempNmsfr == NULL ) 
    tempNmsfr = (double *) malloc( (Nv+1) * sizeof(double) ) ;
  if ( rhsvals   == NULL ) 
    rhsvals   = (double *) malloc( (Nv+1) * sizeof(double) ) ;

  /* temperature */
  tempNmsfr[0] = NV_Ith_S(y,0) ;

#ifdef ALLSPEC
  /* here Nv = 1+(Nspec) */
  Ns = Nv-1 ;
  for ( iscal = 1; iscal < Nv; iscal++ )
    tempNmsfr[iscal] = NV_Ith_S(y,iscal) ;
#else
  /* here Nvars_ = 1+(Nspec-1) */
  Ns = Nv ;
  sumY = 0.0 ;
  for ( iscal = 1; iscal < Nv; iscal++)
  {
    tempNmsfr[iscal] = NV_Ith_S(y,iscal) ;
    sumY += tempNmsfr[iscal] ;
  }
  tempNmsfr[Nv] = 1.0-sumY ;
#endif

  /* get src term from tchem library */
  TC_getSrc ( tempNmsfr, Ns+1, rhsvals ) ;

  for ( iscal = 0; iscal < Nv ; iscal++)
  {
    NV_Ith_S(ydot,iscal) = rhsvals[iscal] ;
//    printf("%d: %e\n",iscal,rhsvals[iscal]);
  }
//  exit(1);
  return ;

}

int chemrhswrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{

  chemrhs(t,y,ydot,f_data) ;
  return ( 0 );

}

#ifdef USEJAC

void chemjac(realtype t, N_Vector y, DlsMat Jac, double *jactmp) 
{

  int i, j, iJac, Nv, Ns ;
  double sumY ;
  unsigned int useJacAnl = 1 ; /* use analytic Jacobian */

  Nv = NV_LENGTH_S(y) ;
  if ( tempNmsfr == NULL ) 
    tempNmsfr = (double *) malloc( (Nv+1) * sizeof(double) ) ;

  /* temperature */
  tempNmsfr[0] = NV_Ith_S(y,0) ;

#ifdef ALLSPEC
  /* Nvars_ = T+(Nspec) */
  Ns = Nv-1 ;
  for ( i = 1; i < Nv; i++)
    tempNmsfr[i] = NV_Ith_S(y,i) ;

  TC_getJacTYN ( tempNmsfr, Ns, jactmp, useJacAnl ) ;

#else
  /* species, Nvars_ = T+(Nspec-1) */
  Ns = Nv ;
  sumY = 0.0 ;
  for ( i = 1; i < Nv; i++)
  {
    tempNmsfr[i] = NV_Ith_S(y,i) ;
    sumY += tempNmsfr[i] ;
  }
  tempNmsfr[Nv] = 1.0-sumY ;

  TC_getJacTYNm1 ( tempNmsfr, Ns, jactmp, useJacAnl ) ;
#endif

  for ( j = 0 ; j < Nv ; j++)
  {
    realtype *colJ = DENSE_COL(Jac,j) ;
    for ( i = 0, iJac = j*Nv ; i < Nv ; i++,iJac++ )
      colJ[i] = jactmp[iJac] ;
  }

  return ;
  
}

int chemjacwrapper(int N, realtype t, 
                   N_Vector y, N_Vector fy, DlsMat J, void *udata,
		   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) 
{

  double *jactmp = (double*) udata ;
  chemjac( t, y, J, jactmp) ;

  return ( 0 );

}

#endif

int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{

  int Nvars ;
  realtype temper;
  double *udataloc = (double *)user_data;

  Nvars = NV_LENGTH_S(y) ;
  temper = NV_Ith_S(y,0) ;

  gout[0] = temper-udataloc[Nvars*Nvars] ;

  return(0);
}



int Check_CVflag(void *flagvalue, char *funcname, int opt)
{

  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if ( ( opt == 0 ) && ( flagvalue == NULL ) )
  {
    printf("CVODE_ERROR: %s failed - returned NULL pointer\n",funcname) ;
    exit(1) ;
  }

  /* Check if flag < 0 */
  else if ( opt == 1 )
  {
    errflag = (int *) flagvalue;
    if ( *errflag < 0 )
    {
      printf("CVODE_ERROR: %s failed with flag = %d\n", funcname, *errflag );
      exit(1) ;
    }
  }

  return ( 0 ) ;

}
