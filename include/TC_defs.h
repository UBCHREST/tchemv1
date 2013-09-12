#ifndef _TCdefsHSeen_
#define _TCdefsHSeen_

/*! #
 * \file TC_defs.h
 * \brief Definitions of variables names used by the library.
*/ 
#include "copyright.h"
#include <stdio.h>
#include "TC_params.h"

/* -------------------------------------------- 
   numbers, lengths, etc.
   -------------------------------------------- */
/** \var static int TC_maxSpecInReac_
 *  \ingroup maxpar
 *  \brief Maximum number of species in a reaction */
static int TC_maxSpecInReac_ ;
/** \var static int TC_maxTbInReac_
 *  \ingroup maxpar
    \brief Max # of third-body efficiencies in a reaction */
static int TC_maxTbInReac_   ;
/** \var static int TC_nNASAinter_
 *  \ingroup maxpar
    \brief # of temperature regions for thermo fits */
static int TC_nNASAinter_    ;
/** \var static int TC_nCpCoef_
 *  \ingroup maxpar
    \brief # of polynomial coefficients for thermo fits */
static int TC_nCpCoef_       ;
/** \var static int TC_nArhPar_
 *  \ingroup maxpar
    \brief # of Arrhenius parameters */
static int TC_nArhPar_       ;
/** \var static int TC_nLtPar_
 *  \ingroup maxpar
    \brief # of parameters for Landau-Teller reactions */
static int TC_nLtPar_        ;
/** \var static int TC_nFallPar_
 *  \ingroup maxpar
    \brief # of parameters for pressure-dependent reactions */
static int TC_nFallPar_      ;
/** \var static int TC_nJanPar_
 *  \ingroup maxpar
    \brief # of parameters for Jannev-Langer fits (JAN) */
static int TC_nJanPar_       ;
/** \var static int TC_maxOrdPar_
 *  \ingroup maxpar
    \brief # of parameters for arbitrary order reactions */
static int TC_maxOrdPar_     ;
/** \var static int TC_nFit1Par_
 *  \ingroup maxpar
    \brief # of parameters for FIT1 fits */
static int TC_nFit1Par_      ;

/** \var static int TC_Nvars_
 *  \ingroup nospec
    \brief # of variables = no. of species + 1 */
static int TC_Nvars_         ; 
/** \var static int TC_Nvjac_
 *  \ingroup nospec
    \brief # of lines/cols in the Jacobian = no. of species + 3 */
static int TC_Nvjac_         ;
/** \var static int TC_Nelem_
 *  \ingroup nospec
    \brief # of chemical elements */
static int TC_Nelem_         ; 
/** \var static int TC_Nspec_
 *  \ingroup nospec
    \brief # of species */
static int TC_Nspec_         ;
/** \var static int TC_nIonSpec_
 *  \ingroup nospec
    \brief # of ion species */
static int TC_nIonSpec_      ; 
/** \var static int TC_electrIndx_
 *  \ingroup nospec
    \brief Index of the electron species */
static int TC_electrIndx_    ;
/** \var static int TC_nIonEspec_
 *  \ingroup nospec
    \brief # of ion species excluding the electron species */
static int TC_nIonEspec_     ;

/** \var static int TC_Nreac_
 *  \ingroup noreac
    \brief # of reactions */
static int TC_Nreac_         ;
/** \var static int TC_nRevReac_
 *  \ingroup noreac
    \brief # of reactions with REV given */
static int TC_nRevReac_      ;
/** \var static int TC_nFallReac_
 *  \ingroup noreac
    \brief # of pressure-dependent reactions */
static int TC_nFallReac_     ;
/** \var static int TC_nThbReac_
 *  \ingroup noreac
    \brief # of reactions using third-body efficiencies */
static int TC_nThbReac_      ;
/** \var static int TC_nLtReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions */
static int TC_nLtReac_       ;
/** \var static int TC_nRltReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions with RLT given */
static int TC_nRltReac_      ;
/** \var static int TC_nHvReac_
 *  \ingroup noreac
    \brief # of reactions with HV */
static int TC_nHvReac_       ;
/** \var static int TC_nJanReac_
 *  \ingroup noreac
    \brief # of reactions with JAN fits */
static int TC_nJanReac_      ;
/** \var static int TC_nFit1Reac_
 *  \ingroup noreac
    \brief # of reactions with FIT1 fits */
static int TC_nFit1Reac_     ;
/** \var static int TC_nExciReac_
 *  \ingroup noreac
    \brief # of reactions with EXCI */
static int TC_nExciReac_     ;
/** \var static int TC_nMomeReac_
 *  \ingroup noreac
    \brief # of reactions with MOME */
static int TC_nMomeReac_     ;
/** \var static int TC_nXsmiReac_
 *  \ingroup noreac
    \brief # of reactions with XSMI */
static int TC_nXsmiReac_     ;
/** \var static int TC_nTdepReac_
 *  \ingroup noreac
    \brief # of reactions with TDEP */
static int TC_nTdepReac_     ;
/** \var static int TC_nRealNuReac_
 *  \ingroup noreac
    \brief # of reactions with non-int stoichiometric coefficients */
static int TC_nRealNuReac_   ;
/** \var static int TC_nOrdReac_
 *  \ingroup noreac
    \brief # of reactions with arbitrary order */
static int TC_nOrdReac_      ;

/* --------------------------------------------
      elements & species -> names, masses, and element count 
   -------------------------------------------- */
/** \var static char TC_sNames_
    \brief species names, name of species i stored at LENGTHOFSPECNAME*i */
static char   *TC_sNames_    ;
/** \var static char TC_eNames_
    \brief species names, name of species i stored at LENGTHOFELEMNAME*i */
static char   *TC_eNames_    ; 
/** \var static double TC_sMass_
    \brief array of species molar masses */
static double *TC_sMass_     ;
/** \var static double TC_eMass_
    \brief array of element molar masses */
static double *TC_eMass_     ; 
/** \var static int TC_elemcount_
    \brief no. of atoms of element j in each species i at (i*TC_Nelem_+j)*/
static int  *TC_elemcount_ ;
  
/* --------------------------------------------
     species charges, number of temp fits, phase 
   -------------------------------------------- */
/** \var static int TC_sCharge_
    \brief species electrical charges */
static int *TC_sCharge_ ;
/** \var static int TC_sTfit_
    \brief no. of temperature fits for thermodynamic properties */
static int *TC_sTfit_   ;
/** \var static int TC_sPhase_
    \brief species phase id */
static int *TC_sPhase_  ; 

/* temperature limits for NASA polynomials and their coefficients */
static double *TC_Tlo_,*TC_Tmi_,*TC_Thi_ ; 
static double *TC_cppol_ ;
  
/* list of non-electron ion species */
static int *TC_sNion_ ; 
  
/* Reaction data */
/* is reaction reversible ? */
static int *TC_isRev_ ; 

/* no. of reac+prod, no. of reac only, no. of prod only, stoichiom.coef. indicator */
static int *TC_reacNrp_, *TC_reacNreac_, *TC_reacNprod_, *TC_reacScoef_; 

/* stoichiometric coeffs and reactants and product indices */
static int *TC_reacNuki_, *TC_reacSidx_ ;
static double *TC_reacNukiDbl_ ;

/* real stoichiometric coeffs */
static int    *TC_reacRnu_ ;
static double *TC_reacRealNuki_ ;

/* list of reactions with given reverse Arrhenius parameters */
static int *TC_reacRev_ ;

/* Arrhenius parameters, forward and reverse */
static double *TC_reacArhenFor_, *TC_reacArhenRev_ ;

/* Pressure dependent reactions */
static int *TC_reacPfal_, *TC_reacPtype_, *TC_reacPlohi_, *TC_reacPspec_ ; 
static double *TC_reacPpar_ ;

/* Third-body data */
static int *TC_reacTbdy_, *TC_reacTbno_, *TC_specTbdIdx_ ;
static double *TC_specTbdEff_ ;

/* Reactions with arbitrary orders  */
static int    *TC_reacAOrd_, *TC_specAOidx_ ;
static double *TC_specAOval_ ;

/* Radiation wavelength data */
static int    *TC_reacHvIdx_ ;
static double *TC_reacHvPar_ ;

/* Other work arrays */
static double *TC_kfor, *TC_krev, *TC_kforP, *TC_krevP, *TC_ropFor, *TC_ropRev, *TC_rop;
static double *TC_cpks, *TC_hks, *TC_gk, *TC_gkp ;
static double *TC_PrDer, *TC_Crnd, *TC_CrndDer, *TC_dFfac, *TC_omg, *TC_omgP, *TC_jacFull ;
static double *TC_Xconc,*TC_scalIn, *TC_Mconc ;
static int    *TC_sigNu, *TC_NuIJ ;
static double *TC_sigRealNu, *TC_RealNuIJ ;
static double *qfr ;

/* variables */
static double TC_reltol,TC_abstol ;
static double TC_rhoref_, TC_pref_, TC_Tref_, TC_Wref_, TC_Daref_, TC_omgref_, TC_cpref_, TC_href_ ;
static double TC_timref_ ;
static int    TC_isInit_, TC_tab_, TC_nonDim_, TC_RVset_ ;

/* tables */
static int    TC_Ntab_ ;
static double TC_delT_, TC_odelT_ ;
static double *TC_cptab,   *TC_cpPtab,  *TC_htab,     *TC_gktab, *TC_gkPtab ; 
static double *TC_kfortab, *TC_krevtab, *TC_kforPtab, *TC_krevPtab ;

/* thermodynamic pressure */
static double TC_pressure_, TC_prescgs_ ;

/* universal gas constant */
static double TC_Runiv_, TC_Rcal_, TC_Rcgs_ ;

/* Backup Arrhenius factors */
static int  TC_ArhenForChg_, TC_ArhenRevChg_ ;
static double *TC_reacArhenForSave_, *TC_reacArhenRevSave_ ;

/* kc coefficients - depend on pressure only */
static double *TC_kc_coeff ;

/* Backups for reactional removal */
static int TC_initRemoved ;
static int TC_NreacBackup_, TC_NrevBackup_, TC_NfalBackup_, TC_NthbBackup_ ;
static int TC_nRealNuReacBackup_, TC_nOrdReacBackup_ ;
static int *TC_isRevBackup_, *TC_reacNrpBackup_, *TC_reacNreacBackup_, *TC_reacNprodBackup_ ;
static int *TC_reacNukiBackup_, *TC_reacSidxBackup_, *TC_reacScoefBackup_ ;
static double *TC_reacArhenForBackup_, *TC_reacArhenRevBackup_, *TC_reacNukiDblBackup_ ;
static int *TC_reacPfalBackup_, *TC_reacPtypeBackup_, *TC_reacPlohiBackup_, *TC_reacPspecBackup_ ;
static int *TC_reacTbdyBackup_, *TC_reacTbnoBackup_, *TC_specTbdIdxBackup_ ;
static double *TC_reacPparBackup_, *TC_specTbdEffBackup_;
static int *TC_reacRnuBackup_, *TC_reacAOrdBackup_, *TC_specAOidxBackup_ ;
static double *TC_reacRealNukiBackup_, *TC_sigRealNuBackup, *TC_RealNuIJBackup;
static double *TC_specAOvalBackup_, *TC_kc_coeffBackup ;
static int *TC_sigNuBackup, *TC_NuIJBackup ;

/* functions */
int  TC_kmodint_(char *mechfile,int *lmech,char *thermofile,int *lthrm) ;

/*
             _____ ____     _       _ _         
            |_   _/ ___|   (_)_ __ (_) |_   ___ 
              | || |       | | '_ \| | __| / __|
              | || |___    | | | | | | |_ | (__ 
              |_| \____|___|_|_| |_|_|\__(_)___|
                      |_____|                   
      
  
*/

int TC_makeSpace() ;

int TC_createTables( double delT ) ;

void TC_errorMSG(int msgID, char const *func, int var1, int var2) ;
void TC_errorINI(FILE *errfile, char const *msg) ;

/*
        _____ ____          _   _ _           
       |_   _/ ___|   _   _| |_(_) |___   ___ 
         | || |      | | | | __| | / __| / __|
         | || |___   | |_| | |_| | \__ \| (__ 
         |_| \____|___\__,_|\__|_|_|___(_)___|
                 |_____|                      
 */
static double fastIntPow(double val, int exponent);

/*
          _____ ____                          
         |_   _/ ___|    _ __ _ __        ___ 
           | || |       | '__| '__|      / __|
           | || |___    | |  | |     _  | (__ 
           |_| \____|___|_|  |_|    (_)  \___|
                   |_____|                    
*/


/* molar rate of progress [kmol/(m3.s)] */
int TC_getRopsLocal(double *scal) ;

/* molar reaction rates [kmol/(m3.s)] */
int TC_getReacRates(double *scal, int Nvars, double *omega) ;

int TC_getgk(double t1,double t_1, double tln) ;
int TC_getgkFcn(double t1,double t_1, double tln) ;
int TC_getgkTab(double t1) ;

int TC_getgkp(double t1,double t_1, double tln) ;
int TC_getgkpFcn(double t1,double t_1, double tln) ;
int TC_getgkpTab(double t1) ;

double TC_getSumNuGk(int i, double *gkLoc) ;
double TC_getSumRealNuGk(int i,int ir, double *gkLoc) ;

int TC_get3rdBdyConc(double *concX, double *concM) ;

int TC_getkForRev(double t1, double t_1, double tln) ;
int TC_getkForRevFcn(double t_1, double tln) ;
int TC_getkForRevTab(double t1) ;

int TC_getkForRevP(double t1, double t_1) ;
int TC_getkForRevPFcn(double t_1) ;
int TC_getkForRevPTab(double t1) ;

int TC_getRateofProg(double *concX) ;
int TC_getRateofProgDer(double *concX, int ireac, int ispec, double *qfr);
int TC_getCrnd(double t1,double t_1,double tln,double *concX, double *concM) ;
int TC_getCrndDer(int ireac, int *itbdy, int *ipfal, 
    double t1,double t_1,double tln,double *concX, double *concM) ;

/*
      _____ ____    _   _                                          
     |_   _/ ___|  | |_| |__   ___ _ __ _ __ ___   ___         ___ 
       | || |      | __| '_ \ / _ \ '__| '_ ` _ \ / _ \       / __|
       | || |___   | |_| | | |  __/ |  | | | | | | (_) |  _  | (__ 
       |_| \____|___\__|_| |_|\___|_|  |_| |_| |_|\___/  (_)  \___|
               |_____|                                             
*/
int TC_getCpSpecMsFcn(double t ,double *cpi) ;
int TC_getCpSpecMsTab(double t1,double *cpi) ;
int TC_getCpSpecMlFcn(double t ,double *cpi) ;
int TC_getCpSpecMlTab(double t1,double *cpi) ;

int TC_getCpSpecMs1Fcn(double t, int i, double *cpi) ;
int TC_getCpSpecMl1Fcn(double t, int i, double *cpi) ;

int TC_getCpMixMsP(double *scal,int Nvars,double *cpmix) ;
int TC_getCpSpecMsP(double t,int Nspec,double *cpi) ;

int TC_getCpSpecMsPFcn(double t ,double *cpi) ;
int TC_getCpSpecMsPtab(double t1,double *cpi) ;

int TC_getCpSpecMs1PFcn(double t, int i, double *cpi) ;

int TC_getHspecMsFcn(double t ,double *hi) ;
int TC_getHspecMsTab(double t1,double *hi) ;
int TC_getHspecMlFcn(double t ,double *hi) ;
int TC_getHspecMlTab(double t1,double *hi) ;

/* other defs */
#define TMAX  3500.0
#define TMIN  290.0
#define DTMID 10.0 /* buffer around middle temperature for NASA polynomials */

#endif

