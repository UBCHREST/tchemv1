#ifndef TCkmodintHSeen
#define TCkmodintHSeen

#include "copyright.h"

/* Length of filenames */
#define lenfile 100

/* Lengths of character strings */
#define lenstr01 1000
#define lenstr02 100
#define lenstr03 10
#define lenread  256

/* Reallocation factors for element, species, and reaction lists */
#define nelemalloc 1
#define nspecalloc 1
#define nreacalloc 5

/* -------------------------------------------------------------------------
                           Structure definitions
   ------------------------------------------------------------------------- */

/**  
 * \brief Entry in the table of periodic elements 
 */
typedef struct 
{
  char    name[LENGTHOFELEMNAME] ;
  double  mass                   ;
} elemtable ;

/**  
 * \brief Element data 
 */
typedef struct
{
  char   name[LENGTHOFELEMNAME] ;
  int    hasmass                ;
  double mass                   ;
} element ;

/** 
 * \brief Species data 
 */
typedef struct
{
  char   name[LENGTHOFSPECNAME]          ;
  int    hasthermo                       ;
  int    hasmass                         ;
  double mass                            ;
  int    charge                          ;
  int    phase                           ;
  int    numofelem                       ;
  int    elemindx[NUMBEROFELEMINSPEC]    ;
  int    elemcontent[NUMBEROFELEMINSPEC] ;
  double nasapoltemp[3]                  ;
  double nasapolcoefs[14]                ;
} species ;

/** 
 * \brief Reaction data 
 */
typedef struct
{
  int isdup     ; /* Duplicate reactions */
  int isreal    ; /* Real stoichiometric coefficients */
  int isrev     ; /* Is reversible */
  int isfall    ; /* Is pressure dependent */
  int specfall  ; /* Type of concetration used for pressure dependency */
  int isthrdb   ; /* Uses third body */
  int nthrdb    ; /* Number of third body coefficients */
  int iswl      ; /* Is photon activation */
  int isbal     ; /* Is balanced */
  int iscomp    ; /* Is complete */
  int inreac    ; /* Number of reactants */
  int inprod    ; /* Number of products */
  int ismome    ; /* Is a momentum-transfer collision */
  int isxsmi    ; /* Is a ion momentum-transfer collision cross-section */
  int isford    ; /* Arbitrary reaction order - forward coeffs */
  int isrord    ; /* Arbitrary reaction order - reverse coeffs */

  int islowset  ; /*  LOW  parameters set */
  int ishighset ; /*  HIGH parameters set */
  int istroeset ; /*  TROE parameters set */
  int issriset  ; /*  SRI  parameters set */
  int isrevset  ; /*  REV  parameters set */
  int isltset   ; /*  LT   parameters set */
  int isrltset  ; /*  RLT  parameters set */
  int ishvset   ; /*  HV   wavelength set */
  int istdepset ; /*  TDEP species    set */
  int isexciset ; /*  EXCI parameter  set */
  int isjanset  ; /*  JAN  parameters set */
  int isfit1set ; /*  FIT1 parameters set */
  
  int    spec [2*NSPECREACMAX] ; /* List of species */
  int    nuki [2*NSPECREACMAX] ; /* Stoichiometric coefficients (integer) */
  double rnuki[2*NSPECREACMAX] ; /* Stoichiometric coefficients (real) */

  double arhenfor[3] ; /* Arrhenius parameters (forward) */
  double arhenrev[3] ; /* Arrhenius parameters (reverse) */
  
  int    ithrdb[NTHRDBMAX] ; /* Indices of third-body species */
  double rthrdb[NTHRDBMAX] ; /* Enhanced efficiencies of third-body species */
  
  double fallpar[8] ; /* Parameters for pressure-dependent reactions */
  double ltpar  [2] ; /* Landau-Teller parameters */
  double rltpar [2] ; /* Reverse Landau-Teller parameters */
  
  double hvpar ; /* Radiation wavelength */
  
  char aunits[4] ; /* pre-exponential factor units */
  char eunits[4] ; /* activation energy units */
  
  int  tdeppar ; /* Species # for temperature dependence */

  double excipar ; /* EXCI (energy loss) parameter */

  double optfit[9] ; /* Optional rate-fit parameters */

  /* Arbitrary reaction order */
  int    arbspec[4*NSPECREACMAX] ; /* List of species for arbitrary reac order */
  double arbnuki[4*NSPECREACMAX] ; /* and list of coefficients */
		
} reaction ;

/* -------------------------------------------------------------------------
                           Elements' functions
   ------------------------------------------------------------------------- */

void setperiodictable ( elemtable *periodictable,int *Natoms,int iflag) ;

void checkeleminlist (char *elemname,element *listelem,int *Nelem, int *ipos) ;

int  getelements (char *linein, char *singleword,element **listelemaddr,
                  int *Nelem, int *Nelemmax, int *iread, int *ierror) ;

void resetelemdata (element *currentelem) ;

int  setelementmass (element *listelem,int *Nelem, elemtable *periodictable,int *Natoms, 
                     int *ierror) ;

/* -------------------------------------------------------------------------
                           Species' functions
   ------------------------------------------------------------------------- */

int getspecies ( char *linein,char *singleword,species **listspecaddr,
                 int *Nspec,int *Nspecmax,int *iread,int *ierror ) ;

void resetspecdata ( species *currentspec ) ;

int setspecmass ( element *listelem,int *Nelem,species *listspec,int *Nspec, 
                  int *ierror ) ;

/* -------------------------------------------------------------------------
                           Thermo functions
   ------------------------------------------------------------------------- */
int checkthermo(species *listspec,int *Nspec) ;

void checkspecinlist(char *specname,species *listspec,int *Nspec, int *ipos) ;

int getthermo(char *linein,char *singleword,FILE *mechin,FILE *thermoin,
              element *listelem,int *Nelem,species *listspec,int *Nspec,double *Tglobal,
              int *ithermo,int *iread,int *ierror) ;

/* -------------------------------------------------------------------------
                           Reactions' functions
   ------------------------------------------------------------------------- */
void resetreacdata(reaction *currentreac, char *aunits, char *eunits) ;

void checkunits (char *linein,char *singleword,char *aunits,char *eunits) ;

int getreacline (char *linein,char *singleword, species *listspec,int *Nspec,
                 reaction *listreac,int *Nreac, int *ierror) ;

int getreacauxl (char *linein,char *singleword, species *listspec,int *Nspec,
                 reaction *listreac, int *Nreac, int *ierror) ;

int getreactions (char *linein,char *singleword, species *listspec,int *Nspec,
                  reaction *listreac, int *Nreac, char *aunits,char *eunits,
                  int *ierror) ;

int rescalereac (reaction *listreac,int *Nreac ) ;

int verifyreac (element *listelem,int *Nelem, species *listspec,int *Nspec,
                reaction *listreac,int *Nreac,int *ierror) ;

int kmodsum(element *listelem,int *Nelem, species *listspec,int *Nspec,
                 reaction *listreac,int *Nreac,
                 int *nIonEspec,int *electrIndx,int *nIonSpec,int *maxSpecInReac,int *maxTbInReac,int *maxOrdPar,
                 int *nFallPar,int *maxTpRange,int *nLtReac,int *nRltReac,int *nFallReac,
                 int *nThbReac,int *nRevReac,int *nHvReac,int *nTdepReac,int *nJanReac,int *nFit1Reac,
                 int *nExciReac,int *nMomeReac,int *nXsmiReac,int *nRealNuReac,int *nOrdReac,
                 int *nNASAinter,int *nCpCoef,int *nNASAfit,int *nArhPar,
                 int *nLtPar,int *nJanPar,int *nFit1Par) ;
/* -------------------------------------------------------------------------
                           I/O functions
   ------------------------------------------------------------------------- */
int out_formatted   (element *listelem,int *Nelem, species *listspec,int *Nspec,
                    reaction *listreac,int *Nreac, char *aunits,char *eunits,
                    FILE *fileascii) ;
int out_unformatted (element *listelem,int *Nelem, species *listspec,int *Nspec,
                    reaction *listreac,int *Nreac, char *aunits,char *eunits,
                    FILE *filelist,int *ierror) ;
int out_mathem      (element *listelem,int *Nelem, species *listspec,int *Nspec,
		    reaction *listreac,int *Nreac, char *aunits,char *eunits)  ;


/* -------------------------------------------------------------------------
                           Error functions
   ------------------------------------------------------------------------- */
void errormsg (int ierror) ;

/* -------------------------------------------------------------------------
                           Character string functions
   ------------------------------------------------------------------------- */
int elimleads  ( char *linein ) ;
int elimends   ( char *linein ) ;
int elimspaces ( char *linein ) ;
int elimcomm   ( char *linein ) ;
int tab2space  ( char *linein ) ;
int extractWordLeft  ( char *linein,char *oneword ) ;
int extractWordRight ( char *linein,char *oneword ) ;
int extractWordLeftauxline ( char *linein,char *oneword,char *twoword,
                             int *inum,int *ierror ) ; 
int extractWordLeftNoslash ( char *linein,char *oneword ) ;

int extractdouble ( char *wordval,double *dvalues,int *inum,int *ierror ) ;

void wordtoupper ( char *linein,char *oneword,int Npos) ;

void cleancharstring ( char *linein,int *len1) ;

int charfixespc(char *singleword,int *len1) ;

int checkstrnum ( char *singleword,int *len1,int *ierror) ;

int findnonnum ( char *specname,int *ipos) ;

#endif
