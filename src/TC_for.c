/*
   TC_chg.c
 */
int tcchgarhenfor_( int *ireac, int *ipos, double *newval )
{
  int ans = 0 ;
  ans =  TC_chgArhenFor( *ireac, *ipos, *newval ) ;
  return ( ans ) ;
}

int tcchgarhenforback_( int *ireac, int *ipos)
{
  int ans = 0 ;
  ans =  TC_chgArhenForBack( *ireac, *ipos) ;
  return ( ans ) ;
}

int tcchgarhenrev_( int *ireac, int *ipos, double *newval )
{
  int ans = 0 ;
  ans =  TC_chgArhenRev( *ireac, *ipos, *newval ) ;
  return ( ans ) ;
}

int tcchgarhenrevback_( int *ireac, int *ipos)
{
  int ans = 0 ;
  ans =  TC_chgArhenRevBack( *ireac, *ipos) ;
  return ( ans ) ;
}


/*
   TC_edit.c
 */
void tcreset_()
{
  TC_reset() ;
  return ;
}
/*
   TC_init.c
 */
int tcinitchem_(char *mechfile,int *lmech,char *thermofile,int *lthrm,
                int *tab, double *delT)
{
  int ans = 0 ;
  /* add a NULL at the ends of mechfile and thermofile */
  memset(&mechfile  [*lmech], 0, 1) ;
  memset(&thermofile[*lthrm], 0, 1) ;
  ans = TC_initChem(mechfile,thermofile,*tab,*delT);
  return ( ans ) ;
}

void tcsetrefval_(double *rhoref, double *pref, double *Tref, double *Wref, 
                  double *Daref, double *omgref, double *cpref, double *href, 
                  double *timref)
{
  TC_setRefVal(*rhoref, *pref, *Tref, *Wref, *Daref, 
	       *omgref, *cpref, *href, *timref);
  return ;
}

void tcsetnondim_() 
{
  TC_setNonDim()  ;
  return ;
}

void tcsetdim_() 
{
  TC_setDim()  ;
  return ;
}

void tcsetthermopres_(double *pressure)
{
  TC_setThermoPres(*pressure) ;
  return ;
}

/*
   TC_mlms.c
 */
int tcdndgetms2cc_(double *scal,int *Nvars,double *concX)
{
  int ans = 0 ;
  ans = TCDND_getMs2Cc( scal, (*Nvars), concX) ;
  return ( ans ) ;
}

int tcgetms2cc_(double *scal,int *Nvars,double *concX)
{
  int ans = 0 ;
  ans = TC_getMs2Cc( scal, (*Nvars), concX) ;
  return ( ans ) ;
}

int tcdndgetml2ms_(double *Xspec,int *Nspec,double *Yspec)
{
  int ans = 0 ;
  ans = TCDND_getMl2Ms( Xspec, (*Nspec), Yspec) ;
  return ( ans ) ;
}

int tcgetml2ms_(double *Xspec,int *Nspec,double *Yspec)
{
  int ans = 0 ;
  ans = TC_getMl2Ms( Xspec, (*Nspec), Yspec) ;
  return ( ans ) ;
}

int tcdndgetms2ml_(double *Yspec,int *Nspec,double *Xspec)
{
  int ans = 0 ;
  ans = TCDND_getMs2Ml( Yspec, (*Nspec), Xspec) ;
  return ( ans ) ;
}

int tcgetms2ml_(double *Yspec,int *Nspec,double *Xspec)
{
  int ans = 0 ;
  ans = TC_getMs2Ml( Yspec, (*Nspec), Xspec) ;
  return ( ans ) ;
}

int tcdndgetms2wmix_(double *Yspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TCDND_getMs2Wmix(Yspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcgetms2wmix_(double *Yspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TC_getMs2Wmix(Yspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcdndgetml2wmix_(double *Xspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TCDND_getMl2Wmix(Xspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcgetml2wmix_(double *Xspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TC_getMl2Wmix(Xspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

/*
   TC_rr.c
 */
int tcgetnreac_() { return ( TC_getNreac() ) ; }   

int tcgetstoicoef_( int *Nspec, int *Nreac, double *stoicoef )
{
  int ans = 0 ;
  ans = TC_getStoiCoef( *Nspec, *Nreac, stoicoef ) ;
  return ( ans );
}

int tcgetstoicoefreac_( int *Nspec, int *Nreac, int *ireac, int *idx, 
                       double *stoicoef )
{
  int ans = 0 ;
  ans = TC_getStoiCoefReac( *Nspec, *Nreac, *ireac, *idx, stoicoef ) ;
  return ( ans );
}

int tcgetarhenfor_( int *ireac, int *ipos, double *val )
{
  int ans = 0 ;
  ans = TC_getArhenFor( *ireac, *ipos, val ) ;
  return ( ans );
}

int tcgetarhenrev_( int *ireac, int *ipos, double *val )
{
  int ans = 0 ;
  ans = TC_getArhenRev( *ireac, *ipos, val ) ;
  return ( ans );
}

int tcdndgetty2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTY2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetty2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTY2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgetty2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTY2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetty2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTY2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgettxc2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTXC2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgettxc2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTXC2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgettxc2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTXC2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgettxc2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTXC2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetrops_(double *scal, int *Nvars, double *datarop) 
{
  int ans = 0 ;
  ans = TC_getRops( scal, *Nvars, datarop) ;
  return ( ans );
}

int tcgetrfrb_(double *scal, int *Nvars, double *dataRfrb)
{
  int ans = 0 ;
  ans = TC_getRfrb(scal, *Nvars, dataRfrb) ;
  return ( ans );
}
/*
   TC_spec.c
 */
int tcgetnspec_() { return ( TC_getNspec() ) ; }   

int tcgetnelem_() { return ( TC_getNelem() ) ; }   

int tcgetnvars_() { return ( TC_getNvars() ) ; }   

int tcgetsnames_(int *Nspec,char *snames)
{
  int i, ans = 0 ;
  ans = TC_getSnames( (*Nspec), snames);
  for ( i=0; i<LENGTHOFSPECNAME * TC_Nspec_; i++)
    if ( snames[i] == 0 ) snames[i] = 32;
  return ( ans ) ;
}

int tcgetsnamelen_() { return ( TC_getSnameLen() ) ; } 

int tcgetspos_(const char *sname, const int *slen)
{
  int ans = 0 ;
  ans = TC_getSpos( sname, *slen) ;
  return ( ans ) ;
}
 
int tcgetsmass_(int *Nspec,double *Wi)
{
  int ans = 0 ;
  ans = TC_getSmass( *Nspec, Wi) ;
  return ( ans ) ;
}
/*
   TC_src.c
 */
int tcdndgetsrc_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TCDND_getSrc(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsrc_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TC_getSrc(scal,*Nvars,omega) ;
  return ( ans );
}

int tcdndgetsrccons_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TCDND_getSrcCons(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsrccons_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TC_getSrcCons(scal,*Nvars,omega) ;
  return ( ans );
}

int tcdndgetjactynm1anl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TCDND_getJacTYNm1anl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjactynm1anl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacTYNm1anl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcdndgetjactynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TCDND_getJacTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjactynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcdndgetjactynm1_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TCDND_getJacTYNm1(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjactynm1_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacTYNm1(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcdndgetjactyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TCDND_getJacTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjactyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjacrptyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacRPTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjacrptynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacRPTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjacrptynnum_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacRPTYNnum(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

/*
   TC_thermo.c
 */
int tcdndgetrhomixms_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TCDND_getRhoMixMs( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

int tcgetrhomixms_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TC_getRhoMixMs( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

int tcdndgetrhomixml_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TCDND_getRhoMixMl( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

int tcgetrhomixml_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TC_getRhoMixMl( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

int tcdndgettmixms_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TCDND_getTmixMs( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

int tcgettmixms_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TC_getTmixMs( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

int tcdndgettmixml_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TCDND_getTmixMl( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

int tcgettmixml_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TC_getTmixMl( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

int tcdndgetms2cpmixms_(double *scal,int *Nvars,double *cpmix) 
{
  int ans = 0 ;
  ans = TCDND_getMs2CpMixMs( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

int tcgetms2cpmixms_(double *scal,int *Nvars,double *cpmix) 
{
  int ans = 0 ;
  ans = TC_getMs2CpMixMs( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

int tcdndgetms2cvmixms_(double *scal,int *Nvars,double *cvmix) 
{
  int ans = 0 ;
  ans = TCDND_getMs2CvMixMs( scal, (*Nvars), cvmix) ;
  return ( ans ) ;
}

int tcgetms2cvmixms_(double *scal,int *Nvars,double *cvmix) 
{
  int ans = 0 ;
  ans = TC_getMs2CvMixMs( scal, (*Nvars), cvmix) ;
  return ( ans ) ;
}

int tcdndgetml2cpmixml_(double *scal,int *Nvars,double *cpmix) 
{
  int ans = 0 ;
  ans = TCDND_getMl2CpMixMl( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

int tcgetml2cpmixml_(double *scal,int *Nvars,double *cpmix) 
{
  int ans = 0 ;
  ans = TC_getMl2CpMixMl( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

int tcdndgetcpspecms_(double *t,int *Nspec,double *cpi) 
{
  int ans = 0 ;
  ans = TCDND_getCpSpecMs( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

int tcgetcpspecms_(double *t,int *Nspec,double *cpi) 
{
  int ans = 0 ;
  ans = TC_getCpSpecMs( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

int tcdndgetcpspecml_(double *t,int *Nspec,double *cpi) 
{
  int ans = 0 ;
  ans = TCDND_getCpSpecMl( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

int tcgetcpspecml_(double *t,int *Nspec,double *cpi) 
{
  int ans = 0 ;
  ans = TC_getCpSpecMl( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

int tcdndgetms2hmixms_(double *scal,int *Nvars,double *hmix) 
{
  int ans = 0 ;
  ans = TCDND_getMs2HmixMs( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

int tcgetms2hmixms_(double *scal,int *Nvars,double *hmix) 
{
  int ans = 0 ;
  ans = TC_getMs2HmixMs( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

int tcdndgetml2hmixml_(double *scal,int *Nvars,double *hmix) 
{
  int ans = 0 ;
  ans = TCDND_getMl2HmixMl( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

int tcgetml2hmixml_(double *scal,int *Nvars,double *hmix) 
{
  int ans = 0 ;
  ans = TC_getMl2HmixMl( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

int tcdndgethspecms_(double *t,int *Nspec,double *hi) 
{
  int ans = 0 ;
  ans = TCDND_getHspecMs( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

int tcgethspecms_(double *t,int *Nspec,double *hi) 
{
  int ans = 0 ;
  ans = TC_getHspecMs( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

int tcdndgethspecml_(double *t,int *Nspec,double *hi) 
{
  int ans = 0 ;
  ans = TCDND_getHspecMl( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

int tcgethspecml_(double *t,int *Nspec,double *hi) 
{
  int ans = 0 ;
  ans = TC_getHspecMl( *t, (*Nspec), hi) ;
  return ( ans ) ;
}


