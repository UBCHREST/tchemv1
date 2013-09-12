#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <assert.h>
#include "math.h"

#include "StiffInteg.h"
#include "TC_interface.h"
//#define DEBUG

int main()
{

  /* Integer parameters */
  int Nspec    ;
  int NiterMax ;
  int oFreq    ;

  /* Double parameters */
  double Tini, Temp_id, deltat, deltatMax, tEnd ;
  double CVrelt, CVsmall ;
  double deltaTemp,pressure,pfac ;

  /* Character array parameters */
  int sizecMax     = 100 ;
  char *mechfile   = new char[sizecMax] ; 
  char *thermofile = new char[sizecMax] ;

  /* Booleab parameters */
  bool withTab ;

  /* Species IC */
  std::vector<std::string> SpecName ;
  std::vector<double>      SpecMsFr ;

  /* Set default values */
  NiterMax  = 20000     ;
  oFreq     = 100       ;
  Tini      = 1000.0    ; /* [K]  */
  Temp_id   = 1500.0    ; /* [K]  */
  deltat    = 1.e-9     ; /* [s]  */
  deltatMax = 1.e-6     ; /* [s]  */
  tEnd      = 0.1e0     ; /* [s]  */
  deltaTemp = 10.0      ; /* [K]  */
  pressure  = 1.01325e5 ; /* [Pa] */
  pfac      = 1.0       ; /* [ ]  */

  CVrelt  = 1.e-12 ;
  CVsmall = 1.e-20 ;

  strcpy(mechfile  ,"chem.inp" ) ; 
  strcpy(thermofile,"therm.dat") ;

  withTab = false ; /* no tabulation by default */

  std::string keyName, sName ;
  double dtemp ;
  std::cin >> keyName ;

  while ( ( keyName != "end" ) && (keyName != "END") )
  {
    // integer parameters
    if ( keyName == "NiterMax"  ) std::cin >> NiterMax ;
    if ( keyName == "oFreq"     ) std::cin >> oFreq    ;

    // double parameters
    if ( keyName == "Tini"      ) std::cin >> Tini      ;
    if ( keyName == "Temp_id"   ) std::cin >> Temp_id   ;
    if ( keyName == "deltat"    ) std::cin >> deltat    ;
    if ( keyName == "deltatMax" ) std::cin >> deltatMax ;
    if ( keyName == "tEnd"      ) std::cin >> tEnd      ;
    if ( keyName == "CVrelt"    ) std::cin >> CVrelt    ;
    if ( keyName == "CVsmall"   ) std::cin >> CVsmall   ;
    if ( keyName == "deltaTemp" ) std::cin >> deltaTemp ;
    if ( keyName == "pfac"      ) std::cin >> pfac      ;

    // character string
    if ( keyName == "mech"     ) 
    {
      std::cin >> sName  ;   strcpy(mechfile,sName.c_str() ) ;
    }
    if ( keyName == "thermo"   ) 
    {
      std::cin >> sName  ;   strcpy(thermofile,sName.c_str() ) ;
    }
    if ( keyName == "withTab" ) std::cin >> withTab ;
    if ( keyName == "spec"   ) 
    {
      std::cin >> sName  ; SpecName.push_back(sName) ;
      std::cin >> dtemp  ; SpecMsFr.push_back(dtemp) ;
    }

    std::cin >> keyName ;

  }

  std::cout<<"------------------------------------------"<<std::endl<<std::flush ;
  std::cout<<"Run parameters : "<<std::endl<<std::flush ;
  std::cout<<std::endl<<std::flush ;
  std::cout<<"    NiterMax  = "<<NiterMax <<std::endl<<std::flush ;
  std::cout<<"    oFreq     = "<<oFreq    <<std::endl<<std::flush ;
  std::cout<<std::endl<<std::flush ;
  std::cout<<"    Tini      = "<<Tini     <<std::endl<<std::flush ;
  std::cout<<"    Temp_id   = "<<Temp_id  <<std::endl<<std::flush ;
  std::cout<<"    deltat    = "<<deltat   <<std::endl<<std::flush ;
  std::cout<<"    deltatMax = "<<deltatMax<<std::endl<<std::flush ;
  std::cout<<"    tEnd      = "<<tEnd     <<std::endl<<std::flush ;
  std::cout<<"    pfac      = "<<pfac     <<std::endl<<std::flush ;
  std::cout<<std::endl<<std::flush ;
  std::cout<<"    Kinetic model : "<<mechfile  <<std::endl<<std::flush ; 
  std::cout<<"    Thermo props  : "<<thermofile<<std::endl<<std::flush ;
  std::cout<<"    Create tables : "<<withTab   <<std::endl<<std::flush ;
  std::cout<<std::endl<<std::flush ;
  for (int i = 0 ; i<SpecName.size() ; i++)
    if ( fabs(SpecMsFr[i]) > 1.e-15 )
      std::cout<<"    X_"<<SpecName[i]<<" = "<<SpecMsFr[i]<<std::endl<<std::flush ;
  std::cout<<"------------------------------------------"<<std::endl<<std::flush ;

  /*
                           _               
                  ___  ___| |_ _   _ _ __  
                 / __|/ _ \ __| | | | '_ \ 
                 \__ \  __/ |_| |_| | |_) |
                 |___/\___|\__|\__,_| .__/ 
                                    |_|    

  */
  /* Initialize TC library */
  TC_initChem( mechfile, thermofile, (int) withTab, 0.2) ; 
  pressure *= pfac ;
  TC_setThermoPres(pressure) ;

  Nspec = TC_getNspec() ;
  std::cout<<"Nspec   = "<<Nspec <<std::endl<<std::flush ;

  double *scal = new double[Nspec+1] ;

  /* Instantiate stiff integrator */
  StiffInteg intgcvode( scal, Nspec, deltaTemp ) ;

/*
          _ __ _   _ _ __  
         | '__| | | | '_ \ 
         | |  | |_| | | | |
         |_|   \____|_| |_|
*/
  /* Set initial conditions */
  std::cout<<"Set initial conditions "<<std::endl<<std::flush ;
  scal[0] = Tini ;
  for (int i = 1 ; i<Nspec+1 ; i++) scal[i] = 0.0;

  for (int i = 0 ; i<SpecName.size() ; i++)
  {  
    int ispec = TC_getSpos( SpecName[i].c_str(), SpecName[i].size() ) ;
    if (ispec < 0 ) 
    {
      std::cout<<"Error : Could not find species "<<SpecName[i]<<std::endl ;
      std::terminate() ;
    }
    else
      std::cout<<"Index of species "<<std::setw(16)<<SpecName[i]<<" is "<<ispec<<std::endl ;
    scal[ispec+1] = SpecMsFr[i] ;
  }

  for (int i = 0 ; i < Nspec+1 ; i++) 
    if ( fabs(scal[i]) > 1.e-15 )
      std::cout<<i<<" "<<scal[i]<<std::endl;

  /* convert mole to mass fractions */
  double *msfr = new double[Nspec] ;
  TC_getMl2Ms( &(scal[1]), Nspec, msfr ) ;
  for (int i = 1 ; i < Nspec+1 ; i++) scal[i] = msfr[i-1];
  delete []msfr ;

  /* set tolerances for cvode */
  intgcvode.setRelT ( CVrelt  ) ; 
  intgcvode.setSmall( CVsmall ) ;

  /* send ignition threshold to cvode */
  intgcvode.setTempThresh( Temp_id ) ;

  /* Integrate */
  intgcvode.compute( tEnd, &deltat, deltatMax, NiterMax, oFreq ) ;

  /* Clean-up */
  delete [] mechfile   ;
  delete [] thermofile ;
  delete [] scal       ;

  /* Reset TC library */
  TC_reset() ;
  return ( 0 ) ;

}

