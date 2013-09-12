#include "stdlib.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace boost::algorithm;
using namespace std;
using std::ifstream;

#include <cstring>

const int MAX_CHARS_PER_LINE = 1024;

void getTextFile(char *filein, vector<string> &ftext) ;
void getSection(vector<string> mech, string sname, int *sst, int *sen) ;
void getReacLines(vector<string> mech, int rst, int ren, vector<int> &reacIdx) ;
void checkSpecInList(string sline, vector<string> specDelList, vector<string> &validSpec) ;
bool checkValidReac(vector<string> mechfile, int ireac, vector<string> specDelList) ;
string delSpecFromReac(string mline, vector<string> specDelList) ;
bool isSpecInVecStr(vector<string> rcomp, vector<string> specDelList) ;
size_t firstCap(string strin);

void TC_reduce(char *mechIn, char *slist, char *mechOut) {

  vector<string> mechfile, specDelList;
  /* read mechanism from file */
  getTextFile(mechIn, mechfile) ;
  cout << "Done reading mechanism" << endl;
  /* read species list */
  getTextFile(slist,  specDelList) ;
  cout << "Done reading list of species" << endl;
  
  vector<int> reacIdx;
  int ElemStart, ElemEnd, SpecStart, SpecEnd, ReacStart, ReacEnd;
  getSection(mechfile,string("ELEM"),&ElemStart,&ElemEnd);
  getSection(mechfile,string("SPEC"),&SpecStart,&SpecEnd);
  getSection(mechfile,string("REAC"),&ReacStart,&ReacEnd);
  getReacLines(mechfile,ReacStart,ReacEnd,reacIdx);

#ifdef VERBOSE
  cout<<"Elements  :"<<ElemStart<<"->"<<ElemEnd<<endl<<flush;
  cout<<"Species   :"<<SpecStart<<"->"<<SpecEnd<<endl<<flush;
  cout<<"Reactions :"<<ReacStart<<"->"<<ReacEnd<<endl<<flush;
  cout<<"Done itemizing reactions"<<endl<<flush;
#endif


  ofstream fileout;
  fileout.open(mechOut);

  /* write elements */
  fileout << "ELEMENTS" <<endl;
  for ( int i = ElemStart; i<=ElemEnd; i++ )
    fileout << mechfile[i] << endl;
  fileout << "END" <<endl;
  
  /* write species */
  fileout << "SPECIES" <<endl;
  for ( int i = SpecStart; i<=SpecEnd; i++ ) {
    vector<string> validSpec;
    checkSpecInList(mechfile[i],specDelList,validSpec);
    for ( int j = 0; j < (int) validSpec.size(); j++ )
      fileout << validSpec[j] <<" " ;
    fileout << endl;
  }
  fileout << "END" <<endl;
  
  /* write reactions */
  fileout << "REACTIONS" <<endl;
  for ( int i = 0; i< (int) reacIdx.size()-1; i++ ) {

    vector<string> validSpec;
    if ( checkValidReac(mechfile,reacIdx[i],specDelList) ) {

      fileout << mechfile[reacIdx[i]] << endl;
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
	string cleanLine = delSpecFromReac(mechfile[j],specDelList);
	if ( cleanLine.size() > 0 )
	  fileout << cleanLine << endl ;
      }

    }

  } 
  fileout << "END" <<endl;
  
  /* done */
  fileout.close() ;

  return ;
  
}

void getTextFile(char *filein, vector<string> &ftext) {

  ifstream fin ;
  fin.open(filein); 
  if (!fin.good()) {
    cout << "getTextFile() Could not open "<<filein<<endl<<flush ;
    exit(0);
  }

  while (!fin.eof())
  {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    string bufstr(buf) ;
    
    /* Trim extra spaces */
    trim(bufstr);
    size_t found=bufstr.find(string("!"));
    if ( found != string::npos ) {
      bufstr.erase (found, bufstr.length());
    }
    if ( bufstr.size() == 0 ) continue ;
    transform(bufstr.begin(), bufstr.end(),bufstr.begin(), ::toupper);
    ftext.push_back(bufstr);
  }
  fin.close();
  return ;

}

void getSection(vector<string> mech, string sname, int *sst, int *sen) {
  
  bool foundSection = false;
  int i = 0 ;
  *sst = *sen = 0;
  while ( ( i < (int) mech.size() ) && ( *sen == 0 ) ){
     
    size_t sFound=mech[i].find(sname) ;
    if ( sFound != string::npos) {
      foundSection = true ;
      *sst = i+1;
    }
    else if ( foundSection ) {
      size_t eFound=mech[i].find(string("END")) ;
      if ( eFound != string::npos) {
	foundSection = false ;
	*sen = i-1;
      }
    }
    i++ ;
  }

  return ;
  
}
 
void getReacLines(vector<string> mech, int rst, int ren, vector<int> &reacIdx) {

  if ( ( rst == 0 ) && ( ren == 0) ) return ;

  for ( int i = rst; i <= ren; i++ ) {
    size_t rFound=mech[i].find(string("=")) ;
    if ( rFound != string::npos) reacIdx.push_back(i) ;
  }
  reacIdx.push_back(ren+1);

  return ;

}

void checkSpecInList(string sline, vector<string> specDelList, vector<string> &validSpec) {
  
  split( validSpec, sline, is_any_of(" ,/<=>"), token_compress_on );

  for ( int j = 0; j < (int) specDelList.size(); j++ )
    for ( int i = 0; i < (int) validSpec.size(); i++ )
      if ( validSpec[i] == specDelList[j] ) {
	validSpec.erase(validSpec.begin()+i);
        i--;
      }
    
  return ;

}

bool checkValidReac(vector<string> mechfile, int ireac, vector<string> specDelList) {

  vector<string> rcomp;

  split( rcomp, mechfile[ireac], is_any_of(" <=>"), token_compress_on );
  string react=rcomp[0], prod=rcomp[1];
  rcomp.clear();


  string stmp;
  size_t findTrdB;

  /* Check third-body in reactants */
  findTrdB=react.find(string("(+"));
  if ( findTrdB != string::npos ) {
    stmp=react.substr(findTrdB) ;
    react.erase(react.begin()+findTrdB,react.end());
    stmp.erase(stmp.begin(),stmp.begin()+2);
    stmp.erase(stmp.end()-1,stmp.end());
    rcomp.push_back(stmp);
    if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
    stmp.clear();
    rcomp.clear();
  } 

  /* Check third-body in products */
  findTrdB=prod.find(string("(+"));
  if ( findTrdB != string::npos ) {
    stmp=prod.substr(findTrdB) ;
    prod.erase(prod.begin()+findTrdB,prod.end());
    stmp.erase(stmp.begin(),stmp.begin()+2);
    stmp.erase(stmp.end()-1,stmp.end());
    rcomp.push_back(stmp);
    if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
    stmp.clear();
    rcomp.clear();
  } 

  split( rcomp, react, is_any_of("+"), token_compress_on );
  if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
  rcomp.clear();

  split( rcomp, prod, is_any_of("+"), token_compress_on );
  if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
  rcomp.clear();

  return (true) ;

}    

bool isSpecInVecStr(vector<string> rcomp, vector<string> specDelList) {

  for (int i = 0 ; i < (int) rcomp.size(); i++) {
    if (firstCap(rcomp[i]) > rcomp[i].size() ) {
      cout<<"error parsing reaction file for "<<endl<<rcomp[i]<<endl<<flush ;
    }
    string spec=rcomp[i].substr(firstCap(rcomp[i]));
    for  (int j = 0; j < (int) specDelList.size(); j++) {
      if (specDelList[j] == spec) return ( true ) ;
    }
  }

  return ( false );

}

string delSpecFromReac(string mline, vector<string> specDelList) {

  vector<string> rcomp;

  split( rcomp, mline, is_any_of(" /"), token_compress_on );

  bool foundmatch=false;
  for (int i = 0 ; i < (int) rcomp.size()-1; i++) {
    for  (int j = 0; j < (int) specDelList.size(); j++ ) 
      if (specDelList[j] == rcomp[i]) foundmatch=true;
    //cout<<"\""<<rcomp[i]<<"\" "<<foundmatch<<endl ;
  }


  if ( foundmatch ) {
    assert((rcomp.size()-1)%2==0);
    string rstr;
    for (int i = 0 ; i < (int) rcomp.size()/2; i++) {
      foundmatch=false;
      for  (int j = 0; j < (int) specDelList.size(); j++ ) 
	if (specDelList[j] == rcomp[2*i]) foundmatch=true;
      if ( !foundmatch ) {
        rstr.append(rcomp[2*i].c_str()) ;
        rstr.append(" /") ;
        rstr.append(rcomp[2*i+1].c_str()) ;
        rstr.append("/ ") ;
      }
    }
    //cout << rstr << endl ;
    return (rstr);
  }
  else
    return ( mline );

}    

size_t firstCap(string strin){

  size_t i=0;
  while (strin[i])
  {
    if (isupper(strin[i])) return (i);
    i++;
  }
  return (strin.size()+1);

}

