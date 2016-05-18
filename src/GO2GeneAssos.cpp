#include <Rcpp.h>
using namespace Rcpp;
#include <string>
using namespace std;

#ifdef __APPLE__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include <vector>
#include <set>

// [[Rcpp::export]]
Rcpp::List go2g_full(Rcpp::List& go2g,
                     Rcpp::List& go_offsp) {
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,std::vector<std::string> > cgo2g;
  #else
    typedef std::tr1::unordered_map<std::string,std::vector<std::string> > cgo2g;
  #endif
  cgo2g go2g_c;
  {
    Rcpp::CharacterVector goNames(go2g.names());
    for(int i=0;i<goNames.size();i++){
      const std::vector<std::string> genes = Rcpp::as<std::vector<std::string> >( go2g[i] );
      go2g_c.insert( std::make_pair( goNames[i], genes ) );
    }
  }
  
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,std::vector<std::string> > cgooffsp;
  #else
    typedef std::tr1::unordered_map<std::string,std::vector<std::string> > cgooffsp;
  #endif
  cgooffsp gooffsp_c;
  {
    Rcpp::CharacterVector goNames(go_offsp.names());
    for(int i=0;i<goNames.size();i++){
      const std::vector<std::string> offs = Rcpp::as<std::vector<std::string> >( go_offsp[i] );
      gooffsp_c.insert( std::make_pair( goNames[i], offs ) );
    }
  }
  
  Rcpp::List go2gfull;
  
  Rcpp::CharacterVector gonames(go2g.names());
  int golen=gonames.size();
  for(int i=0;i<golen;i++){
    std::set<std::string> genes;
    std::vector<std::string> offs=gooffsp_c[ Rcpp::as<std::string>(gonames[i]) ];
    int count=offs.size();
    for(int j=0;j<count;j++){
      std::vector<std::string> gtemp=go2g_c[ offs[j] ];
      genes.insert(gtemp.begin(),gtemp.end());
    }
    std::vector<std::string> gtemp=go2g_c[Rcpp::as<std::string>(gonames[i])];
    genes.insert(gtemp.begin(),gtemp.end());
    
    go2gfull[Rcpp::as<std::string>(gonames[i])]=genes;
  }
  
  return(go2gfull);
}

