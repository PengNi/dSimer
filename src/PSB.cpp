#include <Rcpp.h>
using namespace Rcpp;
#include <string>

#ifdef __APPLE__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

using namespace std;

#include <vector>

// go_jaccardindex must be a symmetric matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix PSB_termsim_cpp(Rcpp::NumericMatrix& go_jaccardindex,
  Rcpp::List& IC){
  
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,float> informationcontent;
  #else
    typedef std::tr1::unordered_map<std::string,float> informationcontent;
  #endif
  
  informationcontent InfCon;
  {
    Rcpp::CharacterVector inames(IC.names());
    for(int i=0;i<inames.size();i++){
      InfCon.insert(std::make_pair(inames[i],Rcpp::as<float>( IC[i] )));
    }
  }  
  
  Rcpp::NumericMatrix sim(go_jaccardindex.nrow(),go_jaccardindex.ncol());
  Rcpp::List dimnames = go_jaccardindex.attr("dimnames");
  sim.attr("dimnames")=dimnames;
  
  Rcpp::CharacterVector rownames(dimnames[0]);
  std::vector<std::string> rowns=Rcpp::as< std::vector<std::string> >(rownames);
  
  
  for(std::size_t i=0;i<rowns.size();i++){
    for(std::size_t j=i;j<rowns.size();j++){
      float tmpval=go_jaccardindex(i,j)*((InfCon[rowns[i]]+InfCon[rowns[j]])/2);
      sim(i,j)=tmpval;
      sim(j,i)=tmpval;
    }
  }
  
  return(sim);
}












