#include <Rcpp.h>
#include <string>
using namespace std;

#ifdef __APPLE__
  #include <unordered_map>
#else
  #include <tr1/unordered_map>
#endif

using namespace Rcpp;
#include <algorithm>
#include <vector>
#include <set>
#include <iterator>


// [[Rcpp::export]]
Rcpp::NumericMatrix BOG_simmat_cpp(
  Rcpp::List& d2g,
  int totalGeneNum,
  Rcpp::List& GeneNumList){
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,std::vector<std::string> > Disease2Gene;
  #else
    typedef std::tr1::unordered_map<std::string,std::vector<std::string> > Disease2Gene;
  #endif
  
  Rcpp::CharacterVector dnames(d2g.names());
  Disease2Gene dgMap;
  {
    Rcpp::CharacterVector dNames(d2g.names());
    for(int i=0;i<dNames.size();i++){
      const std::vector<std::string> genes = Rcpp::as<std::vector<std::string> >( d2g[i] );
      dgMap.insert( std::make_pair( dNames[i], genes ) );
    }
  }
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,int> D2GNum;
  #else
    typedef std::tr1::unordered_map<std::string,int> D2GNum;
  #endif
  D2GNum d2gnum;
  {
    Rcpp::CharacterVector dNames(GeneNumList.names());
    for(int i=0;i<dNames.size();i++){
      d2gnum.insert(std::make_pair(dNames[i],Rcpp::as<int>( GeneNumList[i] )));
    }
  }
  
  int totalnum=totalGeneNum;
  
  //std::cout<<totalnum;
  
  int dlen=dnames.size();
  
  //std::cout<<dlen;
  Rcpp::NumericMatrix sim(dlen,dlen);
  sim.attr("dimnames")=Rcpp::Rcpp_list2(dnames,dnames);
  
  for(int i=0;i<dlen;i++){
    for(int j=i;j<dlen;j++){
      
      float denominator= (float) (d2gnum[ Rcpp::as<std::string>(dnames[i]) ]
                                  *d2gnum[ Rcpp::as<std::string>(dnames[j]) ])/
                                                  (totalnum*totalnum);
      
      std::vector<std::string> ga=dgMap[ Rcpp::as<std::string>(dnames[i]) ];
      std::vector<std::string> gb=dgMap[ Rcpp::as<std::string>(dnames[j]) ];
      
      
      std::sort(ga.begin(), ga.end());
      std::sort(gb.begin(), gb.end());
      std::vector<std::string> s;

      set_intersection(ga.begin(), ga.end(), gb.begin(), gb.end(), 
                        std::back_inserter(s));
      
      std::set<std::string> g(ga.begin(),ga.end());
      g.insert(gb.begin(),gb.end());
      
      float jc = (float)s.size()/g.size();
      
      if(denominator==0){
        sim(i,j)=0;
        sim(j,i)=0;
      }else{
        jc /= denominator;
        sim(i,j)=jc;
        sim(j,i)=jc;
      }
      
    }
  }
  return(sim);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix BOG_normat_cpp(Rcpp::CharacterVector D1,
  Rcpp::CharacterVector D2,
  Rcpp::NumericMatrix& simmat,
  Rcpp::List& matnameloc,
  Rcpp::List& maxsim,
  Rcpp::List& IC){
  
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,int> matnamelocation;
  #else
    typedef std::tr1::unordered_map<std::string,int> matnamelocation;
  #endif
    
    matnamelocation MatNameLoc;
    {
      Rcpp::CharacterVector mnames(matnameloc.names());
      for(int i=0;i<mnames.size();i++){
        MatNameLoc.insert(std::make_pair(mnames[i],Rcpp::as<int>( matnameloc[i] )));
      }
    }
    
    #ifdef __APPLE__
      typedef std::unordered_map<std::string,float> maxsimvalue;
    #else
      typedef std::tr1::unordered_map<std::string,float> maxsimvalue;
    #endif
    
    maxsimvalue MaxSimVal;
    {
      Rcpp::CharacterVector mnames(maxsim.names());
      for(int i=0;i<mnames.size();i++){
        MaxSimVal.insert(std::make_pair(mnames[i],Rcpp::as<float>( maxsim[i] )));
      }
    }
    
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
    
    Rcpp::NumericMatrix sim(D1.size(),D2.size());
    sim.attr("dimnames")=Rcpp::Rcpp_list2(D1,D2);
    
    
    std::vector<std::string> ds1=Rcpp::as< std::vector<std::string> > (D1);
    std::vector<std::string> ds2=Rcpp::as< std::vector<std::string> > (D2);
    
    if(ds1==ds2){
      //std::cout<<"ds1==ds2"<<std::endl;
      for(int i=0;i<D1.size();i++){
        for(int j=i;j<D1.size();j++){
          std::string d1=Rcpp::as<std::string>(D1[i]);
          std::string d2=Rcpp::as<std::string>(D1[j]);
          
          float tmpval = (float) simmat(MatNameLoc[d1],MatNameLoc[d2])/
                                 ((MaxSimVal[d1]+MaxSimVal[d2])/2);
          tmpval *=((InfCon[d1]+InfCon[d2])/2);
          sim(i,j)=tmpval;
          sim(j,i)=tmpval;
        }
      }
    }else{
       for(int i=0;i<D1.size();i++){
        for(int j=0;j<D2.size();j++){
          std::string d1=Rcpp::as<std::string>(D1[i]);
          std::string d2=Rcpp::as<std::string>(D2[j]);
          
          float tmpval = (float) simmat(MatNameLoc[d1],MatNameLoc[d2])/((MaxSimVal[d1]+MaxSimVal[d2])/2);
          tmpval *=((InfCon[d1]+InfCon[d2])/2);
          sim(i,j)=tmpval;
        }
      }
    }
    return(sim);
  
  }










