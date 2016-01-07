
#include <Rcpp.h>
#include <string>
#include <tr1/unordered_map>
using namespace Rcpp;
using namespace std;



// [[Rcpp::export]]
Rcpp::NumericMatrix funsim_cpp(
  Rcpp::CharacterVector d1,
  Rcpp::CharacterVector d2,
  Rcpp::List D2GList,
  Rcpp::List LLSn){
  
  typedef std::string gene_id;
  typedef std::string disease_id;
  typedef std::vector<std::string> gene_set;
  typedef std::tr1::unordered_map<disease_id, gene_set> Disease2Gene;
  
  Disease2Gene dgMap;
  {
    Rcpp::CharacterVector dNames(D2GList.names());
    for(int i=0;i<dNames.size();i++){
      const std::vector<std::string> genes = Rcpp::as<std::vector<std::string> >( D2GList[i] );
      
      dgMap.insert( std::make_pair( dNames[i], genes ) );
    }
  }
  
  typedef std::string gene_pair;
  typedef std::tr1::unordered_map<gene_pair,double> gpLLSn;
  gpLLSn gpSim;
  {
    Rcpp::CharacterVector gps(LLSn.names());
    for(int i=0;i<gps.size();i++){
      gpSim.insert(std::make_pair(gps[i],Rcpp::as<double>(LLSn[i])));
    }
  }
  
  
  int d1len=d1.size();
  int d2len=d2.size();
  Rcpp::NumericMatrix sim(d1len,d2len);
  sim.attr("dimnames")=Rcpp::Rcpp_list2(d1,d2);
  
  int count=0;
  
  for(int i=0;i<d1len;i++){
    const std::string dname1=(std::string) d1[i];
    const gene_set gs1=dgMap.find(dname1)->second;
    
    if(gs1.size()==0){
      for(int j=0;j<d2len;j++){
        if(dname1.compare((std::string)d2[j])==0){
          sim(i,j)=1;
          //cout<<sim(i,j)<<"what"<<endl;
        }else{
          sim(i,j)=0;
        }
      }
    }else{
      
      for(int j=0;j<d2len;j++){
        const std::string dname2=(std::string)d2[j];
        const gene_set gs2=dgMap.find(dname2)->second;
        
        if(dname2.compare(dname1)==0){
          sim(i,j)=1;
          //cout<<sim(i,j)<<dname1<<dname2<<endl;
        }else if(gs2.size()==0){
          sim(i,j)=0;
        }else{
          double score=0;
          
          for(std::size_t i1=0;i1<gs1.size();i1++){
            double max=0;
            for(std::size_t j1=0;j1<gs2.size();j1++){
              if((gs1[i1]).compare(gs2[j1])==0){
                max=1;
                break;
              }else{
                
                count++;
                
                std::string key1=(gs1[i1])+"|"+(gs2[j1]);
                std::string key2=(gs2[j1])+"|"+(gs1[i1]);
                double tempval=gpSim[key1];
                if(tempval){
                  if(tempval>max){
                    max=tempval;
                  }
                }else{
                  tempval=gpSim[key2];
                  if(tempval){
                    if(tempval>max){
                      max=tempval;
                    }
                  }
                }
              }
            }
            score += max;
          }
          
          for(std::size_t j1=0;j1<gs2.size();j1++){
            double max=0;
            for(std::size_t i1=0;i1<gs1.size();i1++){
              //use * or not
              if((gs1[i1]).compare(gs2[j1])==0){
                max=1;
                break;
              }else{
                
                count++;
                
                std::string key1=(gs1[i1])+"|"+(gs2[j1]);
                std::string key2=(gs2[j1])+"|"+(gs1[i1]);
                double tempval=gpSim[key1];
                if(tempval){
                  if(tempval>max){
                    max=tempval;
                  }
                }else{
                  tempval=gpSim[key2];
                  if(tempval){
                    if(tempval>max){
                      max=tempval;
                    }
                  }
                }
              }
            }
            score += max;
          }
          sim(i,j)=score/(gs1.size()+gs2.size());
        }
      }
      //std::cout<<i<<std::endl;
    }
    
  }
  
  return(sim);
}












