#include <Rcpp.h>
#include <string>
using namespace Rcpp;
using namespace std;

#ifdef __APPLE__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include <numeric>
#include <math.h> 
#include <algorithm>
#include <vector>
#include <set>
#include <iterator>

const int sigdim=73;

#ifdef __APPLE__
float toposim_onepair(std::vector<std::string>& gs1,
                      std::vector<std::string>& gs2,
                      std::unordered_map< std::string,std::vector<float> >& graphlet_sig,
                      std::vector<float>& weight_v,
                      float& weight_sum){
  
  int g1size=gs1.size(),g2size=gs2.size();
  
  float sigsim[g1size][g2size];
  std::unordered_map< std::string,int > gs1loc;
  {
    for(int i=0;i<g1size;i++){
      gs1loc.insert(std::make_pair(gs1[i],i));
    }
  }
  std::unordered_map< std::string,int > gs2loc;
  {
    for(int i=0;i<g2size;i++){
      gs2loc.insert(std::make_pair(gs2[i],i));
    }
  }
  
  for(int i=0;i<g1size;i++){
    std::unordered_map< std::string,std::vector<float> >::const_iterator it1=
      graphlet_sig.find(gs1[i]);
    if(it1==graphlet_sig.end()){
      for(int j=0;j<g2size;j++){
        if(gs1[i]==gs2[j]){
          sigsim[i][j]=1;
        }else{
          sigsim[i][j]=0;
        }
      }
    }else{
      for(int j=0;j<g2size;j++){
        if(gs1[i]==gs2[j]){
          sigsim[i][j]=1;
        }else{
          std::unordered_map< std::string,std::vector<float> >::const_iterator it2=
            graphlet_sig.find(gs2[j]);
          if(it2==graphlet_sig.end()){
            sigsim[i][j]=0;
          }else{
            std::vector<float> uv1=it1->second;
            std::vector<float> uv2=it2->second;
            float res=0;
            for(int ii=0;ii<sigdim;ii++){
              res += weight_v[ii]*abs(log(uv1[ii]+1)-log(uv2[ii]+1))/log(max(uv1[ii],uv2[ii])+2);
            }
            sigsim[i][j] =1-res/weight_sum;
          }
        }
      }
    }
  }
  
  /*for(int i=0;i<g1size;i++){
  for(int j=0;j<g2size;j++){
  std::cout<<sigsim[i][j]<<" ";
  }
  std::cout<<std::endl;
}*/
  
  
  float max1=0;
  
  for(int i=0;i<g1size;i++){
    for(int j=0;j<g2size;j++){
      if(max1<sigsim[i][j]){
        max1=sigsim[i][j];
      }
    }
  }
  
  float max2=0;
  std::sort(gs1.begin(), gs1.end());
  std::sort(gs2.begin(), gs2.end());
  std::vector<std::string> s;
  
  set_intersection(gs1.begin(), gs1.end(), gs2.begin(), gs2.end(), 
                   std::back_inserter(s));
  int ssizeM1=s.size()-1;
  if(ssizeM1>0){
    for(int i=0;i<ssizeM1;i++){
      for(std::size_t j=i+1;j<s.size();j++){
        float tmp=sigsim[gs1loc[s[i]]][gs2loc[s[j]]];
        if(max2<tmp){
          max2=tmp;
        }
      }
    }
  }
  
  //std::cout<<max1<<" "<<max2<<std::endl;
  return((max1+max2)/2);
  }
#else
float toposim_onepair(std::vector<std::string>& gs1,
                      std::vector<std::string>& gs2,
                      std::tr1::unordered_map< std::string,std::vector<float> >& graphlet_sig,
                      std::vector<float>& weight_v,
                      float& weight_sum){
  
  int g1size=gs1.size(),g2size=gs2.size();
  
  float sigsim[g1size][g2size];
  std::tr1::unordered_map< std::string,int > gs1loc;
  {
    for(int i=0;i<g1size;i++){
      gs1loc.insert(std::make_pair(gs1[i],i));
    }
  }
  std::tr1::unordered_map< std::string,int > gs2loc;
  {
    for(int i=0;i<g2size;i++){
      gs2loc.insert(std::make_pair(gs2[i],i));
    }
  }
  
  for(int i=0;i<g1size;i++){
    std::tr1::unordered_map< std::string,std::vector<float> >::const_iterator it1=
      graphlet_sig.find(gs1[i]);
    if(it1==graphlet_sig.end()){
      for(int j=0;j<g2size;j++){
        if(gs1[i]==gs2[j]){
          sigsim[i][j]=1;
        }else{
          sigsim[i][j]=0;
        }
      }
    }else{
      for(int j=0;j<g2size;j++){
        if(gs1[i]==gs2[j]){
          sigsim[i][j]=1;
        }else{
          std::tr1::unordered_map< std::string,std::vector<float> >::const_iterator it2=
            graphlet_sig.find(gs2[j]);
          if(it2==graphlet_sig.end()){
            sigsim[i][j]=0;
          }else{
            std::vector<float> uv1=it1->second;
            std::vector<float> uv2=it2->second;
            float res=0;
            for(int ii=0;ii<sigdim;ii++){
              res += weight_v[ii]*abs(log(uv1[ii]+1)-log(uv2[ii]+1))/log(max(uv1[ii],uv2[ii])+2);
            }
            sigsim[i][j] =1-res/weight_sum;
          }
        }
      }
    }
  }
  
  /*for(int i=0;i<g1size;i++){
  for(int j=0;j<g2size;j++){
  std::cout<<sigsim[i][j]<<" ";
  }
  std::cout<<std::endl;
}*/
  
  
  float max1=0;
  
  for(int i=0;i<g1size;i++){
    for(int j=0;j<g2size;j++){
      if(max1<sigsim[i][j]){
        max1=sigsim[i][j];
      }
    }
  }
  
  float max2=0;
  std::sort(gs1.begin(), gs1.end());
  std::sort(gs2.begin(), gs2.end());
  std::vector<std::string> s;
  
  set_intersection(gs1.begin(), gs1.end(), gs2.begin(), gs2.end(), 
                   std::back_inserter(s));
  int ssizeM1=s.size()-1;
  if(ssizeM1>0){
    for(int i=0;i<ssizeM1;i++){
      for(std::size_t j=i+1;j<s.size();j++){
        float tmp=sigsim[gs1loc[s[i]]][gs2loc[s[j]]];
        if(max2<tmp){
          max2=tmp;
        }
      }
    }
  }
  
  //std::cout<<max1<<" "<<max2<<std::endl;
  return((max1+max2)/2);
  }
#endif

// [[Rcpp::export]]
Rcpp::NumericMatrix Sun_topology_cpp(Rcpp::CharacterVector& D1,
                                     Rcpp::CharacterVector& D2,
                                     Rcpp::List& d2g,
                                     Rcpp::NumericMatrix& graphlet_signature,
                                     Rcpp::NumericVector& weight){
  
  std::vector<float> weight_v=Rcpp::as< std::vector<float> >(weight);
  
  float weight_sum=std::accumulate(weight_v.begin(), weight_v.end(),(double) 0);
  
  
  #ifdef __APPLE__
    typedef std::unordered_map<std::string,std::vector<std::string> > Disease2Gene;
  #else
    typedef std::tr1::unordered_map<std::string,std::vector<std::string> > Disease2Gene;
  #endif
  
  Disease2Gene dgMap;
  {
    Rcpp::CharacterVector dNames(d2g.names());
    for(int i=0;i<dNames.size();i++){
      const std::vector<std::string> genes = Rcpp::as<std::vector<std::string> >( d2g[i] );
      dgMap.insert( std::make_pair( dNames[i], genes ) );
    }
  }
  
  #ifdef __APPLE__
    typedef std::unordered_map< std::string,std::vector<float> > graphlet_sig;
  #else
    typedef std::tr1::unordered_map< std::string,std::vector<float> > graphlet_sig;
  #endif
  graphlet_sig gletsig;
  {
    Rcpp::List dimnames = graphlet_signature.attr("dimnames");
    Rcpp::CharacterVector gNames(dimnames[0]);
    for(int i=0;i<gNames.size();i++){
      Rcpp::NumericMatrix::Row r=graphlet_signature.row(i);
      Rcpp::NumericVector rvec=r;
      const std::vector<float> u=Rcpp::as< std::vector<float> >(rvec);
      gletsig.insert( std::make_pair(gNames[i], u ));
    }
  }
  
  
  Rcpp::NumericMatrix sim(D1.size(),D2.size());
  sim.attr("dimnames")=Rcpp::Rcpp_list2(D1,D2);
  
  std::vector<std::string> ds1=Rcpp::as< std::vector<std::string> > (D1);
  std::vector<std::string> ds2=Rcpp::as< std::vector<std::string> > (D2);
  if(ds1==ds2){
    for(int i=0;i<D1.size();i++){
      for(int j=i;j<D1.size();j++){
        if(D1[i]==D1[j]){
          sim(i,j)=1;
          sim(j,i)=1;
        }else{
          std::vector<std::string> gtmp1=dgMap[ Rcpp::as<std::string>(D1[i]) ];
          std::vector<std::string> gtmp2=dgMap[ Rcpp::as<std::string>(D1[j]) ];
          
          float resval=toposim_onepair(gtmp1,gtmp2,gletsig,weight_v,weight_sum);
          sim(i,j)=resval;
          sim(j,i)=resval;
        }
      }
    }
  }else{
    for(int i=0;i<D1.size();i++){
      for(int j=0;j<D2.size();j++){
        if(D1[i]==D1[j]){
          sim(i,j)=1;
          sim(j,i)=1;
        }else{
          std::vector<std::string> gtmp1=dgMap[ Rcpp::as<std::string>(D1[i]) ];
          std::vector<std::string> gtmp2=dgMap[ Rcpp::as<std::string>(D2[j]) ];
          sim(i,j)=toposim_onepair(gtmp1,gtmp2,gletsig,weight_v,weight_sum);
        }
      }
    }
  }
  
  return(sim);
}











