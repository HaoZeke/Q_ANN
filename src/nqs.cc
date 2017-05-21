/*
############################ COPYRIGHT NOTICE ##################################

Code provided by G. Carleo and M. Troyer, written by G. Carleo, December 2016.

Permission is granted for anyone to copy, use, modify, or distribute the
accompanying programs and documents for any purpose, provided this copyright
notice is retained and prominently displayed, along with a complete citation of
the published version of the paper:
 ______________________________________________________________________________
| G. Carleo, and M. Troyer                                                     |
| Solving the quantum many-body problem with artificial neural-networks        |
|______________________________________________________________________________|

The programs and documents are distributed without any warranty, express or
implied.

These programs were written for research purposes only, and are meant to
demonstrate and reproduce the main results obtained in the paper.

All use of these programs is entirely at the user's own risk.

################################################################################
*/


#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cassert>

class Nqs{

  //Neural-network weights
  std::vector<std::vector<std::complex<double> > > W_;

  //Neural-network visible bias
  std::vector<std::complex<double> > a_;

  //Neural-network hidden bias
  std::vector<std::complex<double> > b_;

  //Number of hidden units
  int nh_;

  //Number of visible units
  int nv_;

  //look-up tables
  std::vector<std::complex<double> > Lt_;

  //Useful quantities for safe computation of ln(cosh(x))
  const double log2_;

public:

  Nqs(std::string filename):log2_(std::log(2.)){
    LoadParameters(filename);
  }

  //computes the logarithm of the wave-function
  inline std::complex<double> LogVal(const std::vector<int> & state)const{

    std::complex<double> rbm(0.,0.);

    for(int v=0;v<nv_;v++){
      rbm+=a_[v]*double(state[v]);
    }

    for(int h=0;h<nh_;h++){

      std::complex<double> thetah=b_[h];

      for(int v=0;v<nv_;v++){
        thetah+=double(state[v])*(W_[v][h]);
      }

      rbm+=Nqs::lncosh(thetah);
    }

    return rbm;
  }

  //computes the logarithm of Psi(state')/Psi(state)
  //where state' is a state with a certain number of flipped spins
  //the vector "flips" contains the sites to be flipped
  //look-up tables are used to speed-up the calculation
  inline std::complex<double> LogPoP(const std::vector<int> & state,const std::vector<int> & flips)const{

    if(flips.size()==0){
      return 0.;
    }

    std::complex<double> logpop(0.,0.);

    //Change due to the visible bias
    for(const auto & flip : flips){
      logpop-=a_[flip]*2.*double(state[flip]);
    }

    //Change due to the interaction weights
    for(int h=0;h<nh_;h++){
      const std::complex<double> thetah=Lt_[h];
      std::complex<double> thetahp=thetah;

      for(const auto & flip : flips){
        thetahp-=2.*double(state[flip])*(W_[flip][h]);
      }
      logpop+= ( Nqs::lncosh(thetahp)-Nqs::lncosh(thetah) );
    }

    return logpop;
  }

  inline std::complex<double> PoP(const std::vector<int> & state,const std::vector<int> & flips)const{
    return std::exp(LogPoP(state,flips));
  }

  //initialization of the look-up tables
  void InitLt(const std::vector<int> & state){
    Lt_.resize(nh_);

    for(int h=0;h<nh_;h++){
      Lt_[h]=b_[h];
      for(int v=0;v<nv_;v++){
        Lt_[h]+=double(state[v])*(W_[v][h]);
      }
    }

  }

  //updates the look-up tables after spin flips
  //the vector "flips" contains the indices of sites to be flipped
  void UpdateLt(const std::vector<int> & state,const std::vector<int> & flips){
    if(flips.size()==0){
      return;
    }

    for(int h=0;h<nh_;h++){
      for(const auto & flip : flips){
        Lt_[h]-=2.*double(state[flip])*W_[flip][h];
      }
    }
  }

  //loads the parameters of the wave-function from a given file
  void LoadParameters(std::string filename){

    std::ifstream fin(filename.c_str());

    if(!fin.good()){
      std::cerr<<"# Error : Cannot load from file "<<filename<<" : file not found."<<std::endl;
      std::abort();
    }

    fin>>nv_;
    fin>>nh_;

    if(!fin.good() || nv_<0 || nh_<0){
      std::cerr<<"# Trying to load from an invalid file.";
      std::cerr<<std::endl;
      std::abort();
    }

    a_.resize(nv_);
    b_.resize(nh_);
    W_.resize(nv_,std::vector<std::complex<double> > (nh_));

    for(int i=0;i<nv_;i++){
      fin>>a_[i];
    }
    for(int j=0;j<nh_;j++){
      fin>>b_[j];
    }
    for(int i=0;i<nv_;i++){
      for(int j=0;j<nh_;j++){
        fin>>W_[i][j];
      }
    }

    if(!fin.good()){
      std::cerr<<"# Trying to load from an invalid file.";
      std::cerr<<std::endl;
      std::abort();
    }

    std::cout<<"# NQS loaded from file "<<filename<<std::endl;
    std::cout<<"# N_visible = "<<nv_<<"  N_hidden = "<<nh_<<std::endl;
  }

  //ln(cos(x)) for real argument
  //for large values of x we use the asymptotic expansion
  inline double lncosh(double x)const{
    const double xp=std::abs(x);
    if(xp<=12.){
      return std::log(std::cosh(xp));
    }
    else{
      return xp-log2_;
    }
  }

  //ln(cos(x)) for complex argument
  //the modulus is computed by means of the previously defined function
  //for real argument
  inline std::complex<double> lncosh(std::complex<double> x)const{
    const double xr=x.real();
    const double xi=x.imag();

    std::complex<double> res=Nqs::lncosh(xr);
    res +=std::log( std::complex<double>(std::cos(xi),std::tanh(xr)*std::sin(xi)) );

    return res;
  }

  //total number of spins
  //equal to the number of visible units
  inline int Nspins()const{
    return nv_;
  }

};
