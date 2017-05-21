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
#include <vector>
#include <complex>

//Anti-ferromagnetic Heisenberg model in 1d
class Heisenberg2d{

  //number of spins
  const int nspins_;

  //size of the square lattice
  const int l_;

  //option to use periodic boundary conditions
  const bool pbc_;

  //coupling constant
  const double jz_;

  //Nearest neighbors on the square lattice
  std::vector<std::vector<int> > nn_;

  std::vector<std::vector<int> > bonds_;
public:

  Heisenberg2d(int nspins,double jz,bool pbc=true):nspins_(nspins),pbc_(pbc),jz_(jz),l_(std::sqrt(nspins_)){
    InitLattice();
    std::cout<<"# Using the 2d Heisenberg model with J_z = "<<jz_<<std::endl;
  }

  void InitLattice(){
    if(l_*l_!=nspins_){
      std::cerr<<"# Error , the number of spins is not compabitle with a square lattice "<<std::endl;
      std::abort();
    }

    nn_.resize(nspins_,std::vector<int>(4));

    //Defining the nearest-neighbors
    for(int i=0;i<nspins_;i++){
      nn_[i][0]=(pbc_)?PbcH(i-1,i):H(i-1,i);
      nn_[i][1]=(pbc_)?PbcH(i+1,i):H(i+1,i);
      nn_[i][2]=(pbc_)?PbcV(i-l_):V(i-l_);
      nn_[i][3]=(pbc_)?PbcV(i+l_):V(i+l_);
    }

    for(int i=0;i<nspins_;i++){
      for(int k=0;k<4;k++){
        int j=nn_[i][k];
        if(i<j){
          bonds_.push_back(std::vector<int>{i,j});
        }
      }
    }

  }


  //Finds the non-zero matrix elements of the hamiltonian
  //on the given state
  //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0
  //state' is encoded as the sequence of spin flips to be performed on state
  void FindConn(const std::vector<int> & state,std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){

    mel.resize(1);
    flipsh.resize(1);

    //computing interaction part Sz*Sz
    mel[0]=0.;

    for(int i=0;i<bonds_.size();i++){
      mel[0]+=double(state[bonds_[i][0]]*state[bonds_[i][1]]);
    }

    mel[0]*=jz_;

    //Looks for possible spin flips
    for(int i=0;i<bonds_.size();i++){
      const int si=bonds_[i][0];
      const int sj=bonds_[i][1];

      if(state[si]!=state[sj]){
        mel.push_back(-2);
        flipsh.push_back(std::vector<int>({si,sj}));
      }
    }

  }

  int MinFlips()const{
    return 2;
  }


  //Small functions to set up the lattice
  //Horizontal Pbc
  inline int PbcH(int nn,int s)const{
    if(s%l_==0 && nn==(s-1))
      return (s+l_-1);
    else if((s+1)%l_==0 && nn==(s+1))
      return (s-l_+1);
    else
      return nn;
  }

  //Vertical Pbc
  inline int PbcV(int nn)const{
    if(nn>=nspins_)
      return (nn-nspins_);
    else if(nn<0)
      return (nspins_+nn);
    else
      return nn;
  }

  //Horizontal without pbc
   inline int H(int nn,int s)const{
    if(s%l_==0 && nn==(s-1))
      return (-1);
    else if((s+1)%l_==0 && nn==(s+1))
      return (-1);
    else
      return nn;
  }


  //Vertical without pbc
  inline int V(int nn)const{
    if(nn>=nspins_)
      return -1;
    else if(nn<0)
      return -1;
    else
      return nn;
  }


};
