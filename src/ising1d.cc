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

//Transverse-field Ising model in 1d
class Ising1d{

  //number of spins
  const int nspins_;

  //value of the transverse field
  const double hfield_;

  //option to use periodic boundary conditions
  const bool pbc_;

  //pre-computed quantities
  std::vector<std::complex<double> > mel_;
  std::vector<std::vector<int> > flipsh_;

public:

  Ising1d(int nspins,double hfield,bool pbc=true):nspins_(nspins),hfield_(hfield),pbc_(pbc){
    Init();
  }

  void Init(){
    mel_.resize(nspins_+1,0.);
    flipsh_.resize(nspins_+1);

    for(int i=0;i<nspins_;i++){
      mel_[i+1]=-hfield_;
      flipsh_[i+1]=std::vector<int>(1,i);
    }
    std::cout<<"# Using the 1d Transverse-field Ising model with h = "<<hfield_<<std::endl;
  }

  //Finds the non-zero matrix elements of the hamiltonian
  //on the given state
  //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0
  //state' is encoded as the sequence of spin flips to be performed on state
  void FindConn(const std::vector<int> & state,std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){

    mel.resize(nspins_+1);
    flipsh.resize(nspins_+1);

    //assigning pre-computed matrix elements and spin flips
    mel=mel_;
    flipsh=flipsh_;

    //computing interaction part Sz*Sz
    mel[0]=0.;

    for(int i=0;i<(nspins_-1);i++){
      mel[0]-=double(state[i]*state[i+1]);
    }

    if(pbc_){
      mel[0]-=double(state[nspins_-1]*state[0]);
    }

  }

  int MinFlips()const{
    return 1;
  }

};
