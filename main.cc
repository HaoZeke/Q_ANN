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

#include "src/nqs_paper.hh"

int main(int argc, char *argv[]){

  auto opts=ReadOptions(argc,argv);

  //Definining the neural-network wave-function
  Nqs wavef(opts["filename"]);

  int nsweeps=std::stod(opts["nsweeps"]);
  int nspins=wavef.Nspins();

  //Problem hamiltonian inferred from file name
  std::string model=opts["model"];

  bool printastes=opts.count("filestates");

  int seed=std::stoi(opts["seed"]);

  if(model=="Ising1d"){
    double hfield=std::stod(opts["hfield"]);
    Ising1d hamiltonian(nspins,hfield);

    //Defining and running the sampler
    Sampler<Nqs,Ising1d> sampler(wavef,hamiltonian,seed);
    if(printastes){
      sampler.SetFileStates(opts["filestates"]);
    }
    sampler.Run(nsweeps);
  }
  else if(model=="Heisenberg1d"){
    double jz=std::stod(opts["jz"]);
    Heisenberg1d hamiltonian(nspins,jz);

    //Defining and running the sampler
    Sampler<Nqs,Heisenberg1d> sampler(wavef,hamiltonian,seed);
    if(printastes){
      sampler.SetFileStates(opts["filestates"]);
    }
    sampler.Run(nsweeps);
  }
  else if(model=="Heisenberg2d"){
    double jz=std::stod(opts["jz"]);
    Heisenberg2d hamiltonian(nspins,jz);

    //Defining and running the sampler
    Sampler<Nqs,Heisenberg2d> sampler(wavef,hamiltonian,seed);
    if(printastes){
      sampler.SetFileStates(opts["filestates"]);
    }
    sampler.Run(nsweeps);
  }
  else{
    std::cerr<<"#The given input file does not correspond to one of the implemented problem hamiltonians";
    std::abort();
  }

}
