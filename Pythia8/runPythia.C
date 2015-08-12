#include "Pythia8/Pythia.h"
#ifdef __PYTHIA_8210__
#include "Pythia8Plugins/HepMC2.h"
#else
#include "Pythia8/Pythia8ToHepMC.h"
#endif

#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <iostream>

#include <map> 

using namespace Pythia8; 
template < class T> bool contains(const std::vector<T>& container, const T& value){
  return std::find(container.begin(),container.end(),value)!=container.end();
}

const Particle& find_ptcl(const int pid,const Event& evt){
  return evt[pid];
}
void print_family_tree(const Particle& ptcl, const Event& evt){
  vector<int> daughters = ptcl.daughterList();
  bool isCopy = (ptcl.daughter1() == ptcl.daughter2() && ptcl.daughter2() > 0);
  if(!isCopy){
    cout <<" ("<<ptcl.id();
  }
  if(daughters.size() > 0){
    for(vector<int>::const_iterator d=daughters.begin(); d!=daughters.end(); ++d ){
      print_family_tree(find_ptcl(*d, evt),evt);
    }
  }
  if(!isCopy){
    cout <<")";
  }
  return;
}

int main(int argc, char* argv[]) {
  char tune[25];
  char seed[25];

  //tune = "5";
  //seed = "1370";
  strncpy(tune,"5",25);
  strncpy(seed,"1370",25);
  //TODO this has gotten big enough to require real command line
  //options, setup and option parser and do it right.  Someone's going
  //to get headaches over this very soon

  //TODO add seeds back as an option, required for running on condor!
  if(argc < 3)
    {
      cerr<<"Usage: "<<argv[0]<<" [config file] [out file] [conf lines ...]"<<endl;
      return 1;
    }
  ifstream is(argv[1]);  
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }
  cout<<"Configuring PYTHIA from "<<argv[1]<<endl;
  cout<<"Writing output to "<<argv[2]<<endl;
  vector<std::string> commandStrings;
  for(int i=3; i < argc; i++){
    commandStrings.push_back(string(argv[i]));
    cout <<" Got Pythia command: "<<commandStrings.back()<<endl;
  }

  cout <<"Warning, command line arguments aren't type-checked, don't be stupid." <<endl;
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);
  char processline[128];
  Pythia pythia;
  pythia.readFile(argv[1]);
  int    nRequested    = 100;
  nRequested = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");
  // Set the seed
  pythia.readString("Random:setSeed = on");
  sprintf(processline,"Random:seed = %s",seed);
  pythia.readString(processline);
  // Set the tune
  sprintf(processline,"Tune:pp = %s",tune);
  pythia.readString(processline);

  for(std::vector<std::string>::const_iterator cmd = commandStrings.begin(); cmd != commandStrings.end(); ++cmd){
    pythia.readString(*cmd);
  }
  pythia.init();
  
  int iAbort = 0;
  int nGenerated=0;
  std::vector<int> requestedPdgId;//(3,0);
  while(nGenerated < nRequested) {
    requestedPdgId.clear();
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) {
	cout << " Aborted since reached end of Les Houches Event File\n"; 
	break; 
      }
	  
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
    for(int i=0; i < pythia.event.size(); i++){
      if(abs(pythia.event[i].id())==443 || abs(pythia.event[i].id())==13){
	requestedPdgId.push_back(pythia.event[i].id());
      }
    }
    if(!(contains(requestedPdgId, 443) && 
	 contains(requestedPdgId, 13)  &&
	 contains(requestedPdgId, -13))){
      continue;
    }
    /*
    for(int i=0; i < pythia.event.size(); i++){
      if(abs(pythia.event[i].id())==9940003){
	print_family_tree(pythia.event[i],pythia.event);
	cout << endl  << endl;
	break;
      }
    }
    */
    // print_family_tree( pythia.event[1], pythia.event);

    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    ToHepMC.fill_next_event( pythia, hepmcevt );
    ascii_io << hepmcevt ;
    delete hepmcevt;
    nGenerated++;
  }

  //pythia.stat();
  return 0;
}
