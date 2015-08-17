// -*- C++ -*-

//System 
#include <map>
#include <algorithm>
#include <cmath>

// Rivet 
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
//#include "Rivet/RivetAIDA.hh"
// YODA
// #include "YODA/Histo2D.h"

//Projections
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

// fastjet contrib
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"

// Local Includes
#include "BOOSTFastJets.h"
using fastjet::contrib::Nsubjettiness;
using fastjet::contrib::Njettiness;
typedef std::map<std::string,Rivet::Histo1DPtr> BookedHistos;
template <typename lvec> static void dump4vec(lvec four_mom){
  std::cout<<"( "<<four_mom.pt()<<" [GeV], "<<four_mom.eta()<<", "<<four_mom.phi()<<", "<<four_mom.m()<<" [GeV])"<<std::endl;
}
/*
static void getAngle(const double& eta, const double& phi, 
		     const double& tmag, const double& ttheta,
		     double& theta){
  double etat=tmag*sin(ttheta);
  double phit=tmag*cos(ttheta);
  double mag=sqrt(pow(eta,2)+pow(phi,2));
  double x = (phit*phi + etat*eta)/(tmag*mag);
  theta=acos( x );
  return;
}
*/
namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  class MC_GENSTUDY_CHARMONIUM : public Analysis {
  public:

    /// Constructor
    MC_GENSTUDY_CHARMONIUM()
      : Analysis("MC_GENSTUDY_CHARMONIUM"),
	jetR(0.4),
	beta(1.0),
	nPtBins(10),
	binWidth(25)
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs(-2.5, 2.5, 500*MeV);
      addProjection(fs, "FS");
      const ChargedFinalState cfs(-2.5, 2.5, 500*MeV);
      addProjection(cfs, "CFS");
      IdentifiedFinalState muonfs(fs);
      muonfs.acceptIdPair(13);
      IdentifiedFinalState nufs(fs);
      nufs.acceptNeutrinos();
      MergedFinalState muAndNu(muonfs,nufs);
      addProjection(muAndNu,"MuonsAndNeutrinos");
      
      VetoedFinalState caloParts(FinalState(-4.2,4.2));
      caloParts.addVetoOnThisFinalState(muAndNu);
      
      ChargedLeptons lfs(fs);
      addProjection(lfs, "LFS");
      FastJets JetProjection(cfs,FastJets::ANTIKT, jetR);
      addProjection(JetProjection,"Jets");

      // Histograms
      _histograms["JetMult"] = bookHisto1D("JetMult",8,-0.5,8.5);

      _histograms["JPsiPt"] = bookHisto1D("JPsiPt" , 50, 0, 250);
      _histograms["JPsiM"] = bookHisto1D("JPsiM" , 50, 2.95, 3.2);
      _histograms["JPsiEta"] = bookHisto1D("JPsiEta" , 25, -4.2, 4.2);


      bookJetHistos("Jet");
      bookJetHistos("JetELo");
      bookJetHistos("JetEHi");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      cutFlow["Nominal"]++;
      const double weight = event.weight();
      applyProjection<FinalState>(event,"CFS");
      const ChargedLeptons& lProj = applyProjection<ChargedLeptons>(event, "LFS");
      if(lProj.chargedLeptons().empty()){
      	vetoEvent;
      }
      cutFlow["Leptons"]++;
      Particles muons;
      foreach(const Particle& lepton,lProj.chargedLeptons()){
      	if(abs(lepton.pid())==13){
      	  muons.push_back(lepton);
      	}
      }
      if(muons.size() < 2){
      	vetoEvent;
      }
      cutFlow["2Muons"]++;
      FourMomentum j_psi;
      find_j_psi(muons,j_psi);
      if(j_psi.mass()==0 ) { 
      	vetoEvent;
      }
      cutFlow["JPsiCand"]++;
      const FastJets& jetProj = applyProjection<FastJets>(event, "Jets");
      const PseudoJets jets = jetProj.pseudoJetsByPt(45*GeV);
      if(jets.size() < 1){
      	vetoEvent;
      }
      cutFlow["GreaterThan1Jets"]++;
      _histograms["JetMult"]->fill(jets.size(),weight);
      if(j_psi.pt() < 20*GeV){
      	vetoEvent;
      }
      cutFlow["JPsiPt"]++;
      //fill j_psi histos
      _histograms["JPsiEta"]->fill(j_psi.eta(),weight);
      _histograms["JPsiPt"]->fill(j_psi.pt(),weight);
      _histograms["JPsiM"]->fill(j_psi.mass(),weight);
      fastjet::PseudoJet charmJet;
      double delR(99.);
      double candDelR(99.);
      foreach(const fastjet::PseudoJet& j, jets){
      	delR=deltaR(FourMomentum(j.e(),j.px(),j.py(),j.pz()), j_psi);
      	if( delR < jetR && delR < candDelR) {
      	  charmJet=j;
      	  candDelR=delR;
      	}
      }
      if(isinf(deltaR(FourMomentum(charmJet.e(),
      				   charmJet.px(),
      				   charmJet.py(),
      				   charmJet.pz()),j_psi))){
      	vetoEvent;
      }
      cutFlow["charmJetMatch"]++;
      // fastjet::PseudoJet parton=jets.at(0);
      // if(parton.m2()==charmJet.m2()){
      // 	parton=jets.at(1);
      // }
      fillJetHistos("Jet",charmJet,j_psi,*jetProj.clusterSeq(),weight); 
      if(charmJet.pt() > 45*GeV && charmJet.pt() < 65*GeV){
      	fillJetHistos("JetELo",charmJet,j_psi,*jetProj.clusterSeq(),weight); 
      }
      else if(charmJet.pt() > 175*GeV){
      	fillJetHistos("JetEHi",charmJet,j_psi,*jetProj.clusterSeq(),weight); 
      }
    }

    /// Finalize
    void finalize() {
      cout<<"Cut flow"<<endl;
      cout<<"|-"<<endl;
      for(std::map<std::string, size_t>::const_iterator cut = cutFlow.begin();
	  cut != cutFlow.end(); ++cut){
	cout<<"| "<<cut->first << " | "<<cut->second<<" |"<<endl;
      }
      cout<<"|-"<<endl;
    }

    //@}


  private:
    void fillJetZ(const char* key,const int zVal, const double z,const double pt,const double weight){
      char histName[25];
      snprintf(histName,25,"%sPtZ%d",key,zVal);
      if(fuzzyEquals(z,zVal*0.1,0.1)){
	_histograms[histName]->fill(pt,weight);
      }
    }
    void bookJetHistos(const string& key){
      const double ptMax = (key=="Jet") ? 250. : 450; 
      _histograms[key+"Pt"]		 = bookHisto1D(key+"Pt" , 50, 0, ptMax);
      _histograms[key+"M"]		 = bookHisto1D(key+"M" , 50, 0, 40);
      _histograms[key+"Eta"]		 = bookHisto1D(key+"Eta" , 25, -4.2, 4.2);

      // _histograms[key+"PtZ3"]		 = bookHisto1D(key+"PtZ3" , 50, 0, ptMax);
      // _histograms[key+"PtZ5"]		 = bookHisto1D(key+"PtZ5" , 50, 0, ptMax);
      // _histograms[key+"PtZ8"]		 = bookHisto1D(key+"PtZ8" , 50, 0, ptMax);

      _histograms[key+"DeltaR"]		 = bookHisto1D(key+"DeltaR",50,0,jetR+0.1);
      _histograms[key+"Z"]		 = bookHisto1D(key+"Z",50,0,1.10);
      _histograms[key+"Dipolarity"]      = bookHisto1D(key+"Dipolarity" ,50,0.0,2);

      //Jet Pull to various axes
      _histograms[key+"PMag"]		 = bookHisto1D(key+"PMag" ,50,0,0.06);
      _histograms[key+"PTheta"]		 = bookHisto1D(key+"PTheta" ,50,-PI,PI);

      // _histograms[key+"PMagJPsi"]	 = bookHisto1D(key+"PMagJPsi" ,50,0,0.06);
      // _histograms[key+"PThetaJPsi"]      = bookHisto1D(key+"PThetaJPsi" ,50,0.,PI);

      // _histograms[key+"PMagPtn"]	 = bookHisto1D(key+"PMagPtn" ,50,0,0.06);
      // _histograms[key+"PThetaPtn"]       = bookHisto1D(key+"PThetaPtn" ,50,0.,PI);

      _histograms[key+"PtclMult"]	 = bookHisto1D(key+"PtclMult",41,-0.5,40.5);
      // N sub-jettiness
      _histograms[key+"NSJTau1"]	 = bookHisto1D(key+"NSJTau1", 40, -0.005, 1.005);
      _histograms[key+"NSJTau2"]	 = bookHisto1D(key+"NSJTau2", 40, -0.005, 1.005);
      _histograms[key+"NSJTau3"]	 = bookHisto1D(key+"NSJTau3", 40, -0.005, 1.005);

      _histograms[key+"NSJTau21"]	 = bookHisto1D(key+"NSJTau21", 40, -0.005, 1.25);
      _histograms[key+"NSJTau32"]	 = bookHisto1D(key+"NSJTau32", 40, -0.005, 1.25);

      // char histName[25];
      // for(int i=0; i < nPtBins; i++) {
      // 	//WARNING THIS IS COMPILER DEPENDENT!
      // 	sprintf(histName,(key+"Z_pt%d_%d").c_str(),int((ptMax/nPtBins)*i),int((ptMax/nPtBins)*(i+1)));
      // 	_histograms[string(histName)] = bookHisto1D(histName,50,0,1.10);
      // }
    }
    void fillJetHistos(const string& key, const fastjet::PseudoJet& jet, 
		       /* const PseudoJet& parton ,*/ const FourMomentum& j_psi,
		       const fastjet::ClusterSequence& clusterSeq, const double weight){
      // cout<<key<<endl;
      // dump4vec(jet);
      // dump4vec(parton);
	
      //cache variables for faster access.
      const double jet_pt(jet.pt());
      const double z(jet_pt > 0 ? j_psi.pt()/jet_pt : -1. );
      //kinematics
      _histograms[key+"DeltaR"]->fill(deltaR(j_psi,FourMomentum(jet.e(),jet.px(),jet.py(),jet.pz())),weight);
      _histograms[key+"Pt"]->fill(jet_pt,weight);
      _histograms[key+"M"]->fill(jet.m(),weight);
      _histograms[key+"Eta"]->fill(jet.eta(),weight);
      
      //substructure
      _histograms[key+"Z"]->fill(z,weight);
      _histograms[key+"PtclMult"]->fill(jet.constituents().size(),weight);
      std::vector<double> pull=JetPull(jet);
      _histograms[key+"PMag"]->fill(pull.at(0),weight);
      _histograms[key+"PTheta"]->fill(pull.at(1),weight);
      // if(pull.at(0)!=0){
      // 	double theta=-99.;
      // 	getAngle(j_psi.eta(),j_psi.phi(),pull.at(0),pull.at(1),theta);
      // 	_histograms[key+"PThetaJPsi"]->fill(theta,weight);
      // 	getAngle(parton.eta(),parton.phi(),pull.at(0),pull.at(1),theta);
      // 	_histograms[key+"PThetaPtn"]->fill(theta,weight);
      // }

      _histograms[key+"Dipolarity"]->fill(Dipolarity(jet),weight);
      vector<double> tau_vals(3,-1.);
      char histName[25];
      for(int i=1; i < 4; i++ ){
      	Nsubjettiness nSubJCalc(i,Njettiness::kt_axes, beta, jetR, jetR);
      	tau_vals[i-1]=nSubJCalc(jet);
      	snprintf(histName,25,"%sNSJTau%d",key.c_str(),i);
      	_histograms[histName]->fill(tau_vals[i-1],weight);
      }
      _histograms[key+"NSJTau21"]->fill(tau_vals[1] != 0 ? 
      					tau_vals[1]/tau_vals[0] : -1,weight);
      _histograms[key+"NSJTau32"]->fill(tau_vals[1] != 0 ? 
      					tau_vals[2]/tau_vals[1] : -1,weight);
    }
    void find_j_psi(Particles& muons,FourMomentum& j_psi) {
      FourMomentum cand;
      double deltaM=10000.;
      const double j_psi_m = 3.096916;
      foreach(const Particle& mu1, muons){
      	foreach(const Particle& mu2, muons){
      	  cand=mu1.momentum()+mu2.momentum();
      	  if(mu1.pid()*mu2.pid() < 0 && fabs(cand.mass()-j_psi_m) < deltaM ){
      	    j_psi=cand;
      	    deltaM=fabs(cand.mass()-j_psi_m);
      	  }
      	}
      }
      //reject if ~2*width of J/psi
      // if(deltaM > 0.2){
      // 	j_psi=FourMomentum(0,0,0,0);
      // }
      return;
    }
    const double jetR;
    const double beta;
    /// @name Histograms
    //@{
    BookedHistos _histograms;
    //@}
    std::map<std::string, size_t> cutFlow;
    const int nPtBins;
    const int binWidth;
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_GENSTUDY_CHARMONIUM);
}
