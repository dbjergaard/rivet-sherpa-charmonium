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
#include "Rivet/Projections/UnstableFinalState.hh"
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
namespace Rivet {
  /// Generic analysis looking at various distributions of final state particles
  class MC_GENSTUDY_JPSI_VEC : public Analysis {
  public:

    /// Constructor
    MC_GENSTUDY_JPSI_VEC()
      : Analysis("MC_GENSTUDY_JPSI_VEC"),
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
      addProjection(UnstableFinalState(), "UFS");
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
      // _histograms["JetMult"] = bookHisto1D("JetMult",8,-0.5,8.5);

      _histograms["JPsiPt"] = bookHisto1D("JPsiPt" , 50, 0, 35);
      _histograms["JPsiM"] = bookHisto1D("JPsiM" , 50, 2.95, 3.2);
      _histograms["JPsiEta"] = bookHisto1D("JPsiEta" , 25, -4.2, 4.2);
      _histograms["JPsiRap"] = bookHisto1D("JPsiRap" , 25, -4.2, 4.2);
      _histograms["PdgId"] = bookHisto1D("PdgId",5,0,4);

      pdgIDMap[443]=0;
      pdgIDMap[10441]=1;
      pdgIDMap[20443]=2;
      pdgIDMap[445]=3;
      pdgIDMap[100443]=4;
      // bookJetHistos("Jet");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      cutFlow["Nominal"]++;
      const double weight = event.weight();
      applyProjection<FinalState>(event,"CFS");
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");
      //fill j_psi histos
      foreach (const Particle& j_psi, ufs.particles()) {
        if (!(j_psi.abspid() == 443 || // J\psi
	      j_psi.abspid() == 10441 || //Chi_0
	      j_psi.abspid() == 20443 || //Chi_1
	      j_psi.abspid() == 445 ||  //Chi_2
	      j_psi.abspid() == 100443)) continue; //Psi(2S)
	_histograms["PdgId"]->fill(pdgIDMap[j_psi.abspid()]);
	_histograms["JPsiRap"]->fill(j_psi.rapidity(),weight);
	_histograms["JPsiEta"]->fill(j_psi.eta(),weight);
	_histograms["JPsiPt"]->fill(j_psi.pt(),weight);
	_histograms["JPsiM"]->fill(j_psi.mass(),weight);
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

      _histograms[key+"DeltaR"]		 = bookHisto1D(key+"DeltaR",50,0,jetR+0.1);
      _histograms[key+"Z"]		 = bookHisto1D(key+"Z",50,0,1.10);

      _histograms[key+"PtclMult"]	 = bookHisto1D(key+"PtclMult",41,-0.5,40.5);
      // N sub-jettiness
      _histograms[key+"NSJTau1"]	 = bookHisto1D(key+"NSJTau1", 40, -0.005, 1.005);
      _histograms[key+"NSJTau2"]	 = bookHisto1D(key+"NSJTau2", 40, -0.005, 1.005);
      _histograms[key+"NSJTau3"]	 = bookHisto1D(key+"NSJTau3", 40, -0.005, 1.005);

      _histograms[key+"NSJTau21"]	 = bookHisto1D(key+"NSJTau21", 40, -0.005, 1.25);
      _histograms[key+"NSJTau32"]	 = bookHisto1D(key+"NSJTau32", 40, -0.005, 1.25);

    }
    void fillJetHistos(const string& key, const fastjet::PseudoJet& jet, const FourMomentum& j_psi,
		       const fastjet::ClusterSequence& clusterSeq, const double weight){
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
      return;
    }
    const double jetR;
    const double beta;
    /// @name Histograms
    //@{
    BookedHistos _histograms;
    //@}
    std::map<std::string, size_t> cutFlow;
    std::map<int, int> pdgIDMap;
    const int nPtBins;
    const int binWidth;
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_GENSTUDY_JPSI_VEC);
}
