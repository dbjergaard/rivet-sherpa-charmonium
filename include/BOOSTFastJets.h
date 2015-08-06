//-*- C++ -*-

#ifndef RIVET_BOOSTFastJets_HH
#define RIVET_BOOSTFastJets_HH
#include "Rivet/Projections/FastJets.hh"
namespace Rivet{
  /// structs used in angular correlation calculations
  struct ACFparticlepair {
    double deltaR;
    double weight;
  };

  struct ACFpeak {
    double Rval;
    double height;
    double prominence;
    double partialmass;
    int index;
  };

  struct ppsortfunction { bool operator()
  (ACFparticlepair a, ACFparticlepair b) { return (a.deltaR < b.deltaR); } };

  struct ACFfunction {
    double Rval;
    double value;
  };
  /// Calculate Dipolarity of Jet
  double Dipolarity(const fastjet::PseudoJet &j);
  // Calculate Pull of Jet
  ///\vec{t} ===\Sum_{i\in J} |r_i|p_{Ti}/p_{TJ}\vec{r_i}
  std::vector<double> JetPull(const fastjet::PseudoJet &jet, const double ptmin=-1*GeV); 
  std::vector<double> JetPull(const fastjet::PseudoJet &j1, const fastjet::PseudoJet &j2,
			      const double ptmin=-1*GeV);
  /// Calculate JetCharge
  double JetCharge(const Jet& jet, const fastjet::PseudoJet &j, const double k, const double ptmin=-1.);

  fastjet::JetAlgorithm setJetAlgorithm(FastJets::JetAlgName subJetAlgorithm);
  /// Create a filter, run it over specified jet
  /// Butterworth, Davison, Rubin and Salam, arXiv:0802.2470
  fastjet::PseudoJet Filter(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			    int hardest,double subjet_R);

  /// Create a trimmer, run it over specified jet
  /// Krohn, Thaler and Wang, arXiv:0912.1342
  fastjet::PseudoJet Trimmer(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			     double percentage, double subjet_R);

  /// Create a pruner, run it over specified jet
  /// Ellis, Vermilion and Walsh, arXiv:0903.5081
  fastjet::PseudoJet Pruner(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			    double zcut, double Rcut_factor);

  /// Get N=n_jets subjets to be used for finding N-subjettiness
  /// Thaler, Van Tilburg, arXiv:1011.2268
  PseudoJets GetAxes(unsigned int n_jets,
		     PseudoJets& inputJets, FastJets::JetAlgName subjet_def, double subR);

  /// Get the N-subjettiness with respect to the subjet axes.
  /// Thaler, Van Tilburg, arXiv:1011.2268
  double TauValue(double beta, double jet_rad,
		  PseudoJets& particles, PseudoJets& axes);

  /// Update axes towards Tau(y, phi) minimum.
  /// Thaler, Van Tilburg, arxiv:1108.2701
  void UpdateAxes(double beta,
		  PseudoJets& particles, PseudoJets& axes);

  /// Find peaks in Angular Structure Function for the given particles
  /// Jankowiak, Larkowski, arxiv:1104.1646
  /// Based on code by Jankowiak and Larkowski
  /// Normalisation: 0 - error function, else step function
  vector<ACFpeak> ASFPeaks(PseudoJets& particles,
			   unsigned int most_prominent = 0, double minprominence = 0.0,
			   double sigma = 0.06, unsigned int meshsize = 500, unsigned int normalisation = 0);

  /// Return vectors with R values ([0]), unnormalised ASF ([1]),
  /// and the normalisation function depending on the normalisation variable ([2])
  /// Mainly used for average angular structure function.
  /// Jankowiak, Larkowski, arXiv:1201.2688
  vector<vector<double> > ASF(PseudoJets& particles,
                double sigma = 0.06, unsigned int meshsize = 500, unsigned int normalisation = 0);

  double KeyColToRight(int p, vector<ACFpeak> peaks, vector<double> ASF_erf);
  double KeyColToLeft(int p, vector<ACFpeak> peaks, vector<double> ASF_erf);
}
#endif
