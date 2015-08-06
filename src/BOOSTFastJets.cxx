#include "BOOSTFastJets.h"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {
  /// D===1/R_{12} \Sum_{i\in J} p_{Ti}/p_{TJ} R_i
  double Dipolarity(const fastjet::PseudoJet &j) {
    fastjet::PseudoJet jet1,jet2;
    if (not (j.has_parents(jet1,jet2))) return -1;//not ideal axes
    double dipolarity(0.0);
    double sumpt(0.0);

    double deta(jet2.eta() - jet1.eta()), dphi(mapAngleMPiToPi(jet2.phi() - jet2.phi()));  //vector from 1 to 2.
    const double dmag2= deta*deta + dphi*dphi;
    if (dmag2 < 1e-3) return -1;         //no resolution
    const double dmag = sqrt(dmag2);
    deta /= dmag;                        //now it's a unit vector from 1 to 2
    dphi /= dmag;
    double vx,vy,pt,project;
    foreach (const fastjet::PseudoJet & c, j.constituents()) {
      pt = c.perp();
      sumpt += pt;
      vx = c.eta() - jet1.eta();
      vy = mapAngleMPiToPi(c.phi()-jet1.phi());
      project = vx*deta + vy*dphi;
      if (((project > 0) && (project < dmag))) { //nearest distance to segment is perp. projection
	dipolarity += pt * pow(vx*dphi - vy*deta,2);
      } else {
	if (project > 0) { //closer to jet2, so move the origin
	  vx = c.eta() - jet2.eta();
	  vy = mapAngleMPiToPi(c.phi()-jet2.phi());
	}
	dipolarity += pt * (vx*vx + vy*vy);  //nearest distance is radial vector to origin
      }
    }//constit loop
    if (sumpt < 1e-3) return -1;
    return dipolarity/(sumpt*dmag2);
  }

  ///\vec{t} ===\Sum_{i\in J} |r_i|p_{Ti}/p_{TJ}\vec{r_i}
  std::vector<double> JetPull(const fastjet::PseudoJet &jet, 
			      const double ptmin) {
    return JetPull(jet,jet,ptmin);
  }
  std::vector<double> JetPull(const fastjet::PseudoJet &j1, 
			      const fastjet::PseudoJet &j2,
			      const double ptmin) {
    std::vector<double> pull(2,-99.);
    if(!j1.has_valid_cluster_sequence()){
      return pull;
    }
    const PseudoJets parts = j1.constituents();
    const double jetRap = j2.rapidity(), jetPhi = j2.phi();
    double ty=0, tphi=0, tmag=0, ttheta=0, dphi=0;
    if(parts.size() > 1) {
      foreach (const fastjet::PseudoJet& p, parts) {
	//don't generate a large pull for jets at 2pi
	dphi = mapAngleMPiToPi(p.phi()-jetPhi); 
	if(p.pt() > ptmin) { 
	  double ptTimesRmag=sqrt(pow(p.rapidity()-jetRap,2) + pow(dphi,2))*p.pt();
	  ty+=ptTimesRmag*(p.rapidity()-jetRap);
	  tphi+=ptTimesRmag*(dphi);
	}
      }
      //now parametrize \vec{t}=|t|(cos(\theta_t),sin(\theta_t))
      tmag=sqrt(pow(ty,2) + pow(tphi,2))/j2.pt();
      if(tmag>0) {
	ttheta=atan2(tphi,ty);
      }
    }
    pull[0]=tmag;
    pull[1]=ttheta;
    return pull;
  }

  ///Q===\Sum_{i\in J} q_i*p_{Ti}^k/p_{TJ}
  double JetCharge(const Jet& jet, const fastjet::PseudoJet &j, const double k, const double ptmin) {
    double q(0);
    foreach (const Particle& p, jet.particles()) {
      if(p.pt() < ptmin) continue; 
      q += p.charge() * pow(p.pt(),k);
    }
    return q/pow(j.pt(),k);
  }

  fastjet::JetAlgorithm setJetAlgorithm(FastJets::JetAlgName subJetAlgorithm)
  {
    //Do we want to support all enums? This is only a subset...
    switch(subJetAlgorithm) {
    case FastJets::KT:
      return fastjet::kt_algorithm;
    case FastJets::ANTIKT:
      return fastjet::antikt_algorithm;
    case FastJets::CAM:
      return fastjet::cambridge_algorithm;
    case FastJets::DURHAM:
      return fastjet::ee_kt_algorithm;
    default:
      return fastjet::undefined_jet_algorithm;
    }
  }

  fastjet::PseudoJet Filter(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			    int hardest,double subjet_R=0.3)  {
    //assert(clusterSeq);
    //sanity check on the jet
    if (jet.E() <= 0.0 || clusterSeq.constituents(jet).size() == 0) {
      return jet;
    }
    fastjet::Filter filter(fastjet::JetDefinition(setJetAlgorithm(subjet_def), subjet_R), fastjet::SelectorNHardest(hardest));
    return filter(jet);
  }

  fastjet::PseudoJet Trimmer(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			     double percentage, double subjet_R=0.3) {
    //assert(clusterSeq);
    //sanity check on the jet
    if (jet.E() <= 0.0 || clusterSeq.constituents(jet).size() == 0) {
      return jet;
    }
    fastjet::Filter filter(fastjet::JetDefinition(setJetAlgorithm(subjet_def), subjet_R), fastjet::SelectorPtFractionMin(percentage));
    return filter(jet);
  }
  fastjet::PseudoJet Pruner(const fastjet::ClusterSequence& clusterSeq, fastjet::PseudoJet jet, FastJets::JetAlgName subjet_def,
			    double zcut=0.1, double Rcut_factor=0.5) {
    //sanity check on the jet
    //assert(clusterSeq);
    if (jet.E() <= 0.0 || clusterSeq.constituents(jet).size() == 0) {
      return jet;
    }
    //NOTE: Pruner sets the jet algorithm R to the maximum allowed, so setting it
    //to 1 here only to follow jetalg syntax
    //finally declare the pruner and apply it to the jet
    fastjet::Pruner pruner(fastjet::JetDefinition(setJetAlgorithm(subjet_def), 1), zcut, Rcut_factor);
    return pruner(jet);
  }

  PseudoJets GetAxes(unsigned int n_jets, PseudoJets& inputJets, FastJets::JetAlgName subjet_def, double subR) {
    //sanity check
    if (inputJets.size() < n_jets) {
      std::cout << "Not enough input particles." << endl;
      return inputJets;
    }
    //get subjets, return
    fastjet::ClusterSequence sub_clust_seq(inputJets, fastjet::JetDefinition(setJetAlgorithm(subjet_def), subR));
    return sub_clust_seq.exclusive_jets((signed)n_jets);
  }

  double TauValue(double beta, double jet_rad,
		  PseudoJets& particles, PseudoJets& axes) {
    double tauNum = 0.0;
    double tauDen = 0.0;
    if(particles.size() == 0)return 0.0;
    for (unsigned int i = 0; i < particles.size(); i++) {
      // find minimum distance (set R large to begin)
      double minR = 10000.0;
      for (unsigned int j = 0; j < axes.size(); j++) {
	double tempR = sqrt(particles[i].squared_distance(axes[j]));
	if (tempR < minR) minR = tempR;
      }
      //calculate nominator and denominator
      tauNum += particles[i].perp() * pow(minR,beta);
      tauDen += particles[i].perp() * pow(jet_rad,beta);
    }
    //return N-subjettiness
    return tauNum/tauDen;
  }

  void UpdateAxes(double beta,
		  PseudoJets& particles, PseudoJets& axes) {
    vector<int> belongsto;
    //no reason not to use foreach here
    for (unsigned int i = 0; i < particles.size(); i++) {
      // find minimum distance axis
      int assign = 0;
      double minR = 10000.0;
      for (unsigned int j = 0; j < axes.size(); j++) {
	double tempR = sqrt(particles[i].squared_distance(axes[j]));
	if (tempR < minR) {
	  minR = tempR;
	  assign = j;
	}
      }
      belongsto.push_back(assign);
    }
    // iterative step
    double deltaR2, distphi;
    vector<double> ynom, phinom, den;
    ynom.resize(axes.size());
    phinom.resize(axes.size());
    den.resize(axes.size());

    for (unsigned int i = 0; i < particles.size(); i++) {
      distphi = particles[i].phi() - axes[belongsto[i]].phi();
      deltaR2 = particles[i].squared_distance(axes[belongsto[i]]);

      if(fuzzyEquals(deltaR2, 0.,1e-9))continue;

      if (abs(distphi) <= M_PI) phinom.at(belongsto[i]) += particles[i].perp() * particles[i].phi() * pow(deltaR2, (beta-2)/2);
      else if ( distphi > M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (-2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
      else if ( distphi < M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (+2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);

      ynom.at(belongsto[i]) += particles[i].perp() * particles[i].rap() * pow(deltaR2, (beta-2)/2);

      den.at(belongsto[i]) += particles[i].perp() * pow(deltaR2, (beta-2)/2);
    }

    // reset to new axes

    for (unsigned int j = 0; j < axes.size(); j++) {

      if (fuzzyEquals(den[j], 0.,1e-9)) axes.at(j) = axes[j];
      else {
	axes.at(j).reset_momentum_PtYPhiM(axes[j].perp(), ynom[j] / den[j], fmod( 2*M_PI + (phinom[j] / den[j]), 2*M_PI ), axes[j].perp()/2);
      }
    }
  }

  double KeyColToRight(int p, vector<ACFpeak> peaks, vector<double> ASF_erf) {
    int higherpeak = -1;
    double height = peaks[p].height;
    double keycol = height;
    for (unsigned int k = p+1; k < peaks.size(); k++) {
      if (peaks[k].height > height) {
	higherpeak = k;
	break;
      }
    }
    if (higherpeak != -1) {

      int startindex = peaks[p].index;
      int endindex = peaks[higherpeak].index;
      for (int j = startindex+1; j < endindex; j++) {
	if (ASF_erf[j] < keycol) keycol = ASF_erf[j];
      }
    }
    else keycol = 0.;
    return keycol;
  }

  double KeyColToLeft(int p, vector<ACFpeak> peaks, vector<double> ASF_erf) {
    int higherpeak = -1;
    double height = peaks[p].height;
    double keycol = height;
    for (int k = p-1; k >= 0; k--) {

      if (peaks[k].height > height) {
	higherpeak = k;
	break;
      }
    }
    if (higherpeak != -1) {
      int endindex = peaks[p].index;
      int startindex = peaks[higherpeak].index;
      for (int j = startindex+1; j < endindex; j++) {
	if (ASF_erf[j] < keycol) keycol = ASF_erf[j];
      }
    }
    else keycol = 0.;
    return keycol;
  }

  vector<ACFpeak> ASFPeaks(PseudoJets& particles,
			   unsigned int most_prominent, double minprominence,
			   double sigma, unsigned int meshsize, unsigned int normalisation) {
    vector<ACFpeak> peaks;

    //sanity check
    if(particles.size() < 2) {
      cout << "Not enough particles in jet for ACF." << endl;
      return peaks;
    }
    //pair all particles up
    vector<ACFparticlepair> pairs;
    ACFparticlepair dummy;
    for(unsigned int k = 0; k < particles.size(); k++) {
      for(unsigned int j = 0; j < k; j++) {
	double phidist = abs(particles[k].phi_std() - particles[j].phi_std());
	if( phidist > M_PI ) phidist = 2 * M_PI - phidist;
	//dummy.deltaR = sqrt( pow(particles[k].pseudorapidity() - particles[j].pseudorapidity(), 2) +
	//	pow(phidist, 2) );
	dummy.deltaR = sqrt(particles[k].plain_distance(particles[j]));
	dummy.weight = particles[k].perp() * particles[j].perp() * dummy.deltaR * dummy.deltaR;
	pairs.push_back(dummy);
      }
    }

    //sort by delta R
    sort(pairs.begin(), pairs.end(), ppsortfunction());

    double Rmax = pairs[pairs.size() - 1].deltaR;

    double rVal = 0.;
    double xArg = 0.;
    vector<double> ACF(meshsize), ASF_gauss(meshsize), Rvals(meshsize), erf_denom(meshsize), gauss_peak(meshsize);
    ACF[0] = 0.;
    ASF_gauss[0] = 0.;
    Rvals[0] = 0.;
    erf_denom[0] = 0.;
    gauss_peak[0] = 0.;
    double fVal = 0.;
    double eVal = 0.;
    double gVal = 0.;

    //mesh loop
    for (unsigned int k = 1; k < meshsize; k++) {

      rVal = (double)k*Rmax/(meshsize-1);
      Rvals[k] = rVal;

      //reset
      fVal = 0.;
      eVal = 0.;
      gVal = 0.;

      //Loop on pairs.
      for (unsigned int j = 0; j < pairs.size(); j++) {
	//ACF-Add pairs within mesh's deltaR.
	if (pairs[j].deltaR <= rVal) fVal += pairs[j].weight;

	//Smoothing function argument
	xArg = (double)(rVal-pairs[j].deltaR)/sigma;

	//ASF Error Function Denominator: Add pairs weighted by Erf values.
	eVal += (double)pairs[j].weight*0.5*(1.0+erf(xArg));

	//ASF Gaussian Numerator: Add pairs weight by Gaussian values.
	gVal += (double)pairs[j].weight * (exp( -pow ((abs(xArg)),2.0)));

      }//end pair loop
      ACF[k] = fVal;
      erf_denom[k] = eVal;
      gauss_peak[k] = gVal;
      ASF_gauss[k] = (double)gVal*(1/sqrt(M_PI))*rVal/sigma; //Normalized Gaussian value

    }//end mesh loop

    vector<double> ASF_erf(meshsize), ASF(meshsize);
    ASF[0] = 0.;
    ASF[meshsize-1] = 0.;
    ASF_erf[0] = 0.;

    //Total jet mass
    double jetmass = ACF[meshsize-1];

    //Second mesh loop
    for (unsigned int k = 1; k < meshsize; k++) {

      //Compute simple (spiked) ASF
      if (ACF[k-1] == 0.) {
	ASF[k] = 0.;
      }
      else if (k < meshsize - 2)
        {
	  ASF[k] = (log(ACF[k+1]) - log(jetmass*ACF[k-1]))/(log(Rvals[k+1]) - log(Rvals[k-1]));
        }

      //Compute gaussian (smoothed) ASF
      if(normalisation == 0)ASF_erf[k] = (fuzzyEquals(erf_denom[k],0.,1e-9)) ? 0. : ASF_gauss[k]/erf_denom[k];
      else ASF_erf[k] = (fuzzyEquals(ACF[k],0.,1e-9)) ? 0. : ASF_gauss[k]/ACF[k];
      //if(erf_denom[k] == 0.0)ASF_erf[k] = 0.0;
      //else ASF_erf[k] = ASF_gauss[k]/erf_denom[k];

      //Normalize ACF
      ACF[k] /= jetmass;

    }//end mesh loop

    ACFpeak myPeak;
    double lefth =  0.;
    double height = 0.;
    double righth = 0.;

    for (unsigned int k = 1; k < meshsize-1; k++) {
      lefth  = ASF_erf[k-1];
      height = ASF_erf[k];
      righth = ASF_erf[k+1];
      //Found a peak?
      if (lefth < height && height > righth) {
	myPeak.height = height;
	myPeak.Rval   = Rvals[k];
	myPeak.index  = k;
	//Peaks are stored according to
	//index (equilvalent to R values)
	//from lowest to highest
	peaks.push_back(myPeak);
      }
    }
    //Partial mass of peak
    for (unsigned int p = 0; p < peaks.size(); p++)peaks[p].partialmass = (double)sqrt(gauss_peak[peaks[p].index]);
    //peaks[p].partialmass = (double)sqrt(sqrt(M_PI)*sigma*peaks[p].height*ACF[peaks[p].index]*jetmass/peaks[p].Rval);
    //Prominence of peak
    for (unsigned int k = 0; k < peaks.size(); k++) {
      double height = peaks[k].height;
      double leftdescent = height - KeyColToLeft(k, peaks, ASF_erf);
      double rightdescent = height - KeyColToRight(k, peaks, ASF_erf);
      if (leftdescent < rightdescent) peaks[k].prominence = leftdescent;
      else peaks[k].prominence = rightdescent;
    }
    //return all peaks
    if(most_prominent == 0 && fuzzyEquals(minprominence, 0.,1e-9) ) {
      return peaks;
    }
    //return all peaks over given prominence
    else if(most_prominent == 0) {
      vector<ACFpeak> dummyp;
      for(unsigned int i = 0; i < peaks.size(); i++) {
	if(peaks[i].prominence > minprominence) {
	  dummyp.push_back(peaks[i]);
	}
      }
      return dummyp;
    }

    //else return the N most_prominent ones (might return less than N peaks if there aren't enough over minprominence)
    else {
      vector<ACFpeak> dummyp;
      for(unsigned int i = 0; i < peaks.size(); i++) {
	if(peaks[i].prominence > minprominence) {
	  dummyp.push_back(peaks[i]);
	}
      }

      if(dummyp.size() > most_prominent) {
	vector<ACFpeak> newdummy;
	for(unsigned int j = 0; j < most_prominent; j++) {
	  int assign = -1;
	  double prom = 0.;
	  for(unsigned int k = 0; k < dummyp.size(); k++) {
	    if(dummyp[k].prominence > prom) {
	      assign = k;
	      prom = dummyp[k].prominence;
	    }
	  }
	  if(assign != -1)newdummy.push_back(dummyp[assign]);
	  dummyp[assign].prominence = 0.;
	}
	return newdummy;
      }
      else return dummyp;
    }
  }

  vector<vector<double> > ASF(PseudoJets& particles,
			      double sigma, unsigned int meshsize, unsigned int normalisation) {

    vector<vector<double> > functions;

    //sanity check
    if(particles.size() < 2) {
      cout << "Not enough particles in jet for ACF." << endl;
      return functions;
    }
    //pair all particles up
    vector<ACFparticlepair> pairs;
    ACFparticlepair dummy;
    for(unsigned int k = 0; k < particles.size(); k++) {
      for(unsigned int j = 0; j < k; j++) {
	double phidist = abs(particles[k].phi_std() - particles[j].phi_std());
	if( phidist > M_PI ) phidist = 2 * M_PI - phidist;
	//dummy.deltaR = sqrt( pow(particles[k].pseudorapidity() - particles[j].pseudorapidity(), 2) +
	//	pow(phidist, 2) );
	dummy.deltaR = sqrt(particles[k].plain_distance(particles[j]));
	dummy.weight = particles[k].perp() * particles[j].perp() * dummy.deltaR * dummy.deltaR;
	pairs.push_back(dummy);
      }
    }

    //sort by delta R
    sort(pairs.begin(), pairs.end(), ppsortfunction());

    double Rmax = pairs[pairs.size() - 1].deltaR;

    double rVal = 0.;
    double xArg = 0.;
    vector<double> ACF(meshsize), ASF_gauss(meshsize), Rvals(meshsize), erf_denom(meshsize), gauss_peak(meshsize);
    ACF[0] = 0.;
    ASF_gauss[0] = 0.;
    Rvals[0] = 0.;
    erf_denom[0] = 0.;
    gauss_peak[0] = 0.;
    double fVal = 0.;
    double eVal = 0.;
    double gVal = 0.;

    //mesh loop
    for (unsigned int k = 1; k < meshsize; k++) {

      rVal = (double)k*Rmax/(meshsize-1);
      Rvals[k] = rVal;

      //reset
      fVal = 0.;
      eVal = 0.;
      gVal = 0.;

      //Loop on pairs.
      for (unsigned int j = 0; j < pairs.size(); j++) {
	//ACF-Add pairs within mesh's deltaR.
	if (pairs[j].deltaR <= rVal) fVal += pairs[j].weight;

	//Smoothing function argument
	xArg = (double)(rVal-pairs[j].deltaR)/sigma;

	//ASF Error Function Denominator: Add pairs weighted by Erf values.
	eVal += (double)pairs[j].weight*0.5*(1.0+erf(xArg));

	//ASF Gaussian Numerator: Add pairs weight by Gaussian values.
	gVal += (double)pairs[j].weight * (exp( -pow ((abs(xArg)),2.0)));

      }//end pair loop
      ACF[k] = fVal;
      erf_denom[k] = eVal;
      gauss_peak[k] = gVal;
      ASF_gauss[k] = (double)gVal*(1/sqrt(M_PI))*rVal/sigma; //Normalized Gaussian value

    }//end mesh loop
    functions.push_back(Rvals);
    functions.push_back(ASF_gauss);
    if(normalisation == 0) functions.push_back(erf_denom);
    else functions.push_back(ACF);

    return functions;
  }

}
