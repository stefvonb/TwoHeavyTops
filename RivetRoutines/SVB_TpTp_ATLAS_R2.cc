// -*- C++ -*-
// SVB_TpTp_ATLAS_R2.cc
// Created by Stefan von Buddenbrock on 2016/03/10.

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "fastjet/tools/Filter.hh"

#define M_Z 91.2*GeV
#define M_h 125.0*GeV

namespace Rivet {

  class SVB_TpTp_ATLAS_R2 : public Analysis {
  public:

    /// Constructor
    SVB_TpTp_ATLAS_R2()
    : Analysis("SVB_TpTp_ATLAS_R2"),
    _trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)))
    {    }

    /// Book histograms and initialise projections before the run
    void init() {

      /// Open cut final state for finding photon 'underpants'
      FinalState FS(Cuts::open());
      declare(FS, "FS");

      // Project photons for dressing
      IdentifiedFinalState ph_dressing_FS(FS);
      ph_dressing_FS.acceptIdPair(PID::PHOTON);

      // Project bare electrons
      IdentifiedFinalState el_bare_FS(FS);
      el_bare_FS.acceptIdPair(PID::ELECTRON);

      // Project dressed electrons
      DressedLeptons el_dressed_FS(ph_dressing_FS, el_bare_FS, 0.1, Cuts::pT > 7/GeV && Cuts::abseta < 2.47, true, false);
      declare(el_dressed_FS,"EL_DRESSED_FS");

      // Project bare muons
      IdentifiedFinalState mu_bare_FS(FS);
      mu_bare_FS.acceptIdPair(PID::MUON);

      // Project dressed muons
      DressedLeptons mu_dressed_FS(ph_dressing_FS, mu_bare_FS, 0.1, Cuts::pT > 6/GeV && Cuts::abseta < 2.7, true, false);
      declare(mu_dressed_FS,"MU_DRESSED_FS");

      // Build the anti-kT R=0.4 jets, using FinalState particles (vetoing neutrinos)
      FastJets antikt_04_jets(FS, FastJets::ANTIKT, 0.4, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      declare(antikt_04_jets, "ANTIKT_04_JETS");

      TauFinder had_tau_FS(TauFinder::HADRONIC, Cuts::pT > 20.0/GeV && Cuts::abseta < 2.5);
      addProjection(had_tau_FS,"HAD_TAU_FS");

      ChargedFinalState TRACKS_FS(Cuts::pT > 1.0/GeV);
      declare(TRACKS_FS, "TRACKS_FS");

      _c_passed_SR = bookCounter("passed_SR", "passed_SR");
      _c_passed_presel = bookCounter("passed_presel", "passed_presel");

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Particles& baseline_electrons = apply<DressedLeptons>(event, "EL_DRESSED_FS").particlesByPt();
      const Particles& baseline_muons = apply<DressedLeptons>(event, "MU_DRESSED_FS").particlesByPt();
      const Particles& all_taus = apply<TauFinder>(event, "HAD_TAU_FS").particlesByPt();

      const Particles& tracks = apply<ChargedFinalState>(event, "TRACKS_FS").particlesByPt();

      const Jets& baseline_small_R_jets = apply<FastJets>(event, "ANTIKT_04_JETS").jetsByPt(Cuts::pT>20.0/GeV && Cuts::absrap < 4.4);

      FourMomentum cal_invisible(0.0, 0.0, 0.0, 0.0);
      foreach (const Jet& j, baseline_small_R_jets) {cal_invisible -= j.momentum();}
      foreach (const Particle& e, baseline_electrons) {cal_invisible -= e.momentum();}
      foreach (const Particle& m, baseline_muons) {cal_invisible -= m.momentum();}
      double ETmiss = cal_invisible.pT();

      Particles baseline_taus;
      foreach (const Particle& t, all_taus) {
        unsigned int n_associated_tracks = 0;
        foreach (const Particle& track, tracks) {
          if (deltaR(t, track) < 0.2) n_associated_tracks++;
        }
        if (n_associated_tracks == 1 || n_associated_tracks == 3) baseline_taus.push_back(t);
      }

      // Overlap removal
      // Lepton/Jet
      Jets selected_small_R_jets;
      foreach (const Jet& j, baseline_small_R_jets) {
        bool selected = true;
        foreach (const Particle& e, baseline_electrons) {if (deltaR(j, e) < 0.2) selected = false;}
        foreach (const Particle& m, baseline_muons) {
          if (deltaR(j, m) < 0.4) {
            unsigned int n_associated_tracks = 0;
            foreach (const Particle& c, j.constituents()) {
              if (c.pT() > 0.5/GeV && c.threeCharge() != 0) n_associated_tracks++;
            }
            if (n_associated_tracks < 3) selected = false;
          }
        }
        if (selected) selected_small_R_jets.push_back(j);
      }

      // Jet/Lepton
      Particles selected_muons;
      foreach (const Particle& mu, baseline_muons) {
        bool selected = true;
        double cone_radius = 0.04 + (10.0/GeV)/mu.pT();
        if (cone_radius > 0.4) cone_radius = 0.4;
        foreach (const Jet& j, selected_small_R_jets) {
          if (deltaR(mu, j) < cone_radius) selected = false;
        }
        if (selected) selected_muons.push_back(mu);
      }

      Particles selected_electrons;
      foreach (const Particle& e, baseline_electrons) {
        bool selected = true;
        foreach (const Jet& j, selected_small_R_jets) {
          if (deltaR(e, j) < 0.4) selected = false;
        }
        if (selected) selected_electrons.push_back(e);
      }

      // Electron/Tau
      Particles selected_taus_pos, selected_taus_neg;
      foreach (const Particle& tau, baseline_taus) {
        bool selected = true;
        foreach (const Particle& e, selected_electrons) {
          if (deltaR(tau, e) < 0.1) selected = false;
        }
        if (selected) {
          if (tau.threeCharge() > 0) selected_taus_pos.push_back(tau);
          else selected_taus_neg.push_back(tau);
        }
      }

      // Now define signal objects
      Particles signal_electrons, other_baseline_leptons;
      foreach (const Particle& e, selected_electrons) {
        bool signal = true;
        if (e.pT() <  28.0/GeV) signal = false;
        double cone_radius = (10.0/GeV)/e.pT();
        if (cone_radius > 0.2) cone_radius = 0.2;
        double track_pT_sum = -e.pT();
        foreach (const Particle& t, tracks) {
          if (deltaR(t, e) < 0.2) track_pT_sum += t.pT();
        }
        if (track_pT_sum > 0.06*e.pT()) signal = false;
        if (signal) signal_electrons.push_back(e);
        else other_baseline_leptons.push_back(e);
      }

      Particles signal_muons;
      foreach (const Particle& mu, selected_muons) {
        bool signal = true;
        if (mu.pT() <  28.0/GeV) signal = false;
        double cone_radius = (10.0/GeV)/mu.pT();
        if (cone_radius > 0.3) cone_radius = 0.3;
        double track_pT_sum = -mu.pT();
        foreach (const Particle& t, tracks) {
          if (deltaR(t, mu) < 0.2) track_pT_sum += t.pT();
        }
        if (track_pT_sum > 0.06*mu.pT()) signal = false;
        if (signal) signal_muons.push_back(mu);
        else other_baseline_leptons.push_back(mu);
      }

      Jets signal_small_R_jets, signal_small_R_bjets;
      foreach (const Jet& j, selected_small_R_jets) {
        if (j.pT() < 25.0/GeV) continue;
        if (j.abseta() > 2.5) continue;
        signal_small_R_jets.push_back(j);
        if (j.bTagged()) signal_small_R_bjets.push_back(j);
      }

      // Construct large-R jets
      fastjet::JetDefinition large_R_jetdef(fastjet::antikt_algorithm, 1.0);
      fastjet::ClusterSequence large_R_cluster(signal_small_R_jets, large_R_jetdef);
      std::vector<fastjet::PseudoJet> baseline_large_R_jets = sorted_by_pt(large_R_cluster.inclusive_jets(200.0/GeV));
      std::vector<fastjet::PseudoJet> trimmed_large_R_jets;
      foreach (const Jet& j, baseline_large_R_jets) {
        trimmed_large_R_jets.push_back(_trimmer(j));
      }

      // Pre-selection
      if ((signal_electrons.size() + signal_muons.size()) != 1) vetoEvent;
      Particle L = signal_electrons.size() ? signal_electrons[0] : signal_muons[0];
      if (other_baseline_leptons.size() > 0) vetoEvent;
      if (signal_small_R_jets.size() < 4) vetoEvent;
      if (ETmiss < 300.0/GeV) vetoEvent;
      for (unsigned int i = 0; i < 2; i++) {
        if (deltaPhi(cal_invisible, signal_small_R_jets[i]) < 0.4) vetoEvent;
      }
      if (signal_small_R_bjets.size() < 1) vetoEvent;
      double mTW = sqrt(2*L.pT()*ETmiss*(1 - cos(deltaPhi(L, cal_invisible))));
      if (mTW < 30.0/GeV) vetoEvent;

      _c_passed_presel->fill(weight);

      if (ETmiss < 350.0/GeV) vetoEvent;
      if (mTW < 170.0/GeV) vetoEvent;
      // Insert mT2 cut if we really need to, although this relies on b-tagging weights etc.
      // Neglect HTmiss,sig cut since it relies on JES resolution
      if (signal_small_R_jets[0].pT() < 120.0/GeV) vetoEvent;
      if (signal_small_R_jets[1].pT() < 80.0/GeV) vetoEvent;
      if (signal_small_R_jets[2].pT() < 50.0/GeV) vetoEvent;
      double large_R_pT_cut = ETmiss > 450.0/GeV ? 200.0/GeV : 290.0/GeV;
      std::vector<fastjet::PseudoJet> signal_large_R_jets;
      foreach (const fastjet::PseudoJet& j, trimmed_large_R_jets) {
        if (j.perp() > large_R_pT_cut) signal_large_R_jets.push_back(j);
      }
      if (signal_large_R_jets.size() < 2) vetoEvent;
      if (signal_large_R_jets[0].m() < 80.0/GeV) vetoEvent;
      if (signal_large_R_jets[1].m() < 60.0/GeV) vetoEvent;

      _c_passed_SR->fill(weight);

    }

  void finalize() {

  }

private:

  const fastjet::Filter _trimmer;
  CounterPtr _c_passed_presel, _c_passed_SR;

};

// The hook for the plugin system
DECLARE_RIVET_PLUGIN(SVB_TpTp_ATLAS_R2);

}
