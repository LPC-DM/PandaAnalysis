/**
 * \file EnergyCorrelations.h
 * \brief Optimized code to calculate energy correlation functions. Based on code from fj-contrib
 * \author S.Narayanan
 */
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <map>
#include "TMath.h"
#include "TString.h"
#include "PandaCore/Tools/interface/Functions.h"

#ifndef PANDA_ECF_H
#define PANDA_ECF_H

namespace pandaecf {
  /**
   * \brief delta-r-squared metric between two pseudojets
   * @param  j1 first jet
   * @param  j2 second jet
   * @return    \f$dR^2\f$
   */
  double DeltaR2(fastjet::PseudoJet j1, fastjet::PseudoJet j2);

  /**
   * \brief Calculates un-normalized ECFs.
   *
   * If any of the pointers are not provided (or NULL), that value will not be calculated
   * @param beta         angular paramater
   * @param constituents particles with which to calculate the correlations
   * @param n1           pointer to location of N=1 ECF
   * @param n2           pointer to location of N=2 ECF
   * @param n3           pointer to location of N=3 ECF
   * @param n4           pointer to location of N=4 ECF
   */
  void calcECF(double beta, const std::vector<fastjet::PseudoJet> &constituents, 
               double *n1=0, double *n2=0, double *n3=0, double *n4=0);

  /**
   * \brief Just a bunch of floats and bools ot hold different values of normalized ECFs
   */
  class ECFNManager {
  public:
    ECFNManager() {
      flags["3_1"]=true; 
      flags["3_2"]=true; 
      flags["3_3"]=true; 
      flags["4_1"]=true; 
      flags["4_2"]=true; 
      flags["4_3"]=false; 
    }
    ~ECFNManager() {}

    std::map<TString,double> ecfns; //!< maps "N_I" to ECFN
    std::map<TString,bool>   flags; //!< maps "N_I" to flag

    bool doN1=true, doN2=true, doN3=true, doN4=true;

  };

  /**
   * \brief Calculates normalized energy correlation functions
   * @param beta         angular parameter
   * @param constituents particles with which to calculate the correlations
   * @param manager      provides configuration and storage of ECFNs
   * @param useMin       DEPRECATED
   */
  void calcECFN(double beta, const std::vector<fastjet::PseudoJet> &constituents, 
                ECFNManager *manager, bool useMin=true);
}

#endif
