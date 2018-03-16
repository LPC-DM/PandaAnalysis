/**
 * \file EnergyCorrelations.cc
 * \brief Optimized code to calculate energy correlation functions. Based on code from fj-contrib
 * \author S.Narayanan
 */
#include "../interface/EnergyCorrelations.h"
#define PI 3.141592654

double pandaecf::DeltaR2(fastjet::PseudoJet j1, fastjet::PseudoJet j2) 
{
    float dEta2 = (j1.eta()-j2.eta()); dEta2 *= dEta2;

    float dPhi = j1.phi()-j2.phi();
    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;

    return dEta2 + dPhi*dPhi;
}

void pandaecf::calcECF(double beta, const std::vector<fastjet::PseudoJet> &constituents, 
                       double *n1/*=0*/, double *n2/*=0*/, double *n3/*=0*/, double *n4/*=0*/) 
{
  unsigned int nC = constituents.size();
  double halfBeta = beta/2.;

  // if only N=1,2, do not bother caching kinematics
  if (!n3 && !n4) {
    if (n1) { // N=1
      double val=0;
      for (unsigned int iC=0; iC!=nC; ++iC) {    
        val += constituents[iC].perp();
      }
      *n1 = val;
    }
    if (n2) { // N=2
      double val=0;
      for (unsigned int iC=0; iC!=nC; ++iC) {
        fastjet::PseudoJet iconst = constituents[iC];
        for (unsigned int jC=0; jC!=iC; ++jC) {
          fastjet::PseudoJet jconst = constituents[jC];
          val += iconst.perp() * jconst.perp() * pow(pandaecf::DeltaR2(iconst,jconst),halfBeta);
        }
      }
      *n2 = val;
    }
    return;
  }

  // cache kinematics
  std::vector<double> pTs(nC);
  std::vector<std::vector<double>> dRs(nC);
  for (unsigned int iC=0; iC!=nC; ++iC) {
    dRs[iC].resize(iC);
  }

  for (unsigned int iC=0; iC!=nC; ++iC) {
    fastjet::PseudoJet iconst = constituents[iC];
    pTs[iC] = iconst.perp();
    for (unsigned int jC=0; jC!=iC; ++jC) {
      fastjet::PseudoJet jconst = constituents[jC];
      dRs[iC][jC] = pow(pandaecf::DeltaR2(iconst,jconst),halfBeta);
    }
  }
  
  // now we calculate the real ECFs
  if (n1) { // N=1
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {    
      val += pTs[iC];
    } // iC
    *n1 = val;
  }
  if (n2) { // N=2
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        val += pTs[iC] * pTs[jC] * dRs[iC][jC]; 
      } // jC
    } // iC
    *n2 = val;
  }
  if (n3) {
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC]*dRs[iC][jC];
        for (unsigned int kC=0; kC!=jC; ++kC) {
          val += val_ij * pTs[kC] * dRs[iC][kC] * dRs[jC][kC];
        } // kC
      } // jC
    } // iC
    *n3 = val;
  }
  if (n4) {
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC]*dRs[iC][jC];
        for (unsigned int kC=0; kC!=jC; ++kC) {
          double val_ijk = val_ij * pTs[kC] * dRs[iC][kC] * dRs[jC][kC];
          for (unsigned int lC=0; lC!=kC; ++lC) {
            val += val_ijk * pTs[lC] * dRs[iC][lC] * dRs[jC][lC] * dRs[kC][lC];
          } // lC
        } // kC
      } // jC
    } // iC
    *n4 = val;
  }

}

/* TODO: switch from manual sort to std::partial_sort for readability*/
void pandaecf::calcECFN(double beta, const std::vector<fastjet::PseudoJet> &constituents, 
                        pandaecf::ECFNManager *manager, bool useMin/*=true*/) 
{
  unsigned int nC = constituents.size();
  double halfBeta = beta/2.;

  // get the normalization factor
  double baseNorm=0; 
  pandaecf::calcECF(beta,constituents,&baseNorm,0,0,0);

  // cache kinematics
  std::vector<double> pTs(nC);
  std::vector<std::vector<double>> dRs(nC);
  for (unsigned int iC=0; iC!=nC; ++iC) {
    dRs[iC].resize(iC);
  }

  for (unsigned int iC=0; iC!=nC; ++iC) {
    fastjet::PseudoJet iconst = constituents[iC];
    pTs[iC] = iconst.perp();
    for (unsigned int jC=0; jC!=iC; ++jC) {
      fastjet::PseudoJet jconst = constituents[jC];
      dRs[iC][jC] = pow(pandaecf::DeltaR2(iconst,jconst),halfBeta);
    }
  }
  
  // now we calculate the ECFNs
  if (manager->doN1) { // N=1
    manager->ecfns["1_1"] = 1;
    manager->ecfns["1_2"] = 1;
    manager->ecfns["1_3"] = 1;
  }
  if (manager->doN2) { // N=2
    double norm = pow(baseNorm,2);
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        val += pTs[iC] * pTs[jC] * dRs[iC][jC]; 
      } // jC
    } // iC
    val /= norm;
    manager->ecfns["2_1"] = val;
    manager->ecfns["2_2"] = val;
    manager->ecfns["2_3"] = val;
  }

  bool doI1=manager->flags["3_1"];
  bool doI2=manager->flags["3_2"];
  bool doI3=manager->flags["3_3"];
  if (manager->doN3 && (doI1||doI2||doI3)) {
    double norm = pow(baseNorm,3);
    double val1=0,val2=0,val3=0;
    unsigned int nAngles=3;
    std::vector<double> angles(nAngles);

    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC];
        angles[0] = dRs[iC][jC];

        for (unsigned int kC=0; kC!=jC; ++kC) {
          angles[1] = dRs[iC][kC];
          angles[2] = dRs[jC][kC];

          if (doI1||doI2||doI3) {
            double angle_1=999; unsigned int index_1=999;
            for (unsigned int iA=0; iA!=nAngles; ++iA) {
              if (angles[iA]<angle_1) {
                angle_1 = angles[iA];
                index_1 = iA;
              }
            }
            if (doI2||doI3) {
              double angle_2=999; 
              for (unsigned int jA=0; jA!=nAngles; ++jA) {
                if (jA==index_1) continue;
                if (angles[jA]<angle_2) {
                  angle_2 = angles[jA];
                }
              }
              if (doI3) {
                val3 += val_ij * pTs[kC] * angles[0] * angles[1] * angles[2];
              }
              if (doI2)
                val2 += val_ij * pTs[kC] * angle_1 * angle_2;
            }
            if (doI1)
              val1 += val_ij * pTs[kC] * angle_1;
          }

        } // kC
      } // jC
    } // iC
    val1 /= norm; val2 /= norm; val3 /= norm;
    manager->ecfns["3_1"] = val1;
    manager->ecfns["3_2"] = val2;
    manager->ecfns["3_3"] = val3;
  }

  doI1=manager->flags["4_1"];
  doI2=manager->flags["4_2"];
  if (manager->doN4 && (doI1||doI2)) {
    double norm = pow(baseNorm,4);
    double val1=0,val2=0;
    unsigned int nAngles=6;
    std::vector<double> angles(nAngles);

    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC];
        angles[0] = dRs[iC][jC];

        for (unsigned int kC=0; kC!=jC; ++kC) {
          double val_ijk = val_ij * pTs[kC];
          angles[1] = dRs[iC][kC];
          angles[2] = dRs[jC][kC];

          for (unsigned int lC=0; lC!=kC; ++lC) {
            angles[3] = dRs[iC][lC];
            angles[4] = dRs[jC][lC];
            angles[5] = dRs[kC][lC];

            if (doI1||doI2) {
              double angle_1=999; unsigned int index_1=999;
              for (unsigned int iA=0; iA!=nAngles; ++iA) {
                if (angles[iA]<angle_1) {
                  angle_1 = angles[iA];
                  index_1 = iA;
                }
              }
              if (doI2) {
                double angle_2=999; 
                for (unsigned int jA=0; jA!=nAngles; ++jA) {
                  if (jA==index_1) continue;
                  if (angles[jA]<angle_2) {
                    angle_2 = angles[jA];
                  }
                }
                val2 += val_ijk * pTs[lC] * angle_1 * angle_2;
              }
              if (doI1)
                val1 += val_ijk * pTs[lC] * angle_1;
            }
          } // lC
        } // kC
      } // jC
    } // iC
    val1 /= norm; val2 /= norm;
    manager->ecfns["4_1"] = val1;
    manager->ecfns["4_2"] = val2;
    manager->ecfns["4_3"] = 0;
  }

}
