/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2010/05/04 17:14:46 $
// $Source: /cvsroot/openseescomp/CompositePackages/shenSteel01/shenSteel01.cpp,v $

// Documentation: Shen Steel Model
// uniaxialMaterial shenSteel01 $tag $Es $Fy $Fu $eu $kappaBar0 $Ep0i $alpha $a $b $c $omega $zeta $e $f <options>
//
// Required Input Parameters:
//   $tag                   integer tag identifying material
//   $Es                    modulus of elasticity
//   $Fy                    yield stress
//   $Fu                    ultimate stress
//   $eu                    strain at ultimate stress
//   $kappaBar0             size of the initial bounding lines
//   $Ep0i                  slope of the initial bounding lines
//   $alpha, $a, $b, $c     constants in the equation to compute the reduction in the elastic range
//   $omega                 constant in the equation to compute the slope of the bounding line
//   $zeta                  constant in the equation to compute the size of the bounding line
//   $e, $f                 constants in the equation to compute the shape parameter (for the plastic modulus)
//
// Optional Input:
//   -hotRolled $epst $M $Epst
//     $epst                plastic strain at the end of the yield plateau under monotonic loading
//     $M                   constant in the equation to determine if the yield plateau still continues
//     $Epst                plastic modulus at initial hardening under monotonic loading
//   -coldFormed $ep0 $Epst
//     $ep0                 value of the initial plastic strain
//     $Epst                plastic modulus at initial hardening under monotonic loading
//   -coldFormed2 $ep0
//     $ep0                 value of the initial plastic strain
//   -initialStress $sigma0
//     $sigma0              value of the initial stress (must be in the elastic range)
//   -biaxialStress $alphaLat
//     $alphaLat            ratio of stress in the lateral direction to yield stress (must be in the elastic range)
//   -localBuckling $elb $Ksft $alphaFulb $ref
//     $elb                 strain at initiation of local buckling
//     $Ksft                slope of the descending branch of the local buckling model
//     $alphaFulb           ratio of ultimate local buckling stress to value defined in the next variable
//     $ref                 reference stress for alphaFulb (value in the denominator, either Fy or Flb)
//   -localBucklingDegradationEp $rate $limit
//     $rate                rate of degradation of Ep (plastic modulus) after local buckling
//     $limit               limit of degradation of Ep (plastic modulus) after local buckling
//   -localBucklingDegradationKappa $rate $limit
//     $rate                rate of degradation of kappa (size of the elastic range) after local buckling
//     $limit               limit of degradation of kappa (size of the elastic range) after local buckling
//   -localBucklingDegradationFulb $rate $limit
//     $rate                rate of degradation of Fulb (ultimate local buckling stress) after local buckling
//     $limit               limit of degradation of Fulb (ultimate local buckling stress) after local buckling
//
// References:
//   1. C. Shen et al., "Cyclic Behavior of Structural Steels. II: Theory," Journal of Engineering
//      Mechanics 121, no. 11 (1995): 1165-1172.
//   2. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.

#include <elementAPI.h>
#include "shenSteel01.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

// This is only needed for writing things to file for debugging
#include <fstream>
using namespace std;

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

OPS_Export void localInit() {
  OPS_Error("shenSteel01 unaxial material \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign \n", 1);
}

OPS_Export void *OPS_shenSteel01() {
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // Variables to retrieve input
  double dData[20];
  int sDataLength = 40;
  char *sData = new char[sDataLength];
  int numData;

  // Variables to temporarily store data
  // Required (no need to initialize)
  int tag;
  double Es, Fy, Fu, eu, kappaBar0, Ep0i, alpha, a, b, c, omega, zeta, e, f;
  // Option (initialize to default values)
  bool modelYieldPlateau = false;
  double epst = 0.0;
  double Epst = 0.0;
  double M = 0.0;
  double sigma0 = 0.0;
  double alphaLat = 0.0;
  double ep0 = 0.0;
  bool modelLocalBuckling = false;
  bool modelFirstExcursion = false;
  double elb = 0.0;
  double Ksft = 0.0;
  double alphaFulb = 0.0;
  int refFulb = 0;
  bool modelDegradeEp = false;
  double degradeEpRate = 0.0;
  double degradeEpLimit = 0.0;
  bool modelDegradeKappa = false;
  double degradeKappaRate = 0.0;
  double degradeKappaLimit = 0.0;
  bool modelDegradeFulb = false;
  double degradeFulbRate = 0.0;
  double degradeFulbLimit = 0.0;

  // Parse the input
  numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial shenSteel01 tag\n" << endln;
    return 0;
  }

  numData = 14;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input, want: uniaxialMaterial shenSteel01 tag Es Fy Fu eu kappaBar0 Ep0i alpha a b c omega zeta e f <options> \n";
    return 0;
  }

  Es = dData[0];  Fy = dData[1];  Fu = dData[2];  eu = dData[3];
  kappaBar0 = dData[4];  Ep0i = dData[5];
  alpha = dData[6];  a = dData[7];  b = dData[8];  c= dData[9];
  omega = dData[10];  zeta = dData[11];
  e = dData[12]; f = dData[13];

  // Check the data
  if ( Fy <= 0 ) {
    opserr << "WARNING Fy should be a positive value\n";
    return 0;
  }
  if ( Fu < Fy ) {
    opserr << "WARNING Fy greater than Fu\n";
    return 0;
  }
  if ( Es <= 0 ) {
    opserr << "WARNING Es should be a positive value\n";
    return 0;
  }
  if ( eu < Fy/Es ) {
    opserr << "WARNING eu less than ey = Fy/Es \n";
    return 0;
  }

  // Loop through remaining arguments
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
    if ( OPS_GetStringCopy(&sData) != 0 ) {
      opserr << "WARNING invalid input";
      return 0;
    }

    if ( strcmp(sData,"-hotRolled") == 0 ) {
      numData = 3;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -hotRolled $epst $M $Epst \n";
        return 0;
      }
      epst = dData[0];  M = dData[1];  Epst = dData[2];
      modelYieldPlateau = true;
      modelFirstExcursion = true;

    } else if ( strcmp(sData,"-coldFormed") == 0 ) {
      numData = 2;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -coldFormed $ep0 $Epst \n";
        return 0;
      }
      ep0 = dData[0]; Epst = dData[1];
      modelYieldPlateau = false;
      modelFirstExcursion = true;

    } else if ( strcmp(sData,"-coldFormed2") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -coldFormed2 $ep0 \n";
        return 0;
      }
      ep0 = dData[0];
      modelYieldPlateau = false;
      modelFirstExcursion = false;

    } else if ( strcmp(sData,"-initialStress") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -initialStress $sigma0 \n";
        return 0;
      }
      sigma0 = dData[0];
      if ( fabs(sigma0) >= Fy ) {
        opserr << "WARNING invalid input, initial stress (sigma0) is greater than Fy \n";
        return 0;
      }

    } else if ( strcmp(sData,"-biaxialStress") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -biaxialStress $alphaLat \n";
        return 0;
      }
      alphaLat = dData[0];

    } else if ( strcmp(sData,"-localBuckling") == 0 ) {
      numData = 3;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -localBuckling $elb $Ksft $alphaFulb $ref \n";
        return 0;
      }
      elb = dData[0];  Ksft = dData[1];  alphaFulb = dData[2];
      if ( OPS_GetStringCopy(&sData) != 0 ) {
        opserr << "WARNING invalid input";
        return 0;
      }
      if ( strcmp(sData,"Fy") == 0 ) {
        refFulb = 1;
      } else if ( strcmp(sData,"Flb") == 0 ) {
        refFulb = 2;
      } else {
        opserr << "WARNING invalid $ref, options are Fy and Flb, you put " << sData << "\n";
        return 0;
      }
      modelLocalBuckling = true;

    } else if ( strcmp(sData,"-localBucklingDegradationEp") == 0 ) {
      if ( modelLocalBuckling == false ) {
        opserr << "WARNING invalid input, set -localBucklingDegradationEp after -localBuckling \n";
        return 0;
      }
      numData = 2;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -localBucklingDegradationEp $rate $limit \n";
        return 0;
      }
      degradeEpRate = dData[0];  degradeEpLimit = dData[1];
      modelDegradeEp = true;

    } else if ( strcmp(sData,"-localBucklingDegradationKappa") == 0 ) {
      if ( modelLocalBuckling == false ) {
        opserr << "WARNING invalid input, set -localBucklingDegradationKappa after -localBuckling \n";
        return 0;
      }
      numData = 2;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -localBucklingDegradationKappa $rate $limit \n";
        return 0;
      }
      degradeKappaRate = dData[0];  degradeKappaLimit = dData[1];
      modelDegradeKappa = true;

    } else if ( strcmp(sData,"-localBucklingDegradationFulb") == 0 ) {
      if ( modelLocalBuckling == false ) {
        opserr << "WARNING invalid input, set -localBucklingDegradationFulb after -localBuckling \n";
        return 0;
      }
      numData = 2;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -localBucklingDegradationFulb $rate $limit \n";
        return 0;
      }
      degradeFulbRate = dData[0];  degradeFulbLimit = dData[1];
      modelDegradeFulb = true;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  theMaterial = new shenSteel01(tag, Fy, Fu, Es, eu,
      kappaBar0, Ep0i, alpha, a, b, c, omega, zeta, e, f,
      modelYieldPlateau, M, Epst, epst,
      modelFirstExcursion, sigma0, alphaLat, ep0,
      modelLocalBuckling, elb, Ksft, alphaFulb, refFulb,
      modelDegradeEp, degradeEpRate, degradeEpLimit,
      modelDegradeKappa, degradeKappaRate, degradeKappaLimit,
      modelDegradeFulb, degradeFulbRate, degradeFulbLimit);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type shenSteel01 \n";
    return 0;
  }

  return theMaterial;
}


shenSteel01::shenSteel01(int tag, double iA1, double iA2, double iA3, double iA4,
    double iB1, double iB2, double iB3, double iB4, double iB5, double iB6, double iB7, double iB8, double iB9, double iB10,
    bool iC1, double iC2, double iC3, double iC4,
    bool iD1, double iD2, double iD3, double iD4,
    bool iE1, double iE2, double iE3, double iE4, int iE5,
    bool iF1, double iF2, double iF3,
    bool iG1, double iG2, double iG3,
    bool iH1, double iH2, double iH3)
:UniaxialMaterial(tag,MAT_TAG_shenSteel01),
 fy(iA1), fu(iA2), Ee(iA3), eu(iA4),
 Rbso(iB1), Epoi(iB2), alfa(iB3), a(iB4), bb(iB5), c(iB6), w(iB7), ksi(iB8), e(iB9), fE(iB10),
 modelYieldPlateau(iC1), M(iC2), Est(iC3), est(iC4),
 modelFirstExcursion(iD1), initStress(iD2), alphaLat(iD3), ep0(iD4),
 modelLocalBuckling(iE1), localBucklingStrain(iE2), Ksft(iE3), alphaFulb(iE4), refFulb(iE5),
 modelDegradeEp(iF1), degradeEpRate(iF2), degradeEpLimit(iF3),
 modelDegradeKappa(iG1), degradeKappaRate(iG2), degradeKappaLimit(iG3),
 modelDegradeFulb(iH1), degradeFulbRate(iH2), degradeFulbLimit(iH3)
{
  Rlso = fy * ( alfa - a * exp( -bb * ep0 * 100 ) - ( alfa - a - 1 ) * exp( -c * ep0 * 100 ) );
  this->revertToStart();
  if ( initStress <= Tls_n || initStress >= Tls_p ) {
    opserr << "WARNING initial stress is not in the elastic range, setting initial stress to zero \n";
    initStress = 0.0;
  }
}

shenSteel01::shenSteel01(int tag)
  :UniaxialMaterial(tag,MAT_TAG_shenSteel01)
{
  // does nothing
}

shenSteel01::~shenSteel01() {
  // does nothing
}

int shenSteel01::setTrialStrain(double strain, double strainRate) {

  // Return state variables to their latest committed value
  this->revertToLastCommit();

  // Set trial strain, strain increment, and loading direction
  if ( initStress == 0) {
    trialStrain = strain;
  } else {
    trialStrain = strain + initStress/Ee;
  }
  double strainIncrement;
  strainIncrement = trialStrain - committedStrain;
  if ( strainIncrement > 0.0 ) {
    loadingDirection = 1;
  } else {
    loadingDirection = -1;
  }

  // Check if the strain has changed
  if ( trialStrain == committedStrain )
    return 0;

  // Determine status w.r.t. plasticity model (i.e., ignoring local buckling)
  switch ( commitedPlasticityStatus ) {
  case 0 :
    // Elastic
    double stressIfElastic;
    stressIfElastic = committedStress + Ee * strainIncrement;
    if ( stressIfElastic > Tls_p ) {
      // Above the loading surface --> tensile plasticity
      plasticityStatus = 1;
    } else if ( stressIfElastic < Tls_n ) {
      // Below the loading surface --> compressive plasticity
      plasticityStatus = 2;
    } else {

    }
    break;

  case 1 :
    // Tensile Plasticity
    if ( loadingDirection == 1 ) {
      // Continued tensile plasticity
      plasticityStatus = 1;
    } else {
      double stressIfElastic;
      stressIfElastic = committedStress + Ee * strainIncrement;
      if ( stressIfElastic >= Tls_n ) {
        // Within the loading surface --> elastic
        plasticityStatus = 0;
      } else {
        // Outside the loading surface --> compressive plasticity
        plasticityStatus = 2;
      }
    }
    break;

  case 2 :
    // Compressive Plasticity
    if ( loadingDirection == -1 ) {
      // Continued compressive plasticity
      plasticityStatus = 2;
    } else {
      double stressIfElastic;
      stressIfElastic = committedStress + Ee * strainIncrement;
      if ( stressIfElastic <= Tls_p ) {
        // Within the loading surface --> elastic
        plasticityStatus = 0;
      } else {
        // Outside the loading surface --> tensile plasticity
        plasticityStatus = 1;
      }
    }
    break;

  default :
    // It should not arrive here
    opserr << "Trouble in shenSteel01::setTrialStrain (1) \n";
    return -1;
    break;
  }

  // Check for unloading from tensile plasticity and update appropriate variables
  if ( commitedPlasticityStatus == 1 && plasticityStatus != 1 ){
    // Local Buckling Stuff
    if ( localBucklingHistory == 0 ) {
      // Reference strain for first local buckling if unloading from tensile branch
      localBucklingReferenceStrain = committedStrain - committedStress / Ee;
    } else if ( modelDegradeFulb == true && localBucklingHistory == 2 ) {
      // Set localBucklingConstantResidualStress
      // Reduction in constant residual stress based on plastic work
      localBucklingCyclicReduction = ( 1.0 - degradeFulbRate * sqrt ( (W-localBucklingReferenceWork) / fy ) );
      if ( localBucklingCyclicReduction < degradeFulbLimit ) {
        localBucklingCyclicReduction = degradeFulbLimit;
      }
      localBucklingConstantResidualStress = localBucklingCyclicReduction*localBucklingBaseConstantResidualStress;
      localBucklingBoundingStress = localBucklingConstantResidualStress;
    }

    // Virtual Bounding Line
    delta_yForC = (Epo * ep + Cmax_strs) - committedStress;
    if ( delta_yForC <= 0.0) {
      delta_yForC = 0.0;
    }
  }

  // Check for unloading from compressive plasticity and update appropriate variables
  if ( commitedPlasticityStatus == 2 && plasticityStatus != 2 ){
    // Virtual Bounding Line
    delta_yForT = committedStress - (Epo * ep - Cmax_strs);
    if ( delta_yForT <= 0.0) {
      delta_yForT = 0.0;
    }
  }


  // Set trial stress and tangent according to the plasticity model
  Epo = Epoi / ( 1 + w * W );

  switch ( plasticityStatus ) {
  case 0 :
    // Elastic - stress point is within loading surfaces
    trialTangent = Ee;
    trialStress = committedStress + strainIncrement * trialTangent;
    break;

  case 1 :
    // Moving in the tensile or positive direction
    if ( firstExcursionStatus == 1 && lastYieldedIn == 2 ) {
      firstExcursionStatus = 0;
    }
    lastYieldedIn = 1;

    // Compute the stress and tangent based on the plasticity formulation
    if ( commitedPlasticityStatus != 1 ) {
      // If immediately after elastic loading then compute delta_in
      if ( lastYieldedIn != commitedLastYieldedIn ) {
        delta_in = fabs( Epo * ep + Tbs_p - Tls_p);
      }
      double plasticStrainIncrement;
      plasticStrainIncrement = strainIncrement - (Tls_p-committedStress)/Ee;
      plasticityModel(plasticStrainIncrement, Tls_p, Ee, trialStress, trialTangent);
    } else {
      plasticityModel(strainIncrement, committedStress, committedTangent, trialStress, trialTangent);
    }

    // Yield Plateau
    if ( modelYieldPlateau == true && commitedYieldPlateauStatus != 0 && trialStress >= fy ) {
      if ( commitedPlasticityStatus == 0 ) {
        delta_in = fabs( Epo * ep + Tbs_p - fy);
      }

      // Assess if the yield plateau will still be active at the end of the step
      double epIfPlateau, epMaxIfPlateau, epBarIfPlateau, WIfPlateau, plateauTest;
      epIfPlateau = trialStrain - fy/Ee;
      if( epIfPlateau > epmax ) {
        epMaxIfPlateau = epIfPlateau;
      } else {
        epMaxIfPlateau = epmax;
      }
      epBarIfPlateau = epMaxIfPlateau - epmin;
      WIfPlateau = W + (epIfPlateau-Cep)*fy;
      plateauTest = ( epBarIfPlateau/est -1 ) - M * ( WIfPlateau/est/fy -1);


      if ( plateauTest < 0 ) {
        // Yield plateau continues
        yieldPlateauStatus = 1;
        trialStress = fy;
        trialTangent = 0.0;
      } else {
        // If yield plateau vanishes, then set appropriate values
        yieldPlateauStatus = 0;
        double strainIncrementInPlaeau;
        strainIncrementInPlaeau = 0.5*strainIncrement;
        esh = trialStrain-strainIncrementInPlaeau;
        trialStress = fy + (strainIncrement-strainIncrementInPlaeau)*Est;
        trialTangent = Est;
        firstExcursionStress = fy;
      }
    }

    // First Excursion
    if ( modelFirstExcursion == true && commitedFirstExcursionStatus != 0 && trialStress > fy && trialStress >= firstExcursionStress ) {
      if ( esh == 0.0 ) {
        esh = committedStrain;
      }
      // Stress-Strain Relationship for First Excursion
      if (trialStrain < eu) {
        double p = Est*((eu-esh)/(fu-fy));
        if (p > 2.0) {
          trialStress  = fu + (fy-fu)*pow((eu-trialStrain)/(eu-esh),p);
          trialTangent = Est*pow((eu-trialStrain)/(eu-esh),p-1);
        } else {
          double fsh2 = 0.5*(fy+fu);
          double esh2 = esh + (fsh2-fy)/Est;
          if (trialStrain < esh2) {
            trialStress  = fy + Est*(trialStrain-esh);
            trialTangent = Est;
          } else {
            double Est2 = (fu-fsh2)/(eu-esh2);
            trialStress  = fsh2 + Est2*(trialStrain-esh2);
            trialTangent = Est2;
          }
        }
      } else {
        trialStress = fu;
        trialTangent = 0.0;
      }

      // Update First Excursion Stress
      firstExcursionStress = trialStress;
    }

    break;

  case 2 :
    // Moving in the compressive or negative direction
    if ( firstExcursionStatus == 1 && lastYieldedIn == 1 ) {
      firstExcursionStatus = 0;
    }
    lastYieldedIn = 2;


    // Compute the stress and tangent based on the plasticity formulation
    if ( commitedPlasticityStatus != 2 ) {
      // If immediately after elastic loading then compute delta_in
      if ( lastYieldedIn != commitedLastYieldedIn ) {
        delta_in = fabs( Epo * ep + Tbs_n - Tls_n);
      }
      double plasticStrainIncrement;
      plasticStrainIncrement = strainIncrement - (Tls_n-committedStress)/Ee;
      plasticityModel(plasticStrainIncrement, Tls_n, Ee, trialStress, trialTangent);
    } else {
      plasticityModel(strainIncrement, committedStress, committedTangent, trialStress, trialTangent);
    }


    // Yield Plateau
    if ( modelYieldPlateau == true && commitedYieldPlateauStatus != 0 && trialStress <= -fy ) {
      if ( commitedPlasticityStatus == 0 ) {
        delta_in = fabs( Epo * ep + Tbs_n - (-fy));
      }

      // Assess if the yield plateau will still be active at the end of the step
      double epIfPlateau, epMinIfPlateau, epBarIfPlateau, WIfPlateau, plateauTest;
      epIfPlateau = trialStrain + fy/Ee;
      if( epIfPlateau > epmin ) {
        epMinIfPlateau = epIfPlateau;
      } else {
        epMinIfPlateau = epmin;
      }
      epBarIfPlateau = epmax - epMinIfPlateau;
      WIfPlateau = W + (epIfPlateau-Cep)*(-fy);
      plateauTest = ( epBarIfPlateau/est -1 ) - M * ( WIfPlateau/est/fy -1);


      if ( plateauTest < 0 ) {
        // Yield plateau continues
        yieldPlateauStatus = 1;
        trialStress = -fy;
        trialTangent = 0;
      } else {
        // If yield plateau vanishes, then set appropriate values
        yieldPlateauStatus = 0;
        double strainIncrementInPlaeau;
        strainIncrementInPlaeau = 0.5*strainIncrement; //est - (Cep + M*CW/fy)/(1-M);
        esh = trialStrain-strainIncrementInPlaeau;
        trialStress = -fy + (strainIncrement-strainIncrementInPlaeau)*Est;
        trialTangent = Est;
        firstExcursionStress = -fy;
      }
    }

    // First Excursion
    if ( modelFirstExcursion == true && commitedFirstExcursionStatus != 0 && trialStress < -fy && trialStress <= firstExcursionStress ) {
      if ( esh == 0.0 ) {
        esh = committedStrain;
      }

      // Stress-Strain Relationship for First Excursion
      if (trialStrain > -eu) {
        double p = Est*((eu+esh)/(fu-fy));
        if (p > 2.0) {
          trialStress  = -fu - (fy-fu)*pow((-eu-trialStrain)/(-eu-esh),p);
          trialTangent = Est*pow((-eu-trialStrain)/(-eu-esh),p-1);
        } else {
          double fsh2 = -0.5*(fy+fu);
          double esh2 = esh + (fsh2+fy)/Est;
          if (trialStrain > esh2) {
            trialStress  = -fy + Est*(trialStrain-esh);
            trialTangent = Est;
          } else {
            double Est2 = (-fu-fsh2)/(-eu-esh2);
            trialStress  = fsh2 + Est2*(trialStrain-esh2);
            trialTangent = Est2;
          }
        }
      } else {
        trialStress = -fu;
        trialTangent = 0.0;
      }

      // Update First Excursion Stress
      firstExcursionStress = trialStress;

    }

    break;

  default :
    // It should not arrive here
    opserr << "Trouble in shenSteel01::setTrialStrain (2) \n";
    return -1;
    break;
  }


  // Determine status w.r.t. the local buckling model
  if ( modelLocalBuckling ) {
    if ( commitedLocalBucklingHistory == 0 ) {
      // Detect first local buckling using a strain measure (but stress must also be lower than a certain amount)
      double stressForLocalBuckling;
      if ( refFulb == 1 ) {
        stressForLocalBuckling = alphaFulb*(-fy);
      } else if ( refFulb == 2 ) {
        double scratch;
        scratch = localBucklingStrain*Ee;
        if ( scratch < -fy )
          scratch = -fy;
        stressForLocalBuckling = alphaFulb*scratch;
      } else {
        opserr << "shenSteel01::setTrialStrain unknown refFulb: "<<refFulb<<"\n";
        return -1;
      }

      if ( trialStrain <= (localBucklingStrain+localBucklingReferenceStrain) && trialStress <= stressForLocalBuckling ) {
        // Local buckling has initiated

        if ( committedStrain <= (localBucklingStrain+localBucklingReferenceStrain) ) {
          // if it could have buckled before but didn't because of the stress constraint
          localBucklingBaseConstantResidualStress = stressForLocalBuckling;
          localBucklingConstantResidualStress = localBucklingBaseConstantResidualStress;
          localBucklingStatus = 2;
        } else {
          // if it has achieved higher than the constant residual stress using fy

          // determine the new local buckling bounding stress
          double strainBeforeLocalBuckling = (localBucklingStrain+localBucklingReferenceStrain) - committedStrain;
          double stressAtLocalBuckling = committedStress + strainBeforeLocalBuckling*trialTangent;

          if ( stressAtLocalBuckling >= stressForLocalBuckling ) {
            // the stress must have reached at least -1*frs*fy
            stressAtLocalBuckling = stressForLocalBuckling;
          }

          // compute the constant local buckling residual stress
          if ( refFulb == 1 ) {
            localBucklingBaseConstantResidualStress = alphaFulb*(-fy);
          } else if ( refFulb == 2 ) {
            localBucklingBaseConstantResidualStress = alphaFulb*stressAtLocalBuckling;
          } else {
            opserr << "shenSteel01::setTrialStrain unknown refFulb: "<<refFulb<<"\n";
            return -1;
          }
          localBucklingConstantResidualStress = localBucklingBaseConstantResidualStress;

          // strain at which constant residual stress begins
          localBucklingConstantResidualStrain = (localBucklingStrain+localBucklingReferenceStrain) - (stressAtLocalBuckling-localBucklingConstantResidualStress)/Ksft;

          // determine what branch of local buckling the material is in
          if ( trialStrain > localBucklingConstantResidualStrain ) {
            localBucklingStatus = 1;
          } else {
            localBucklingStatus = 2;
          }
        }
      } else {
        localBucklingStatus = 0;
      }
    } else {
      switch ( commitedLocalBucklingStatus ) {
      case 0 :
        // Detect subsequent local buckling using a stress measure

        double localBucklingStressLimit;
        double localBucklingBoundingStrain;
        localBucklingBoundingStrain = localBucklingConstantResidualStrain + (localBucklingBoundingStress-localBucklingConstantResidualStress)/Ksft;
        if ( trialStrain >= localBucklingBoundingStrain ) {
          localBucklingStressLimit = localBucklingBoundingStress;
        } else if ( trialStrain >= localBucklingConstantResidualStrain ) {
          localBucklingStressLimit = localBucklingConstantResidualStress - Ksft*(localBucklingConstantResidualStrain-trialStrain);
        } else {
          localBucklingStressLimit = localBucklingConstantResidualStress;
        }

        if ( trialStress < localBucklingStressLimit ) {
          switch ( commitedLocalBucklingHistory ) {
          case 1:
            // strain at which constant residual stress begins
            if ( trialStrain >= localBucklingBoundingStrain ) {
            double strainAtStartOfLB = committedStrain + (localBucklingBoundingStress-committedStress)/trialTangent;
            localBucklingConstantResidualStrain = strainAtStartOfLB - (localBucklingBoundingStress-localBucklingConstantResidualStress)/Ksft;
            }

            // determine what branch of local buckling the material is in
            if ( trialStrain > localBucklingConstantResidualStrain ) {
              localBucklingStatus = 1;
            } else {
              localBucklingStatus = 2;
            }
            break;

          case 2:
            localBucklingStatus = 2;
            break;

          default :
            // It should not arrive here
            opserr << "Trouble in shenSteel01::setTrialStrain (3) \n";
            return -1;
            break;
          }
        } else {
          localBucklingStatus = 0;
        }
        break;

      case 1 :
        if ( strainIncrement > 0.0 ) {
          localBucklingStatus = 0;
        } else {
          // determine what branch of local buckling the material is in
          if ( trialStrain > localBucklingConstantResidualStrain ) {
            localBucklingStatus = 1;
          } else {
            localBucklingStatus = 2;
          }
        }
        break;

      case 2 :
        if ( strainIncrement > 0.0 ) {
          localBucklingStatus = 0;
        } else {
          localBucklingStatus = 2;
        }
        break;

      default :
        // It should not arrive here
        opserr << "Trouble in shenSteel01::setTrialStrain (4) \n";
        return -1;
        break;
      }
    }

    // Set trial stress and tangent if it is locally buckled
    switch ( localBucklingStatus ) {
    case 0 :
      break;

    case 1 :
      trialStress = localBucklingConstantResidualStress + Ksft*(trialStrain - localBucklingConstantResidualStrain );
      trialTangent = Ksft;
      localBucklingHistory = 1;
      localBucklingBoundingStress = trialStress;
      break;

    case 2 :
      trialStress = localBucklingConstantResidualStress;
      trialTangent = 0.0;
      localBucklingHistory = 2;
      localBucklingBoundingStress = trialStress;
      break;

    default :
      // It should not arrive here
      opserr << "Trouble in shenSteel01::setTrialStrain (5) \n";
      return -1;
      break;
    }
  }

  // Update state variables

  // Plastic Strain and Accumulated Plastic Work
  ep = trialStrain - trialStress/Ee;
  if( ep > epmax )
    epmax = ep;
  if( ep < epmin )
    epmin = ep;
  double ebar_p, dep;
  if( epmax - epmin > ep0) {
    ebar_p = epmax - epmin;
  } else {
    ebar_p = ep0;
  }
  dep = ep - Cep;
  if (trialStress*dep > 0) {
    W = W + trialStress*dep;
  }

  if ( modelLocalBuckling && localBucklingHistory == 2 && commitedLocalBucklingHistory != 2)
    localBucklingReferenceWork = W;

  // Size of the Bounding Surface
  double rho = 0.5*ebar_p;
  Tbs_p =   ( fu + ( Rbso - fu ) * exp ( -ksi*rho*rho ) );
  Tbs_n = - ( fu + ( Rbso - fu ) * exp ( -ksi*rho*rho ) );
  if ( modelLocalBuckling && localBucklingHistory != 0 && Tbs_n < localBucklingBoundingStress ) {
    Tbs_n = localBucklingBoundingStress;
  }

  // Loading Surface
  if ( plasticityStatus != 0 ) {
    double Rls;
    if ( strainIncrement <= 0.0 ) {
      // Moving in the compressive or negative direction
      Tls_n = trialStress;
      Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 ) - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
      if( modelLocalBuckling && modelDegradeKappa == true && localBucklingHistory != 0 ) {
        double reduction = ( 1.0 - degradeKappaRate * sqrt( W / fy ) ); // @todo - maybe have (W-refW)
        if ( reduction < degradeKappaLimit ) reduction = degradeKappaLimit;
        Rls = Rls * reduction;
      }
      Tls_p = trialStress + 2 * Rls;
    } else {
      // Moving in the tensile or positive direction
      Tls_p = trialStress;
      Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 ) - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
      Tls_n = trialStress - 2 * Rls;
    }
  }

  // Make sure the loading surface is not bigger than the bounding surface
  if ( Tls_p > Epo * ep + Tbs_p )
    Tls_p = Epo * ep + Tbs_p;
  if ( Tls_n < Epo * ep + Tbs_n )
    Tls_n = Epo * ep + Tbs_n;

  // Maximum Stress Seen by the Material
  if( fabs(trialStress-Epo*ep) > Cmax_strs )
    Tmax_strs = fabs(trialStress-Epo*ep);

  return 0;
}


void shenSteel01::plasticityModel(double strainIncrement, double initialStress, double initialTangent, double &finalStress, double &finalTangent)
{
  // Initialize Variables
  double delta, h, Ep;

  if ( strainIncrement == 0.0 ) {
    // No Strain Increment
    finalStress = initialStress;
    finalTangent = initialTangent;
    return;

  } else if ( strainIncrement < 0.0) {
    // Compressive Strain Increment

    delta = fabs( Epo * ep + Tbs_n - initialStress);
    h = e * delta + fE;
    if ( initialStress > Epo * ep - Cmax_strs ) {
      // Within the memory surface
      if ( delta_in-delta <= DBL_EPSILON ) {
        Ep = Epo + h * ( delta + delta_yForC ) / 1e-15;
      } else {
        Ep = Epo + h * ( delta + delta_yForC ) / ( delta_in - delta );
      }
    } else {
      // Outside the memory surface
      if ( delta_in-delta <= DBL_EPSILON ) {
        Ep = Epo + h * ( delta ) / 1e-15;
      } else {
        Ep = Epo + h * ( delta ) / ( delta_in - delta );
      }
    }


  } else {
    // Tensile Strain Increment

    delta = fabs( Epo * ep + Tbs_p - initialStress);
    h = e * delta + fE;
    if ( initialStress < Epo * ep + Cmax_strs ) {
      // Within the memory surface
      if ( delta_in-delta <= DBL_EPSILON ) {
        Ep = Epo + h * ( delta + delta_yForT ) / 1e-15;
      } else {
        Ep = Epo + h * ( delta + delta_yForT ) / ( delta_in - delta );
      }
    } else {
      // Outside the memory surface
      if ( delta_in-delta <= DBL_EPSILON ) {
        Ep = Epo + h * ( delta ) / 1e-15;
      } else {
        Ep = Epo + h * ( delta ) / ( delta_in - delta );
      }
    }

    // Local Buckling Effects
    if( modelDegradeEp == true && localBucklingHistory != 0 ) {
      double reduction = ( 1.0 - degradeEpRate * sqrt ( W / fy ) ); // @todo - maybe have (W-refW)
      if ( reduction < degradeEpLimit ) reduction = degradeEpLimit;
      Ep = Ep * reduction;
    }
  }

  if (Ep < 0)
    Ep = 0;

  finalTangent = Ee * Ep / ( Ee + Ep );
  finalStress = initialStress + strainIncrement * finalTangent;

  // Make sure that the loading point does not breach the bounding surface
  if ( finalStress >= Epo * ep + Tbs_p ) {
    finalStress = Epo * ep + Tbs_p;
    finalTangent =  Epo;
  } else if ( finalStress <= Epo * ep + Tbs_n ) {
    finalStress = Epo * ep + Tbs_n;
    finalTangent = Epo;
  }


  return;
}

double shenSteel01::getStrain(void) {
  if (initStress == 0.0) {
    return trialStrain;
  } else {
    return trialStrain - initStress/Ee;
  }
}

double shenSteel01::getStress(void) {
  return trialStress;
}

double shenSteel01::getInitialTangent(void) {
  return Ee;
}

double shenSteel01::getTangent(void) {
  return trialTangent;
}

int shenSteel01::commitState(void) {
  // State Variables
  committedTangent = trialTangent;
  committedStress = trialStress;
  committedStrain = trialStrain;

  // Flags to determine the general state of the material
  commitedPlasticityStatus = plasticityStatus;
  commitedYieldPlateauStatus = yieldPlateauStatus;
  commitedFirstExcursionStatus = firstExcursionStatus;
  commitedLocalBucklingStatus = localBucklingStatus;
  commitedLocalBucklingHistory = localBucklingHistory;
  committedLoadingDirection = loadingDirection;
  commitedLastYieldedIn = lastYieldedIn;

  // Updated at "end" of each step
  Cls_p = Tls_p;
  Cls_n = Tls_n;
  Cbs_p = Tbs_p;
  Cbs_n = Tbs_n;
  Cmem_p = Tmem_p;
  Cmem_n = Tmem_n;
  Cvbs_p = Tvbs_p;
  Cvbs_n = Tvbs_n;
  Cep = ep;
  Cepmin = epmin;
  Cepmax = epmax;
  CW = W;
  Cmax_strs = Tmax_strs;

  // Updated upon transition from elastic to plastic
  Cdelta_in = delta_in;

  // Updated upon unloading from plasticity
  Cdelta_yForT = delta_yForT;
  Cdelta_yForC = delta_yForC;
  committedLocalBucklingReferenceStrain = localBucklingReferenceStrain;
  committedLocalBucklingCyclicReduction = localBucklingCyclicReduction;

  // local buckling state variables
  committedLocalBucklingConstantResidualStrain = localBucklingConstantResidualStrain;
  committedLocalBucklingReferenceWork = localBucklingReferenceWork;
  committedLocalBucklingBaseConstantResidualStress = localBucklingBaseConstantResidualStress;
  committedLocalBucklingConstantResidualStress = localBucklingConstantResidualStress;
  committedLocalBucklingBoundingStress = localBucklingBoundingStress;

  Cesh = esh;
  commitedFirstExcursionStress = firstExcursionStress;

  return 0;
}


int shenSteel01::revertToLastCommit(void) {
  // State Variables
  trialTangent = committedTangent;
  trialStress = committedStress;
  trialStrain = committedStrain;

  // Flags to determine the general state of the material
  plasticityStatus = commitedPlasticityStatus;
  yieldPlateauStatus = commitedYieldPlateauStatus;
  firstExcursionStatus = commitedFirstExcursionStatus;
  localBucklingStatus = commitedLocalBucklingStatus;
  localBucklingHistory = commitedLocalBucklingHistory;
  loadingDirection = committedLoadingDirection;
  lastYieldedIn = commitedLastYieldedIn;

  // Updated at "end" of each step
  Tls_p = Cls_p;
  Tls_n = Cls_n;
  Tbs_p = Cbs_p;
  Tbs_n = Cbs_n;
  Tmem_p = Cmem_p;
  Tmem_n = Cmem_n;
  Tvbs_p = Cvbs_p;
  Tvbs_n = Cvbs_n;
  ep = Cep;
  epmin = Cepmin;
  epmax = Cepmax;
  W = CW;
  Tmax_strs = Cmax_strs;

  // Updated upon transition from elastic to plastic
  delta_in = Cdelta_in;

  // Updated upon unloading from plasticity
  delta_yForT = Cdelta_yForT;
  delta_yForC = Cdelta_yForC;
  localBucklingReferenceStrain = committedLocalBucklingReferenceStrain;
  localBucklingCyclicReduction = committedLocalBucklingCyclicReduction;

  // local buckling state variables
  localBucklingConstantResidualStrain = committedLocalBucklingConstantResidualStrain;
  localBucklingReferenceWork = committedLocalBucklingReferenceWork;
  localBucklingBaseConstantResidualStress = committedLocalBucklingBaseConstantResidualStress;
  localBucklingConstantResidualStress = committedLocalBucklingConstantResidualStress;
  localBucklingBoundingStress = committedLocalBucklingBoundingStress;

  esh = Cesh;
  firstExcursionStress = commitedFirstExcursionStress;

  // Open file for debugging if needed
//  ofstream shenSteel01debug;
//  shenSteel01debug.open("shenSteel01debug.dat2",ios::app);
//  shenSteel01debug << Tls_p << " " << Tls_n << " " << Tbs_p << " " << Tbs_n << " " << W << " " << Epoi / ( 1 + w * W ) << "\n";
//  shenSteel01debug.close();

  return 0;
}


int shenSteel01::revertToStart(void) {
  // State Variables
  trialTangent = Ee;
  trialStress = initStress;
  trialStrain = initStress/Ee;

  // Flags to determine the general state of the material
  plasticityStatus = 0;
  if ( est <= 0.0) {
    yieldPlateauStatus = 0;
  } else {
    yieldPlateauStatus = 1;
  }
  firstExcursionStatus = 1;
  localBucklingStatus = 0;
  localBucklingHistory = 0;
  loadingDirection = 0;
  lastYieldedIn = 0;

  // Updated at "end" of each step
  double alphaZp = 0.5*(alphaLat+pow(4-3*alphaLat*alphaLat,0.5));
  double alphaZn = 0.5*(alphaLat-pow(4-3*alphaLat*alphaLat,0.5));
  Tls_p = alphaZp * Rlso;
  Tls_n = alphaZn * Rlso;
  Tbs_p = Rbso;
  Tbs_n = -Rbso;
  Tmem_p = Tls_p;
  Tmem_n = Tls_n;
  Tvbs_p = 0.0;
  Tvbs_n = 0.0;
  ep = 0.0;
  epmin = 0.0;
  epmax = 0.0;
  W = 0.0;
  Tmax_strs = fabs(initStress);

  // Updated upon transition from elastic to plastic
  delta_in = 0.0;

  // Updated upon unloading from plasticity
  delta_yForT = 0.0;
  delta_yForC = 0.0;
  localBucklingReferenceStrain = 0.0;
  localBucklingCyclicReduction = 1.0;

  // local buckling state variables
  localBucklingConstantResidualStrain = -fy/Ee;
  localBucklingConstantResidualStress = -fy;
  localBucklingBaseConstantResidualStress = -fy;
  localBucklingReferenceWork = 0.0;
  localBucklingBoundingStress = -fy;

  esh = 0.0;
  firstExcursionStress = 0.0;

  this->commitState();

  return 0;
}


UniaxialMaterial * shenSteel01::getCopy(void) {
  shenSteel01 *theCopy = new shenSteel01(this->getTag());

  // Input Variables
  theCopy-> fy = fy;
  theCopy-> fu = fu;
  theCopy-> Ee = Ee;
  theCopy-> eu = eu;
  theCopy-> Rbso = Rbso;
  theCopy-> Epoi = Epoi;
  theCopy-> alfa = alfa;
  theCopy-> a = a;
  theCopy-> bb = bb;
  theCopy-> c = c;
  theCopy-> w = w;
  theCopy-> ksi = ksi;
  theCopy-> e = e;
  theCopy-> fE = fE;
  theCopy-> modelYieldPlateau = modelYieldPlateau;
  theCopy-> M = M;
  theCopy-> Est = Est;
  theCopy-> est = est;
  theCopy-> modelFirstExcursion = modelFirstExcursion;
  theCopy-> initStress = initStress;
  theCopy-> alphaLat = alphaLat;
  theCopy-> ep0 = ep0;
  theCopy-> modelLocalBuckling = modelLocalBuckling;
  theCopy-> localBucklingStrain = localBucklingStrain;
  theCopy-> Ksft = Ksft;
  theCopy-> alphaFulb = alphaFulb;
  theCopy-> refFulb = refFulb;
  theCopy-> modelDegradeEp = modelDegradeEp;
  theCopy-> degradeEpRate = degradeEpRate;
  theCopy-> degradeEpLimit = degradeEpLimit;
  theCopy-> modelDegradeKappa = modelDegradeKappa;
  theCopy-> degradeKappaRate = degradeKappaRate;
  theCopy-> degradeKappaLimit = degradeKappaLimit;
  theCopy-> modelDegradeFulb = modelDegradeFulb;
  theCopy-> degradeFulbRate = degradeFulbRate;
  theCopy-> degradeFulbLimit = degradeFulbLimit;

  // Computed Material Properties
  theCopy-> Rlso = Rlso;

  // State Variables
  theCopy-> trialTangent = trialTangent;
  theCopy-> committedTangent = committedTangent;
  theCopy-> trialStress = trialStress;
  theCopy-> committedStress = committedStress;
  theCopy-> trialStrain = trialStrain;
  theCopy-> committedStrain = committedStrain;

  // Flags to determine the general state of the material
  theCopy-> plasticityStatus = plasticityStatus;
  theCopy-> commitedPlasticityStatus = commitedPlasticityStatus;
  theCopy-> yieldPlateauStatus = yieldPlateauStatus;
  theCopy-> commitedYieldPlateauStatus = commitedYieldPlateauStatus;
  theCopy-> firstExcursionStatus = firstExcursionStatus;
  theCopy-> commitedFirstExcursionStatus = commitedFirstExcursionStatus;
  theCopy-> localBucklingStatus = localBucklingStatus;
  theCopy-> commitedLocalBucklingStatus = commitedLocalBucklingStatus;
  theCopy-> localBucklingHistory = localBucklingHistory;
  theCopy-> commitedLocalBucklingHistory = commitedLocalBucklingHistory;
  theCopy-> loadingDirection = loadingDirection;
  theCopy-> committedLoadingDirection = committedLoadingDirection;
  theCopy-> lastYieldedIn = lastYieldedIn;
  theCopy-> commitedLastYieldedIn = commitedLastYieldedIn;

  // Updated at "end" of each step
  theCopy-> Tls_p = Tls_p;
  theCopy-> Cls_p = Cls_p;
  theCopy-> Tls_n = Tls_n;
  theCopy-> Cls_n = Cls_n;
  theCopy-> Tbs_p = Tbs_p;
  theCopy-> Cbs_p = Cbs_p;
  theCopy-> Tbs_n = Tbs_n;
  theCopy-> Cbs_n = Cbs_n;
  theCopy-> Tmem_p = Tmem_p;
  theCopy-> Cmem_p = Cmem_p;
  theCopy-> Tmem_n = Tmem_n;
  theCopy-> Cmem_n = Cmem_n;
  theCopy-> Tvbs_p = Tvbs_p;
  theCopy-> Cvbs_p = Cvbs_p;
  theCopy-> Tvbs_n = Tvbs_n;
  theCopy-> Cvbs_n = Cvbs_n;
  theCopy-> ep = ep;
  theCopy-> Cep = Cep;
  theCopy-> epmin = epmin;
  theCopy-> Cepmin = Cepmin;
  theCopy-> epmax = epmax;
  theCopy-> Cepmax = Cepmax;
  theCopy-> W = W;
  theCopy-> CW = CW;
  theCopy-> Cmax_strs = Cmax_strs;
  theCopy-> Tmax_strs = Tmax_strs;

  // Updated upon transition from elastic to plastic
  theCopy-> delta_in = delta_in;
  theCopy-> Cdelta_in = Cdelta_in;

  // Updated upon unloading from plasticity
  theCopy-> delta_yForT = delta_yForT;
  theCopy-> Cdelta_yForT = Cdelta_yForT;
  theCopy-> delta_yForC = delta_yForC;
  theCopy-> Cdelta_yForC = Cdelta_yForC;
  theCopy-> localBucklingReferenceStrain = localBucklingReferenceStrain;
  theCopy-> committedLocalBucklingReferenceStrain = committedLocalBucklingReferenceStrain;
  theCopy-> localBucklingCyclicReduction = localBucklingCyclicReduction;
  theCopy-> committedLocalBucklingCyclicReduction = committedLocalBucklingCyclicReduction;

  // local buckling state variables
  theCopy-> localBucklingConstantResidualStrain = localBucklingConstantResidualStrain;
  theCopy-> committedLocalBucklingConstantResidualStrain = committedLocalBucklingConstantResidualStrain;
  theCopy-> localBucklingReferenceWork = localBucklingReferenceWork;
  theCopy-> committedLocalBucklingReferenceWork = committedLocalBucklingReferenceWork;
  theCopy-> localBucklingBaseConstantResidualStress = localBucklingBaseConstantResidualStress;
  theCopy-> committedLocalBucklingBaseConstantResidualStress = committedLocalBucklingBaseConstantResidualStress;
  theCopy-> localBucklingConstantResidualStress = localBucklingConstantResidualStress;
  theCopy-> committedLocalBucklingConstantResidualStress = committedLocalBucklingConstantResidualStress;
  theCopy-> localBucklingBoundingStress = localBucklingBoundingStress;
  theCopy-> committedLocalBucklingBoundingStress = committedLocalBucklingBoundingStress;


  theCopy-> esh = esh;
  theCopy-> Cesh = Cesh;
  theCopy-> firstExcursionStress = firstExcursionStress;
  theCopy-> commitedFirstExcursionStress = commitedFirstExcursionStress;

  return theCopy;
}


int shenSteel01::sendSelf(int cTag, Channel &theChannel) {
  // @todo work on this.
  return -1;
}

int shenSteel01::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker) {
  // @todo work on this.
  return -1;
}

void shenSteel01::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s<<"shenSteel01, tag: "<<this->getTag()<<endln;
    s<<" Es:      "<<Ee<<endln;
    s<<" Fy:      "<<fy<<endln;
    s<<" Fu:      "<<fu<<endln;
    s<<" eu:      "<<eu<<endln;
    s<<" kappaBar0: "<<Rbso<<" Ep0i: "<<Epoi<<" alpha: "<<alfa<<" a: "<<a<<" b: "<<bb<<" c: "<<c<<" omega: "<<w<<" zeta: "<<ksi<<" e: "<<e<<" f: "<<fE<< endln;
    s<<" initStress: "<<initStress<<" alphaLat: "<<alphaLat<<" ep0: "<<ep0<<" Epst: "<<Est<< endln;
    if (modelYieldPlateau == true)
      s<<" modelYieldPlateau: true   M: "<<M<<" epst: "<<est<< endln;
    if (modelLocalBuckling == true)
      s<<" modelLocalBuckling: true  localBucklingStrain: "<<localBucklingStrain<<" Ksft: "<<Ksft<<" alphaFulb: "<<alphaFulb<<" refFulb: "<<refFulb<< endln;
    if (modelDegradeEp == true)
      s<<" modelDegradeEp: true      degradeEpRate: "<<degradeEpRate<<" degradeEpLimit: "<<degradeEpLimit<< endln;
    if (modelDegradeKappa == true)
      s<<" modelDegradeKappa: true   degradeKappaRate: "<<degradeKappaRate<<" degradeKappaLimit: "<<degradeKappaLimit<< endln;
    if (modelDegradeFulb == true)
      s<<" modelDegradeFulb: true    degradeFulbRate: "<<degradeFulbRate<<" degradeFulbLimit: "<<degradeFulbLimit<< endln;
  } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"shenSteel01\", ";
    s << "\"Es\": " << Ee << ", ";
    s << "\"Fy\": " << fy << ", ";
    s << "\"Fu\": " << fu << ", ";
    s << "\"eu\": " << eu << ", ";
    s << "\"kappaBar0\": " << Rbso << ", ";
    s << "\"Ep0i\": " << Epoi << ", ";
    s << "\"alpha\": " << alfa << ", ";
    s << "\"a\": " << a << ", ";
    s << "\"b\": " << bb << ", ";
    s << "\"c\": " << c << ", ";
    s << "\"omega\": " << w << ", ";
    s << "\"zeta\": " << ksi << ", ";
    s << "\"e\": " << e << ", ";
    s << "\"f\": " << fE << ", ";
    s << "\"initStress\": " << initStress << ", ";
    s << "\"alphaLat\": " << alphaLat << ", ";
    s << "\"ep0\": " << ep0 << ", ";
    s << "\"Epst\": " << Est;
    if (modelYieldPlateau == true)
      s << ", \"modelYieldPlateau\": true, \"M\": " << M << ", \"epst\": " << est;
    if (modelLocalBuckling == true)
      s << ", \"modelLocalBuckling\": true, \"localBucklingStrain\": " << localBucklingStrain << ", \"Ksft\": " << Ksft << ", \"alphaFulb\": " << alphaFulb << ", \"refFulb\": " << refFulb;
    if (modelDegradeEp == true)
      s << ", \"modelDegradeEp\": true, \"degradeEpRate\": " << degradeEpRate << ", \"degradeEpLimit\": " << degradeEpLimit;
    if (modelDegradeKappa == true)
      s << ", \"modelDegradeKappa\": true, \"degradeKappaRate\": " << degradeKappaRate << ", \"degradeKappaLimit\": " << degradeKappaLimit;
    if (modelDegradeFulb == true)
      s << ", \"modelDegradeFulb\": true, \"degradeFulbRate\": " << degradeFulbRate << ", \"degradeFulbLimit\": " << degradeFulbLimit;
    s << "}";
  }
  return;
}
