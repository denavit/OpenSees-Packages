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
// $Date: 2010-09-17 23:49:28 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/shenSteelRCFT/shenSteelRCFT.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) shenSteelRCFT.C, revA"

#include <elementAPI.h>
#include "shenSteelRCFT.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <fstream>
using namespace std;

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

OPS_Export void
localInit() 
{
	  OPS_Error("shenSteelRCFT uniaxial material written by Mark D Denavit, University of Illinois at Urbana-Champaign \n", 1);
}

OPS_Export void *
OPS_shenSteelRCFT()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[19];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial shenSteelRCFT tag\n" << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
	opserr << "WARNING invalid input, want: uniaxialMaterial shenSteelRCFT tag Fy Fu Es R epo localBucklingStrain Ksft frs Rbso Epoi alfa a bb c w ksi e fE \n";
    return 0;	
  }

  // Check the data
  if ( dData[0] <= 0 ) {
	  opserr << "WARNING Fy should be a positive value\n";
	  return 0;
  }
  if ( dData[1] < dData[0] ) {
	  opserr << "WARNING Fy greater than Fu\n";
	  return 0;
  }  
  if ( dData[2] <= 0 ) {
	  opserr << "WARNING Es should be a positive value\n";
	  return 0;
  }  
  // @todo check the rest of the data
  
  
  theMaterial = new shenSteelRCFT(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		  dData[6], dData[7], dData[8], dData[9], dData[10],
		  dData[11], dData[12], dData[13], dData[14], dData[15],
		  dData[16], dData[17]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type shenSteelRCFT \n";
    return 0;
  }
 
  return theMaterial;
}




shenSteelRCFT::shenSteelRCFT(int tag, double i1, double i2, double i3, double i4, double i5,
		double i6, double i7, double i8, double i9, double i10, 
		double i11, double i12, double i13, double i14, double i15, 
		double i16, double i17, double i18)
:UniaxialMaterial(tag,MAT_TAG_shenSteelRCFT),
 fy(i1), fu(i2), Ee(i3), R(i4), epo(i5),
 localBucklingStrain(i6), Ksft(i7), frs(i8), Rbso(i9),
 Epoi(i10), alfa(i11), a(i12), bb(i13), c(i14),
 w(i15), ksi(i16), e(i17), fE(i18)
{
  this->revertToStart();
}

shenSteelRCFT::shenSteelRCFT()
:UniaxialMaterial(0,MAT_TAG_shenSteelRCFT),
 fy(0.0), fu(0.0), Ee(0.0), R(0.0), epo(0.0),
 localBucklingStrain(0.0), Ksft(0.0), frs(0.0), Rbso(0.0),
 Epoi(0.0), alfa(0.0), a(0.0), bb(0.0), c(0.0),
 w(0.0), ksi(0.0), e(0.0), fE(0.0)
{
  this->revertToStart();
}

shenSteelRCFT::~shenSteelRCFT()
{
  // does nothing
}

int 
shenSteelRCFT::setTrialStrain(double strain, double strainRate)
{
//    // Open file for debugging if needed
//    ofstream shenSteelRCFTdebug;
//    shenSteelRCFTdebug.open("shenSteelRCFTdebug.dat",ios::app);

	// Return state variables to their latest committed value
	this->revertToLastCommit();	

	// Set trial strain, strain increment, and loading direction
	trialStrain = strain;
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
		opserr << "Trouble in shenSteelRCFT::setTrialStrain (1) \n";
		return -1;
		break;
	}

	// Initialize Variables
	double delta, h, Ep, Epo;
	Epo = Epoi / ( 1 + w * W ); // Slope of the Bounding Surface
	
	// Check for unloading from tensile plasticity and update appropriate variables
	if ( commitedPlasticityStatus == 1 && plasticityStatus != 1 ){
		// Local Buckling Stuff
		if ( localBucklingHistory == 0 ) {
			// Reference strain for first local buckling if unloading from tensile branch
			localBucklingReferenceStrain = committedStrain - committedStress / Ee;
		} else if ( localBucklingHistory == 2 ) {
			// Set localBucklingConstantResidualStress
			// Reduction in constant residual stress based on plastic work
			localBucklingCyclicReduction = ( 1.0 - 5.0 * sqrt ( (W-localBucklingReferenceWork) / fy ) * R );
			if ( localBucklingCyclicReduction < 0.05 ) localBucklingCyclicReduction = 0.05;
			localBucklingConstantResidualStress = localBucklingCyclicReduction*localBucklingBaseConstantResidualStress;
			localBucklingBoundingStress = localBucklingConstantResidualStress;
		}
		
		// Virtual Bounding Line Stuff
		delta_yForC = (Epo * ep + Cmax_strs) - committedStress;
		if ( delta_yForC <= 0.0)
			delta_yForC = 0.0;
	}
	
	// Check for unloading from compressive plasticity and update appropriate variables
	if ( commitedPlasticityStatus == 2 && plasticityStatus != 2 ){	
		delta_yForT = committedStress - (Epo * ep - Cmax_strs);
		if ( delta_yForT <= 0.0)
			delta_yForT = 0.0;
	}

	// Set trial stress and tangent according to the plasticity model
	switch ( plasticityStatus ) {
	case 0 :
		// Elastic - stress point is within loading surfaces
		trialTangent = Ee;
		trialStress = committedStress + strainIncrement * trialTangent;
		break;
		
	case 1 :
		// Moving in the tensile or positive direction
		lastYieldedIn = 1;

		// Compute the stress and tangent
		if ( commitedPlasticityStatus != 1 ) {
			// If immediately after elastic loading then compute delta_in
			if ( lastYieldedIn != commitedLastYieldedIn ) {
				delta_in = fabs( Epo * ep + Tbs_p - Tls_p);
			}

			delta = delta_in;
			h = e * delta + fE * Ee;
			if ( committedStress < Epo * ep + Cmax_strs ) {
				// Within the memory surface
				Ep = Epo + h * ( delta + delta_yForT ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
			} else {
				// Outside the memory surface
				Ep = Epo + h * ( delta ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
			}
			if( localBucklingHistory != 0 ) {
				double reduction = ( 1.0 - 20 * sqrt ( W / fy ) * R );
				if ( reduction < 0.05 ) reduction = 0.05;
				Ep = Ep * reduction;
			}
			trialTangent = Ee * Ep / ( Ee + Ep );

			// Correct the overshooting of the yield surface
			trialStress = Tls_p + ( strainIncrement - (Tls_p-committedStress)/Ee ) * trialTangent;

			// Make sure that the loading point does not breach the bounding surface
			if ( trialStress > Epo * ep + Tbs_p ) {
				trialStress = Epo * ep + Tbs_p;
				trialTangent =  Epo;
			}
		} else {
			// If continued plastic loading
			delta = fabs( Epo * ep + Tbs_p - committedStress);
			h = e * delta + fE * Ee;
			if ( committedStress < Epo * ep + Cmax_strs ) {
				// Within the memory surface
				Ep = Epo + h * ( delta + delta_yForT ) / ( delta_in - delta );
			} else {
				// Outside the memory surface
				Ep = Epo + h * ( delta ) / ( delta_in - delta );
			}
			if( localBucklingHistory != 0 ) {
				double reduction = ( 1.0 - 20.0 * sqrt( W / fy ) * R );
				if ( reduction < 0.05 ) reduction = 0.05;
				Ep = Ep * reduction;
			}
			trialTangent = Ee * Ep / ( Ee + Ep );

			trialStress = committedStress + strainIncrement * trialTangent;

			// Make sure that the loading point does not breach the bounding surface
			if ( trialStress > Epo * ep + Tbs_p ) {
				trialStress = Epo * ep + Tbs_p;
				trialTangent =  Epo;
			}
		}
		break;
		
	case 2 :
		// Moving in the compressive or negative direction
		lastYieldedIn = 2;

		// Compute the stress and tangent
		if ( commitedPlasticityStatus != 2 ) {
			// If immediately after elastic loading then compute delta_in
			if ( lastYieldedIn != commitedLastYieldedIn ) {
				delta_in = fabs( Epo * ep + Tbs_n - Tls_n);
			}
			delta = delta_in;
			h = e * delta + fE * Ee;
			if ( committedStress > Epo * ep - Cmax_strs ) {
				// Within the memory surface
				Ep = Epo + h * ( delta + delta_yForC ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
			} else {
				// Outside the memory surface
				Ep = Epo + h * ( delta ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
			}
			trialTangent = Ee * Ep / ( Ee + Ep );
			// Correct the overshooting of the yield surface
			trialStress = Tls_n + ( strainIncrement - (Tls_n-committedStress) / Ee ) * trialTangent;

			// Make sure that the loading point does not breach the bounding surface
			if ( trialStress <= Epo * ep + Tbs_n ) {
				trialStress = Epo * ep + Tbs_n;
				trialTangent = Epo;
			}
		} else {
			// If continued plastic loading
			delta = fabs( Epo * ep + Tbs_n - committedStress);
			h = e * delta + fE * Ee;
			if ( committedStress > Epo * ep - Cmax_strs ) {
				// Within the memory surface
				Ep = Epo + h * ( delta + delta_yForC ) / ( delta_in - delta );
			} else {
				// Outside the memory surface
				Ep = Epo + h * ( delta ) / ( delta_in - delta );
			}
			trialTangent = Ee * Ep / ( Ee + Ep );
			trialStress = committedStress + strainIncrement * trialTangent;

			// Make sure that the loading point does not breach the bounding surface
			if ( trialStress <= Epo * ep + Tbs_n ) {
				trialStress = Epo * ep + Tbs_n;
				trialTangent = Epo;
			}
		}
		break;
		
	default :
		// It should not arrive here
		opserr << "Trouble in shenSteelRCFT::setTrialStrain (2) \n";
		return -1;
		break;
	}


	// Determine status w.r.t. the local buckling model
	if ( commitedLocalBucklingHistory == 0 ) {
		// Detect first local buckling using a strain measure
		if ( trialStrain <= (localBucklingStrain+localBucklingReferenceStrain) && trialStress <= -1*frs*max(localBucklingStrain*Ee,-fy)) {
			// Local buckling has initiated

			if ( committedStrain <= (localBucklingStrain+localBucklingReferenceStrain) ) {
				// if it could have buckled before but didn't because of the stress constraint
				localBucklingBaseConstantResidualStress = -1*frs*max(localBucklingStrain*Ee,-fy);
				localBucklingConstantResidualStress = localBucklingBaseConstantResidualStress;
				localBucklingStatus = 2;
			} else {
				// if it has achieved higher than the constant residual stress using fy

				// determine the new local buckling bounding stress
				double strainBeforeLocalBuckling = (localBucklingStrain+localBucklingReferenceStrain) - committedStrain;
				double stressAtLocalBuckling = committedStress + strainBeforeLocalBuckling*trialTangent;

				if ( stressAtLocalBuckling >= -1*frs*max(localBucklingStrain*Ee,-fy) ) {
					// the stress must have reached at least -1*frs*fy
					stressAtLocalBuckling = -1*frs*max(localBucklingStrain*Ee,-fy);
				}


				// compute the constant local buckling residual stress
				localBucklingBaseConstantResidualStress = frs*stressAtLocalBuckling;
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
			if ( trialStress < localBucklingBoundingStress ) {
				switch ( commitedLocalBucklingHistory ) {
				case 1:
					// strain at which constant residual stress begins
					{
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
					opserr << "Trouble in shenSteelRCFT::setTrialStrain (3) \n";
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
			opserr << "Trouble in shenSteelRCFT::setTrialStrain (4) \n";
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
		opserr << "Trouble in shenSteelRCFT::setTrialStrain (5) \n";
		return -1;
		break;
	}
	
	
	// Update state variables

	// Plastic Strain and Accumulated Plastic Work
	ep = trialStrain - trialStress/Ee;
	if( ep > epmax )
		epmax = ep;
	if( ep < epmin )
		epmin = ep;
	double ebar_p, dep;
	if( epmax - epmin > epo) {
		ebar_p = epmax - epmin;
	} else {
		ebar_p = epo;
	}
	dep = ep - Cep;
	W = W + fmax(trialStress*dep,0);
	if ( localBucklingHistory == 2 && commitedLocalBucklingHistory != 2)
		localBucklingReferenceWork = W;

	// Size of the Bounding Surface
	Tbs_p =   ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi  * pow ( 0.5 * ebar_p, 2 ) ) );
	Tbs_n = - ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi  * pow ( 0.5 * ebar_p, 2 ) ) );
	if ( localBucklingHistory != 0 && Tbs_n > localBucklingBoundingStress ) {
		Tbs_n = localBucklingBoundingStress;
	}

	// Loading Surface
	if ( plasticityStatus != 0 ) {
		double Rls;
		if ( strainIncrement <= 0.0 ) {
			// Moving in the compressive or negative direction
			Tls_n = trialStress;
			Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 ) - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
			if( localBucklingHistory != 0 ) {
				double reduction = ( 1.0 - 30.0 * sqrt( W / fy ) * R );
				if ( reduction < 0.05 ) reduction = 0.05;
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
	if ( Tls_p > Tbs_p )
		Tls_p = Tbs_p;
	if ( Tls_n < Tbs_n )
		Tls_n = Tbs_n;

	// Maximum Stress Seen by the Material
	if( fabs(trialStress) > Cmax_strs )
		Tmax_strs = fabs(trialStress);

//	if ( this->getTag() == 4 )
//		shenSteelRCFTdebug << "Tag #4 plasticityStatus " << plasticityStatus << "  localBucklingStatus " << localBucklingStatus << "\n";

	return 0;
}

double 
shenSteelRCFT::getStrain(void)
{
  return trialStrain;
}

double 
shenSteelRCFT::getStress(void)
{
  return trialStress;
}

double 
shenSteelRCFT::getInitialTangent(void)
{
	return Ee;
}

double 
shenSteelRCFT::getTangent(void)
{
  return trialTangent;
}

int 
shenSteelRCFT::commitState(void)
{
	committedTangent = trialTangent;
	committedStress = trialStress;
	committedStrain = trialStrain;

	commitedPlasticityStatus = plasticityStatus;
	commitedLocalBucklingStatus = localBucklingStatus;
	commitedLocalBucklingHistory = localBucklingHistory;
	committedLoadingDirection = loadingDirection;
	commitedLastYieldedIn = lastYieldedIn;

	Cls_p = Tls_p;	Cls_n = Tls_n;
	Cbs_p = Tbs_p;	Cbs_n = Tbs_n;
	Cmem_p = Tmem_p;  Cmem_n = Tmem_n;
	Cvbs_p = Tvbs_p;  Cvbs_n = Tvbs_n;
	Cep = ep; 
	Cepmin = epmin;  Cepmax = epmax;
	CW = W; 
	Cmax_strs = Tmax_strs; 
	
	Cdelta_in = delta_in;

	Cdelta_yForT = delta_yForT;
	Cdelta_yForC = delta_yForC;
	committedLocalBucklingReferenceStrain = localBucklingReferenceStrain;
	committedLocalBucklingCyclicReduction = localBucklingCyclicReduction;
	
	committedLocalBucklingConstantResidualStrain = localBucklingConstantResidualStrain; 
	committedLocalBucklingReferenceWork = localBucklingReferenceWork;
	committedLocalBucklingBaseConstantResidualStress = localBucklingBaseConstantResidualStress;
	committedLocalBucklingConstantResidualStress = localBucklingConstantResidualStress; 
	committedLocalBucklingBoundingStress = localBucklingBoundingStress;

	return 0;
}


int 
shenSteelRCFT::revertToLastCommit(void)
{
	trialTangent = committedTangent;
	trialStress = committedStress;
	trialStrain = committedStrain;

	plasticityStatus = commitedPlasticityStatus;
	localBucklingStatus = commitedLocalBucklingStatus;
	localBucklingHistory = commitedLocalBucklingHistory;
	loadingDirection = committedLoadingDirection;
	lastYieldedIn = commitedLastYieldedIn;

	Tls_p = Cls_p;  Tls_n = Cls_n;
	Tbs_p = Cbs_p;  Tbs_n = Cbs_n;
	Tmem_p = Cmem_p;  Tmem_n = Cmem_n;
	Tvbs_p = Cvbs_p;  Tvbs_n = Cvbs_n;
	ep = Cep; 
	epmin = Cepmin;  epmax = Cepmax;
	W = CW; 
	Tmax_strs = Cmax_strs; 
	
	delta_in = Cdelta_in;

	delta_yForT = Cdelta_yForT;
	delta_yForC = Cdelta_yForC;
	localBucklingReferenceStrain = committedLocalBucklingReferenceStrain;
	localBucklingCyclicReduction = committedLocalBucklingCyclicReduction;

	localBucklingConstantResidualStrain = committedLocalBucklingConstantResidualStrain; 
	localBucklingReferenceWork = committedLocalBucklingReferenceWork;
	localBucklingBaseConstantResidualStress = committedLocalBucklingBaseConstantResidualStress;
	localBucklingConstantResidualStress = committedLocalBucklingConstantResidualStress; 
	localBucklingBoundingStress = committedLocalBucklingBoundingStress; 
	
	return 0;
}


int 
shenSteelRCFT::revertToStart(void)
{
	Rlso = fy;

	trialTangent = Ee;
	trialStress = 0.0;
	trialStrain = 0.0;

	plasticityStatus = 0;
	localBucklingStatus = 0;
	localBucklingHistory = 0;
	loadingDirection = 0;
	lastYieldedIn = 0;

	double Rls;
	Rls = fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) );
	Tls_p =  Rls; Tls_n = -Rls;
	Tbs_p = Rbso; Tbs_n = -Rbso;
	Tmem_p = Tls_p; Tmem_n = Tls_n;
	Tvbs_p = 0.0; Tvbs_n = 0.0;
	ep = 0.0;
	epmin = 0.0; epmax = 0.0;
	W = 0.0;
	Tmax_strs = 0.0;

	delta_in = 0.0;

	delta_yForT = 0.0;
	delta_yForC = 0.0;
	localBucklingReferenceStrain = 0.0;
	localBucklingCyclicReduction = 1.0;

	localBucklingConstantResidualStrain = -fy/Ee;
	localBucklingReferenceWork = 0.0;
	localBucklingBaseConstantResidualStress = -fy;
	localBucklingConstantResidualStress = -fy;
	localBucklingBoundingStress = -fy;

	this->commitState();
    
	return 0;
}


UniaxialMaterial *
shenSteelRCFT::getCopy(void)
{
	// @todo work on this.
  
	
  shenSteelRCFT *theCopy = new shenSteelRCFT(this->getTag(),fy,fu,Ee,R,epo,localBucklingStrain,Ksft,
		  frs,Rbso,Epoi,alfa,a,bb,c,w,ksi,e,fE);
  //theCopy->ep = this->ep;
 
  return theCopy;
}


int 
shenSteelRCFT::sendSelf(int cTag, Channel &theChannel)
{
	// @todo work on this.

//  int res = 0;
//  static Vector data(6);
//  data(0) = this->getTag();
//  data(1) = ep;
//  data(2) = E;
//  data(3) = ezero;
//  data(4) = fyp;
//  data(5) = fyn;
//
//  res = theChannel.sendVector(this->getDbTag(), cTag, data);
//  if (res < 0) 
//    opserr << "shenSteelRCFT::sendSelf() - failed to send data\n";
//
//  return res;
	return 1;
}

int 
shenSteelRCFT::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
	// @todo work on this.

//  int res = 0;
//  static Vector data(6);
//  res = theChannel.recvVector(this->getDbTag(), cTag, data);
//  if (res < 0) 
//    opserr << "shenSteelRCFT::recvSelf() - failed to recv data\n";
//  else {
//    this->setTag(data(0));
//    ep    = data(1);
//    E     = data(2);
//    ezero = data(3);
//    fyp   = data(4);
//    fyn   = data(5);  
//  }
//
//  return res;
	return 1;
}

void 
shenSteelRCFT::Print(OPS_Stream &s, int flag)
{
  s << "shenSteelRCFT tag: " << this->getTag() << endln;
  s << " Fy:" << fy << endln;
  s << " Fu:" << fu << endln;
  s << " Es: " << Ee << endln;
  s << " Ksft:" << Ksft << endln;
  return;
}
