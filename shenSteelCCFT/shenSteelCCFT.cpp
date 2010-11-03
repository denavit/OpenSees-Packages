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
// $Date: 2010-05-04 17:14:46 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/shenSteelCCFT/shenSteelCCFT.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) shenSteelCCFT.C, revA"

#include <elementAPI.h>
#include "shenSteelCCFT.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


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
	  OPS_Error("shenSteelCCFT unaxial material \nWritten by Mark D Denavit, University of Illinois at Urbana-Champaign, Copyright 2010 \nUse at your Own Peril\n", 1);
}

OPS_Export void *
OPS_shenSteelCCFT()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[19];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial shenSteelCCFT tag\n" << endln;
    return 0;
  }

  numData = 19;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
	opserr << "WARNING invalid input, want: uniaxialMaterial shenSteelCCFT tag Fy Fu Es R epo alphaHoop localBucklingStrain Ksft Rcrit Rbso Epoi alfa a bb c w ksi e fE \n";
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
  
  
  theMaterial = new shenSteelCCFT(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		  dData[6], dData[7], dData[8], dData[9], dData[10],
		  dData[11], dData[12], dData[13], dData[14], dData[15],
		  dData[16], dData[17], dData[18]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type shenSteelCCFT \n";
    return 0;
  }
 
  return theMaterial;
}




shenSteelCCFT::shenSteelCCFT(int tag, double i1, double i2, double i3, double i4, double i5, 
		double i6, double i7, double i8, double i9, double i10, 
		double i11, double i12, double i13, double i14, double i15, 
		double i16, double i17, double i18, double i19)
:UniaxialMaterial(tag,MAT_TAG_shenSteelCCFT),
 fy(i1), fu(i2), Ee(i3), R(i4), epo(i5),
 alphaHoop(i6), localBucklingStrain(i7), Ksft(i8), Rcrit(i9), Rbso(i10),
 Epoi(i11), alfa(i12), a(i13), bb(i14), c(i15),
 w(i16), ksi(i17), e(i18), fE(i19)
{
  Rlso = fy;
  this->revertToStart();
}

shenSteelCCFT::shenSteelCCFT()
:UniaxialMaterial(0,MAT_TAG_shenSteelCCFT),
 fy(0.0), fu(0.0), Ee(0.0), R(0.0), epo(0.0),
 alphaHoop(0.0), localBucklingStrain(0.0), Ksft(0.0), Rcrit(0.0), Rbso(0.0),
 Epoi(0.0), alfa(0.0), a(0.0), bb(0.0), c(0.0),
 w(0.0), ksi(0.0), e(0.0), fE(0.0)
{
  Rlso = 0.0;
  this->revertToStart();
}

shenSteelCCFT::~shenSteelCCFT()
{
  // does nothing
}

int 
shenSteelCCFT::setTrialStrain(double strain, double strainRate)
{
	// Return state variables to their latest commited value
	backtocommitStatevar();	

	double strainIncrement, stressIncrement, stressIfElastic, dep;

	// Trial strain
	trialStrain = strain;
	strainIncrement = trialStrain - committedStrain;
	
	// Determine loading direction
	if ( strainIncrement >= 0.0 ) {
		loadingDirection = 1;
	} else {
		loadingDirection = 2;
	}
	
	// Determine if a reversal has occured (for computation of virtual bounding line)
	if ( loadingDirection != committedLoadingDirection && committedLoadingDirection != 0 && lastYieldedIn != commitedLastYieldedIn ) {
		// Reversal has occured
		if ( loadingDirection == 1 ) {
			delta_y = committedStress - (Epo * ep - Cmax_strs);			
		} else {
			delta_y = (Epo * ep + Cmax_strs) - committedStress;
		}
		if ( delta_y <= 0.0)
			delta_y = 0.0;

		//opserr<<"Reversal has occured, delta_y = "<<delta_y<<" loadingDirection = "<<loadingDirection<<" committedLoadingDirection = "<<committedLoadingDirection<<"\n";
	}
	
	
	// Determine whether elastic or plastc, is the stress point inside the loading surface
	// Calculate trial stress
	stressIfElastic = committedStress + committedTangent * strainIncrement;
	if( ( stressIfElastic < Tls_p ) && ( stressIfElastic > Tls_n ) ){
		// ELASTIC LOADING
		elastic = 1;
		plastic = 0;
		Trule = 1;
	} else {
		// PLASTIC LOADING 
		elastic = 0;
		plastic = 1;
		Trule = 2;
	}
	
	// Detect first local buckling using a strain measure (skip if local buckling has occured)
	if ( elb == 0 && trialStrain - localBucklingReferenceStrain <= localBucklingStrain ) {
		elastic = 0;
		plastic = 0;
		lb = 1;
		
		// stress at which local buckling occurred which is new local buckling bounding stress
		localBucklingBoundingStress = committedStress;
		
		// make sure the buckling stress is not positive
		if ( localBucklingBoundingStress > 0.0 )
			localBucklingBoundingStress = 0.0;
		
		// compute the constant local buckling residual stress
		if ( R > Rcrit ) {
			localBucklingConstantResidualStress = (Rcrit/R)*localBucklingBoundingStress;
		} else {
			localBucklingConstantResidualStress = localBucklingBoundingStress;
		}
		
		// strain at which constant residual stress begins 
		localBucklingConstantResidualStrain = committedStrain - (localBucklingBoundingStress-localBucklingConstantResidualStress)/Ksft;	
	}
	
	// If locally buckled: 
	//   1) compressive loading is locally buckled 
	//   2) tensile loading is elastic
	if ( Clb == 1 ){
		if ( strainIncrement <= 0.0 ) {
			elastic = 0;
			plastic = 0;
			lb = 1;
		} else {
			elastic = 1;
			plastic = 0;
			lb = 0;
		}
	}
	
	// Determine state if material is elastic
	if ( elastic == 1 ) {
		// ELASTIC LOADING ( Stress point is within loading surfaces )
		
		// New tangent and stress
		trialTangent = Ee;
		stressIncrement = strainIncrement * trialTangent;
		trialStress = committedStress + stressIncrement;
		
		// Update reference strain for first local buckling if unloading from tensile branch
		if ( strainIncrement < 0.0 && Crule == 2 ){
			localBucklingReferenceStrain = committedStrain - committedStress / Ee;
		}
		
		// Check if stress point is outside of the bounding surface
		if ( trialStress <=   Epo * ep + Tbs_n ){
			trialStress =   Epo * ep + Tbs_n;
			//opserr<<"Breach Bounding Surface\n";
		} 
		if ( trialStress >=   Epo * ep + Tbs_p ){
			trialStress =  Epo * ep + Tbs_p;
			//opserr<<"Breach Bounding Surface\n";
		}
		
	// Determine state if material is plastic	
	} else if ( plastic == 1 ) {
		// Plastic Loading ( Stress point is outside of loading surfaces )
		
		// Moving in the compressive or negative direction
		if ( strainIncrement <= 0.0 ) {
			lastYieldedIn = 2;
				
			// Compute the stress and tangent
			if ( Crule == 1) {
				// If immediatley after elastic loading then compute delta_in
				if ( lastYieldedIn != commitedLastYieldedIn ) {
					delta_in = fabs( Epo * ep + Tbs_n - Tls_n);
					//opserr<<" delta_in updated";
				}
				delta = delta_in;
				h = e * delta + fE * Ee;
				if ( committedStress > Epo * ep - Cmax_strs ) { 
					// Within the memory surface
					Ep = Epo + h * ( delta + delta_y ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
				} else {
					// Outside the memory surface
					Ep = Epo + h * ( delta ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
				}
				trialTangent = Ee * Ep / ( Ee + Ep );
				// Correct the overshooting of the yield surface
				trialStress = Tls_n + ( ( stressIfElastic - Tls_n ) / Ee ) * trialTangent;
				// Make sure that the loading point does not breach the bounding surface
				if ( trialStress <= Epo * ep + Tbs_n )
					trialStress = Epo * ep + Tbs_n;
				dep = ( ( stressIfElastic - Tls_n ) / Ee ) * trialTangent / Ep;
			} else {
				// If continuued plastic loading
				delta = fabs( Epo * ep + Tbs_n - committedStress);
				h = e * delta + fE * Ee;
				if ( committedStress > Epo * ep - Cmax_strs ) { 
					// Within the memory surface
					Ep = Epo + h * ( delta + delta_y ) / ( delta_in - delta );
				} else {
					// Outside the memory surface
					Ep = Epo + h * ( delta ) / ( delta_in - delta );				
				}
				trialTangent = Ee * Ep / ( Ee + Ep );
				trialStress = committedStress + strainIncrement * trialTangent;
				// Make sure that the loading point does not breach the bounding surface
				if ( trialStress <= Epo * ep + Tbs_n ) {
					trialStress = Epo * ep + Tbs_n;
					//opserr<<"Breach Bounding Surface\n";
				}
				//opserr<<" cont plastic, Epo = "<<Epo<<" h = "<<h<<" delta "<<delta<<" delta_y "<<delta_y<<" delta_in "<<delta_in<<"\n";
				//opserr<<" cont plastic, delta = "<<delta<<" h = "<<h<<" Ep "<<Ep<<" trialTangent "<<trialTangent<<" trialStress "<<trialStress<<"\n";
			}
			
		// Moving in the tensile or positive direction	
		} else {
			lastYieldedIn = 1;
			
			// Compute the stress and tangent
			if ( Crule == 1) {
				// If immediatley after elastic loading then compute delta_in
				if ( lastYieldedIn != commitedLastYieldedIn ) {
					delta_in = fabs( Epo * ep + Tbs_p - Tls_p);
					//opserr<<" delta_in updated";
				}
				delta = delta_in;
				h = e * delta + fE * Ee;
				if ( committedStress < Epo * ep + Cmax_strs ) { 
					// Within the memory surface
					Ep = Epo + h * ( delta + delta_y ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
				} else {
					// Outside the memory surface
					Ep = Epo + h * ( delta ) / 1e-15; // @todo - better idea for (delta_in-delta) than 0.0000001
				}
				if( elb != 0 ) {
					double reduction = ( 1.0 - 10 * sqrt ( W / fy ) * R );
					if ( reduction < 0.05 ) reduction = 0.05;
					Ep = Ep * reduction;
				}
				trialTangent = Ee * Ep / ( Ee + Ep );
				// Correct the overshooting of the yield surface
				trialStress = Tls_p + ( ( stressIfElastic - Tls_p ) / Ee ) * trialTangent;
				// Make sure that the loading point does not breach the bounding surface
				if ( trialStress > Epo * ep + Tbs_p )
					trialStress = Epo * ep + Tbs_p;
				dep = ( ( stressIfElastic - Tls_p ) / Ee ) * trialTangent / Ep;
			} else {
				// If continuued plastic loading
				delta = fabs( Epo * ep + Tbs_p - committedStress);
				h = e * delta + fE * Ee;
				if ( committedStress < Epo * ep + Cmax_strs ) { 
					// Within the memory surface
					Ep = Epo + h * ( delta + delta_y ) / ( delta_in - delta );
				} else {
					// Outside the memory surface
					Ep = Epo + h * ( delta ) / ( delta_in - delta );				
				}
				if( elb != 0 ) {
					double reduction = ( 1.0 - 10.0 * sqrt( W / fy ) * R );
					if ( reduction < 0.05 ) reduction = 0.05;
					Ep = Ep * reduction;
				}
				trialTangent = Ee * Ep / ( Ee + Ep );
				trialStress = committedStress + strainIncrement * trialTangent;
				// Make sure that the loading point does not breach the bounding surface
				if ( trialStress > Epo * ep + Tbs_p ) {
					trialStress = Epo * ep + Tbs_p;
					//opserr<<"Breach Bounding Surface\n";
				}
			}
		}
		

		// Update state variables
		//ep = ep + dep;
		ep = trialStrain - trialStress/Ee;
		dep = ep - Cep;
		if( ep > epmax )
			epmax = ep;
		if( ep < epmin ) 
			epmin = ep;
		if( epmax - epmin > epo)
			ebar_p = epmax - epmin;

		// Update loading surface
		if ( strainIncrement <= 0.0 ) {
			// Moving in the compressive or negative direction
			Tls_n = trialStress;
			Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 ) - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
			if( elb != 0 ) {
				double reduction = ( 1.0 - 15.0 * sqrt( W / fy ) * R );
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
		
		Tbs_n = - ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi  * pow ( 0.5 * ebar_p, 2 ) ) );
		if ( elb != 0 ) {
			Tbs_n = localBucklingBoundingStress;
		}
		Tbs_p =   ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi  * pow ( 0.5 * ebar_p, 2 ) ) );
				
		// Total work
		W = W + trialStress * dep;
		
		// Slope of the bounding line
		Epo = Epoi / ( 1 + w * W );
		
	} 
	
	// Detect local buckling if it has previously locally buckled
	if ( elb != 0 && trialStress < localBucklingBoundingStress ) {
		lb = 1;
		
		if ( elb == 1 ) {
			// stress at which local buckling occured which is new local buckling bounding stress
			localBucklingBoundingStress = committedStress;
			// make sure the buckling stress is not positive
			if ( localBucklingBoundingStress > 0.0 )
				localBucklingBoundingStress = 0.0;
			
			// residual stress
			if ( localBucklingBoundingStress >= localBucklingConstantResidualStress ) 
				localBucklingConstantResidualStress = localBucklingBoundingStress;
			
			// strain at which constant residual stress begins 
			localBucklingConstantResidualStrain = committedStrain - (localBucklingBoundingStress-localBucklingConstantResidualStress)/Ksft;	
		}
		// set some things
	}
	
	if (lb == 1) {
		if ( elb == 0 || elb == 1 ) {
			if ( trialStrain <= localBucklingConstantResidualStrain ) {
				trialStress = localBucklingConstantResidualStress;
				trialTangent = 0.0001; // @todo - better small tangent
				elb = 2;
			} else {
				trialStress = localBucklingConstantResidualStress + Ksft*(trialStrain - localBucklingConstantResidualStrain );
				trialTangent = Ksft;
				elb = 1;
			}
		} else if ( elb == 2 ) {
			trialStress = localBucklingConstantResidualStress;
			trialTangent = 0.0001; // @todo - better small tangent
		}
		localBucklingBoundingStress = trialStress;
	}
	
	// Update the maximum stress seen by the material
	if( fabs(trialStress) > Cmax_strs )
		Tmax_strs = fabs(trialStress);

	return 0;
}

double 
shenSteelCCFT::getStrain(void)
{
  return trialStrain;
}

double 
shenSteelCCFT::getStress(void)
{
  return trialStress;
}

double 
shenSteelCCFT::getInitialTangent(void)
{
	return Ee;
}

double 
shenSteelCCFT::getTangent(void)
{
  return trialTangent;
}

int 
shenSteelCCFT::commitState(void)
{
  this->commitStatevar();
  return 0;
}	


int 
shenSteelCCFT::revertToLastCommit(void)
{
  this->backtocommitStatevar();
  return 0;
}


int 
shenSteelCCFT::revertToStart(void)
{
  elastic = 1; Celastic = 1;
  plastic = 0; Cplastic = 0;
  lb = 0; Clb = 0;
  elb = 0; Celb = 0;
  loadingDirection = 0; committedLoadingDirection = 0;
  Trule = 1; Crule = 1;
  lastYieldedIn = 0; commitedLastYieldedIn = 0;
 
  Ep = 0.0; CEp = 0.0;
  Epo = 0.0; CEpo = 0.0;
  delta = 0.0; Cdelta = 0.0;
  h = 0.0; Ch = 0.0;
  trialTangent = Ee; committedTangent = Ee;
  trialStress = 0.0; committedStress = 0.0;
  trialStrain = 0.0; committedStrain = 0.0;
  
  Rls = fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) ); CRls = Rls; 
    // Set Biaxial Stress Reduction Constant
  double alphaZp = 0.5*(alphaHoop+pow(4-3*alphaHoop*alphaHoop,0.5));
  double alphaZn = 0.5*(alphaHoop-pow(4-3*alphaHoop*alphaHoop,0.5));
  Tls_p = alphaZp * Rls; Cls_p = Tls_p;
  Tls_n = alphaZn * Rls; Cls_n = Tls_n; 
  //opserr << "alphaHoop: " << alphaHoop << "  alphaZp: " << alphaZp << "  alphaZn: " << alphaZn << "\n";
  
  Tbs_p = Rbso; Cbs_p = Rbso;
  Tbs_n = -Rbso; Cbs_n = -Rbso; 
  Tmem_p = Tls_p; Cmem_p = Tls_p;
  Tmem_n = Tls_n; Cmem_n = Tls_n; 
  Tvbs_p = 0.0; Cvbs_p = 0.0;
  Tvbs_n = 0.0; Cvbs_n = 0.0; 
  ep = 0.0; Cep = 0.0; 
  epmin = 0.0; Cepmin = 0.0;
  epmax = 0.0; Cepmax = 0.0; 
  ebar_p = 0.0; Cebar_p = 0.0;
  W = 0.0; CW = 0.0; 
  Tmax_strs = 0.0; Cmax_strs = 0.0; 
  
  delta_y = 0.0; Cdelta_y = 0.0; 
  
  delta_in = 0.0; Cdelta_in = 0.0;
  
  localBucklingReferenceStrain = 0.0; committedLocalBucklingReferenceStrain = 0.0;
  
  localBucklingConstantResidualStrain = -fy/Ee; committedLocalBucklingConstantResidualStrain = localBucklingConstantResidualStrain; 
  localBucklingConstantResidualStress = -fy; committedLocalBucklingConstantResidualStress = localBucklingConstantResidualStress; 

  localBucklingBoundingStress = -fy; committedLocalBucklingBoundingStress = localBucklingBoundingStress; 
    
  return 0;
}


UniaxialMaterial *
shenSteelCCFT::getCopy(void)
{
	// @todo work on this.
  
	
  shenSteelCCFT *theCopy = new shenSteelCCFT(this->getTag(),fy,fu,Ee,R,epo,alphaHoop,localBucklingStrain,Ksft,
		  Rcrit,Rbso,Epoi,alfa,a,bb,c,w,ksi,e,fE);
  //theCopy->ep = this->ep;
 
  return theCopy;
}


int 
shenSteelCCFT::sendSelf(int cTag, Channel &theChannel)
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
//    opserr << "shenSteelCCFT::sendSelf() - failed to send data\n";
//
//  return res;
	return 1;
}

int 
shenSteelCCFT::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
	// @todo work on this.

//  int res = 0;
//  static Vector data(6);
//  res = theChannel.recvVector(this->getDbTag(), cTag, data);
//  if (res < 0) 
//    opserr << "shenSteelCCFT::recvSelf() - failed to recv data\n";
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
shenSteelCCFT::Print(OPS_Stream &s, int flag)
{
  s << "shenSteelCCFT tag: " << this->getTag() << endln;
  s << " Fy:" << fy << endln;
  s << " Fu:" << fu << endln;
  s << " Es: " << Ee << endln;
  s << " Ksft:" << Ksft << endln;
  return;
}


void
shenSteelCCFT::commitStatevar(void)
{
	  Celastic = elastic;
	  Cplastic = plastic;
	  Clb = lb;
	  Celb = elb;
	  committedLoadingDirection = loadingDirection;
	  Crule = Trule;
	  commitedLastYieldedIn = lastYieldedIn;
	 
	  CEp = Ep;
	  CEpo = Epo;
	  Cdelta = delta;
	  Ch = h;
	  committedTangent = trialTangent;
	  committedStress = trialStress;
	  committedStrain = trialStrain;
	  
	  CRls = Rls; 
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
	  Cebar_p = ebar_p;
	  CW = W; 
	  Cmax_strs = Tmax_strs; 
	  
	  Cdelta_y = delta_y; 
	  
	  Cdelta_in = delta_in;
	  
	  committedLocalBucklingReferenceStrain = localBucklingReferenceStrain;
	  
	  committedLocalBucklingConstantResidualStrain = localBucklingConstantResidualStrain; 
	  committedLocalBucklingConstantResidualStress = localBucklingConstantResidualStress; 

	  committedLocalBucklingBoundingStress = localBucklingBoundingStress; 
}	


void
shenSteelCCFT::backtocommitStatevar(void)
{
	  elastic = Celastic;
	  plastic = Cplastic;
	  lb = Clb;
	  elb = Celb;
	  loadingDirection = committedLoadingDirection;
	  Trule = Crule;
	  lastYieldedIn = commitedLastYieldedIn;
	 
	  Ep = CEp;
	  Epo = CEpo;
	  delta = Cdelta;
	  h = Ch;
	  trialTangent = committedTangent;
	  trialStress = committedStress;
	  trialStrain = committedStrain;
	  
	  Rls = CRls; 
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
	  ebar_p = Cebar_p;
	  W = CW; 
	  Tmax_strs = Cmax_strs; 
	  
	  delta_y = Cdelta_y; 
	  
	  delta_in = Cdelta_in;
	  
	  localBucklingReferenceStrain = committedLocalBucklingReferenceStrain;
	  
	  localBucklingConstantResidualStrain = committedLocalBucklingConstantResidualStrain; 
	  localBucklingConstantResidualStress = committedLocalBucklingConstantResidualStress; 

	  localBucklingBoundingStress = committedLocalBucklingBoundingStress; 
}
