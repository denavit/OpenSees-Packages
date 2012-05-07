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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/sakinoSunConcrete04/sakinoSunConcrete04.cpp,v $

// Written: fmk
//
// Description: This file contains the class implementation for
// ElasticMaterial.
//
// What: "@(#) sakinoSunConcrete04.C, revA"

#include <elementAPI.h>
#include "sakinoSunConcrete04.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
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
  OPS_Error("sakinoSunConcrete04 unaxial material \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign\n", 1);
}

// Documentation: Sakino and Sun Concrete Model
// uniaxialMaterial sakinoSunConcrete04 $tag $fcc $ecc $Ec $W <$ft $et> <$beta>
//
// This is a modification of the Concrete04 material to use the Sakino & Sun model in the
// compressive region.
//
// Input Parameters:
//   $tag			integer tag identifying material
//   $fcc			peak compressive stress (input as negative value)
//   $ecc			strain at peak compressive stress (input as negative value)
//   $Ec			initial modulus of elasticity
//   $W				constant in the Sakino & Sun concrete model
//   $ft			maximum tensile strength of concrete
//   $et			ultimate tensile strain of concrete
//   $beta			the exponential curve parameter to define the residual stress (as a factor of $ft) at $etu
//
// References:
//   1. OpenSees (2010). “Concrete04 Material,” Open source software, http://opensees.berkeley.edu.
//   2. K. Sakino and Y. Sun, “Stress-strain curve of concrete confined by rectilinear hoop,”
//      Journal of Structural Construction Engineering 461 (1994): 95–104 (in Japanese).
//   3. K. Sakino et al., “Behavior of Centrally Loaded Concrete-Filled Steel-Tube Short Columns,”
//      Journal of Structural Engineering 130 (2004): 180.
//

OPS_Export void *
OPS_sakinoSunConcrete04()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];

  int    tag;
  double fcc;
  double ecc;
  double Ec;
  double W;
  double ft;
  double et;
  double beta;

  int numData;
  int remainingArgs;
  remainingArgs = OPS_GetNumRemainingInputArgs();

  if ( remainingArgs != 8 && remainingArgs != 7 && remainingArgs != 5 ) {
	  opserr << "WARNING invalid input, want: uniaxialMaterial sakinoSunConcrete04 tag fcc ecc Ec W <ft et> <beta> \n";
	  return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial sakinoSunConcrete04 tag" << endln;
    return 0;
  }
  tag = iData[0];

  numData = remainingArgs-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input, want: uniaxialMaterial sakinoSunConcrete04 tag fcc ecc Ec W <ft et> <beta> \n";
    return 0;
  }

  fcc = dData[0];
  if ( fcc >= 0 ) {
	  opserr << "WARNING fcc should be input as a negative value";
	  return 0;
  }

  ecc = dData[1];
  if ( ecc >= 0 ) {
	  opserr << "WARNING ecc should be input as a negative value";
	  return 0;
  }

  Ec = dData[2];
  if ( Ec <= 0 ) {
	  opserr << "WARNING Ec should be input as a positive value";
	  return 0;
  }

  W = dData[3];
  if ( W <= 0 ) {
	  opserr << "WARNING W should be input as a positive value";
	  return 0;
  }

  if ( remainingArgs == 8 || remainingArgs == 7 ) {
	  ft = dData[4];
	  if ( ft < 0 ) {
		  opserr << "WARNING ft should be input as a positive value";
		  return 0;
	  }

	  et = dData[5];
	  if ( et < 0 ) {
		  opserr << "WARNING et should be input as a positive value";
		  return 0;
	  }
  } else {
	  ft = 0.0;
	  et = 0.0;
  }

  if ( remainingArgs == 8 ) {
	  beta = dData[remainingArgs-2];
	  if ( beta < 0 ) {
		  opserr << "WARNING beta should be input as a positive value";
		  return 0;
	  }
  } else {
	  if ( ft != 0.0 || et != 0.0 ) {
		  beta = 0.1;
	  } else {
		  beta = 0.0;
	  }
  }

  theMaterial = new sakinoSunConcrete04(tag, fcc, ecc, Ec, W, ft, et, beta);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type sakinoSunConcrete04\n";
    return 0;
  }

  return theMaterial;
}


sakinoSunConcrete04::sakinoSunConcrete04
(int tag, double FPC, double EPSC0, double EC0, double W, double FCT, double ETU, double BETA)
  :UniaxialMaterial(tag, MAT_TAG_sakinoSunConcrete04),
   fpc(FPC), epsc0(EPSC0), Ec0(EC0), w(W), fct(FCT), etu(ETU), beta(BETA),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(FCT),
   Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0)
{
  // Make all concrete parameters negative
  if (fpc > 0.0 || epsc0 > 0.0) {
    opserr << "error: negative values required for concrete stress-strain model" << endln;
  }

  if (fct < 0.0) {
    fct = 0.0;
    opserr << "warning: fct less than 0.0 so the tensile response part is being set to 0" << endln;
  }

  Ctangent = Ec0;
  CunloadSlope = Ec0;
  CUtenSlope = Ec0;

  // Set trial values
  this->revertToLastCommit();
}


sakinoSunConcrete04::sakinoSunConcrete04():UniaxialMaterial(0, MAT_TAG_sakinoSunConcrete04),
			 fpc(0.0), epsc0(0.0), Ec0(0.0), w(0.0), fct(0.0), etu(0.0), beta(0.0),
			 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(0.0),
			 CUtenSlope(0.0), Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0)
{
  // Set trial values
  this->revertToLastCommit();

}

sakinoSunConcrete04::~sakinoSunConcrete04()
{
  // Does nothing
}

int sakinoSunConcrete04::setTrialStrain (double strain, double strainRate)
{

  /*// Reset trial history variables to last committed state*/
  TminStrain = CminStrain;
  TmaxStrain = CmaxStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  TUtenSlope = CUtenSlope;
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;

  /* // Set trial strain*/
  if (fct == 0.0 && strain > 0.0) {
    Tstrain = strain;
    Tstress = 0.0;
    Ttangent = 0.0;
    TUtenSlope = 0.0;
    return 0;
  }

  /*// Determine change in strain from last converged state*/
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;

  Tstrain = strain;

  /*// Calculate the trial state given the change in strain  // determineTrialState (dStrain);*/
  TunloadSlope = CunloadSlope;
  TUtenSlope = CUtenSlope;
  if (dStrain <= 0.0) {	  /*// Material can be either in Compression-Reloading	  // or Tension-Unloading state.*/
    if (Tstrain > 0.0) {
      /*// Material is in Tension-Unloading State*/
      Ttangent = TUtenSlope;
      Tstress = Tstrain * TUtenSlope;
    } else {
      /*// Material is in Compression-Reloading State*/
      TminStrain = CminStrain;
      TendStrain = CendStrain;
      TunloadSlope = CunloadSlope;
      CompReload();
    }
  } else {
    /*// Material can be either in Compression-Unloading	  // or Tension-Reloading State.*/
    if (Tstrain >= 0.0) {    /*// Material is in Tension-Reloading State*/
      TmaxStrain = CmaxStrain;
      if (Tstrain < TmaxStrain) {
	Tstress = Tstrain * CUtenSlope;
	Ttangent = CUtenSlope;
	TUtenSlope = CUtenSlope;
      } else {
	TmaxStrain = Tstrain;
	TensEnvelope();
	setTenUnload();
      }
    } else {
      if (Tstrain <= TendStrain) {
	Ttangent = TunloadSlope;
	Tstress = Ttangent * (Tstrain - TendStrain);
      } else {
	Tstress = 0.0;
	Ttangent = 0.0;
      }
    }
  }
  return 0;
}

void sakinoSunConcrete04::CompReload()
{
  if (Tstrain <= TminStrain) {

    TminStrain = Tstrain;

    /*// Determine point on envelope*/
    CompEnvelope ();
    setCompUnloadEnv ();

  }
  else if (Tstrain < TendStrain) {
    Ttangent = TunloadSlope;
    Tstress = Ttangent*(Tstrain-TendStrain);
  }
  else if (Tstrain <= 0.0) {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
}

void sakinoSunConcrete04::CompEnvelope()
{
  if (Tstrain >= -1.0) {

	double V = Ec0*epsc0/fpc;
    double eta = Tstrain/epsc0;
    Tstress = fpc * (V*eta+(w-1)*pow(eta,2)) / (1+(V-2)*eta+w*pow(eta,2));
    Ttangent = fpc / epsc0 * (1-eta)*(V-2*eta+V*eta+2*w*eta) / ( pow( w*pow(eta,2)-2*eta+V*eta+1 , 2) );

  } else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }

}

void sakinoSunConcrete04::setCompUnloadEnv()	{
  double tempStrain = TminStrain;

  if (tempStrain < -1.0)
    tempStrain = -1.0;

  double eta = tempStrain/epsc0;

  double ratio = 0.707*(eta-2.0) + 0.834; // unloading parameter as per Karsan-Jirsa

  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;

  TendStrain = ratio*epsc0;

  double temp1 = TminStrain - TendStrain;

  double temp2 = Tstress/Ec0;

  if (temp1 > -DBL_EPSILON) {	// temp1 should always be negative
    TunloadSlope = Ec0;
  }
  else if (temp1 <= temp2) {
    TendStrain = TminStrain - temp1;
    TunloadSlope = Tstress/temp1;
  }
  else {
    TendStrain = TminStrain - temp2;
    TunloadSlope = Ec0;
  }


  if (Tstrain >= 0.0) {
    /*opserr << "actually made it in here" << endln;*/
    /*TunloadSlope = Ec0;*/
  }

}

void sakinoSunConcrete04::TensReload()
{  TensEnvelope();  setTenUnload();}

void sakinoSunConcrete04::TensEnvelope()
{  double ect = fct / Ec0;
  if (Tstrain <= ect) {
    Tstress = Tstrain * Ec0;
    Ttangent = Ec0;
  } else if (Tstrain > etu) {
    Tstress = 0.0;
    Ttangent = 0.0;
  } else {
    Tstress = fct * pow(beta, (Tstrain - ect) / (etu - ect));
    Ttangent = fct * pow(beta, (Tstrain - ect) / (etu - ect)) * log(beta) / (etu - ect);
  }
}

void sakinoSunConcrete04::setTenUnload(){
  TUtenStress = Tstress;
  TUtenSlope = Tstress / Tstrain;
}
double sakinoSunConcrete04::getStress ()
{     return Tstress;}

double sakinoSunConcrete04::getStrain (){   return Tstrain;}

double sakinoSunConcrete04::getTangent ()
{   return Ttangent;}

int sakinoSunConcrete04::commitState ()
{
  /*// History variables*/
  CminStrain = TminStrain;
  CmaxStrain = TmaxStrain;
  CunloadSlope = TunloadSlope;
  CendStrain = TendStrain;
  CUtenSlope = TUtenSlope;

  /*// State variables*/
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  return 0;
}

int sakinoSunConcrete04::revertToLastCommit ()
{
  /*// Reset trial history variables to last committed state*/
  TminStrain = CminStrain;
  TmaxStrain = CmaxStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  TUtenSlope = CUtenSlope;

  /*// Recompute trial stress and tangent*/
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  return 0;
}

int sakinoSunConcrete04::revertToStart ()
{

  /*// History variables*/
  CminStrain = 0.0;   CmaxStrain = 0.0;   CunloadSlope = Ec0;   CendStrain = 0.0;   CUtenSlope = Ec0;
  /*// State variables*/
  Cstrain = 0.0;   Cstress = 0.0;   Ctangent = Ec0;
  /*// Reset trial variables and state*/   this->revertToLastCommit();      return 0;
}

UniaxialMaterial* sakinoSunConcrete04::getCopy ()
{
  sakinoSunConcrete04* theCopy = new sakinoSunConcrete04(this->getTag(),
				       fpc, epsc0, Ec0, w, fct, etu, beta);

  /*// Converged history variables*/
  theCopy->CminStrain = CminStrain;   theCopy->CmaxStrain = CmaxStrain;   theCopy->CunloadSlope = CunloadSlope;   theCopy->CendStrain = CendStrain;   theCopy->CUtenSlope = CUtenSlope;

  /*// Converged state variables*/   theCopy->Cstrain = Cstrain;   theCopy->Cstress = Cstress;   theCopy->Ctangent = Ctangent;

  return theCopy;
}

int sakinoSunConcrete04::sendSelf (int commitTag, Channel& theChannel)
{
  int res = 0;

  static Vector data(14);
  data(0) = this->getTag();

  /* Material properties*/
  data(1) = fpc;
  data(2) = epsc0;
  data(3) = Ec0;
  data(4) = w;
  data(5) = fct;
  /*// History variables from last converged state*/
  data(6) = CminStrain;
  data(7) = CunloadSlope;
  data(8) = CendStrain;
  data(9) = CUtenSlope;
  data(10) = CmaxStrain;

  /*// State variables from last converged state*/
  data(11) = Cstrain;
  data(12) = Cstress;
  data(13) = Ctangent;
  /*// Data is only sent after convergence, so no trial variables   // need to be sent through data vector*/

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "sakinoSunConcrete04::sendSelf() - failed to send data\n";
  return res;
}

int sakinoSunConcrete04::recvSelf (int commitTag, Channel& theChannel,
			  FEM_ObjectBroker& theBroker)
{
  int res = 0;
  static Vector data(14);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);

  if (res < 0) {
    opserr << "sakinoSunConcrete04::recvSelf() - failed to receive data\n";
    this->setTag(0);
  }
  else {
    this->setTag(int(data(0)));

    /*// Material properties */
    fpc = data(1);
    epsc0 = data(2);
    Ec0 = data(4);
    w = data(3);
    fct = data(5);

    /*// History variables from last converged state*/
    CminStrain = data(6);
    CunloadSlope = data(7);
    CendStrain = data(8);
    CcompStrain = data(9);
    CmaxStrain = data(10);

    /*// State variables from last converged state*/
    Cstrain = data(11);
    Cstress = data(12);
    Ctangent = data(13);

    /*// Set trial state variables*/
    this->revertToLastCommit();
  }

  return res;
}

void sakinoSunConcrete04::Print (OPS_Stream& s, int flag)
{
  s << "sakinoSunConcrete04, tag: " << this->getTag() << endln;
  s << "  fpc: " << fpc << endln;
  s << "  epsc0: " << epsc0 << endln;
  s << "  fct: " << fct << endln;
  s << "  w: " << w << endln;
  s << "  Ec0:  " << Ec0 << endln;
  s << "  etu:  " << etu << endln;
  s << "  beta: " << beta << endln;
}

int
sakinoSunConcrete04::getMaterialType()
{
	return 0;
}
