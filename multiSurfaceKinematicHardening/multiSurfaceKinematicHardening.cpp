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
// $Date: 2010-05-04 17:14:45 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/multiSurfaceKinematicHardening/multiSurfaceKinematicHardening.cpp,v $

// Documentation: Mroz Multi-Surface Kinematic Hardening Material
//
// uniaxialMaterial multiSurfaceKinematicHardening $tag $InputType ...
//
// uniaxialMaterial multiSurfaceKinematicHardening $tag -Direct $E ($center_i $radius_i $hardeningModulus_i)
// uniaxialMaterial multiSurfaceKinematicHardening $tag -StressStrainSymmetric ($stress_i $strain_i) <$finalTangentStiffness>
//
// Required Input Parameters:
//   $tag           integer tag identifying material
//   $E             elastic modulus
//
// Optional Input:
//
// References:
//   1. Mróz, Z. (1967). “On the description of anisotropic workhardening.”
//      Journal of the Mechanics and Physics of Solids, 15(3), 163-175.

#include <elementAPI.h>
#include "multiSurfaceKinematicHardening.h"

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
  OPS_Error("multiSurfaceKinematicHardening unaxial material \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign\n", 1);
}

OPS_Export void *
OPS_multiSurfaceKinematicHardening()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numData;
  int sDataLength = 40;
  char *sData = new char[sDataLength];

  int tag;
  double E;
  double sigma0 = 0.0;
  int numSurfaces;
  double dData2[10];
  double * dData;
  double * startCenters;
  double * radii;
  double * hardeningModulii;

  numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial multiSurfaceKinematicHardening tag \n";
    return 0;
  }


  if ( OPS_GetString(sData, sDataLength) != 0 ) {
    opserr << "WARNING invalid input";
    return 0;
  }

  if ( strcmp(sData,"-initialStress") == 0 ) {
    numData = 1;
    if (OPS_GetDoubleInput(&numData, dData2) != 0) {
      opserr << "WARNING invalid input, want: -initialStress $sigma0 \n";
      return 0;
    }
    sigma0 = dData2[0];

    // Get the next string
    if ( OPS_GetString(sData, sDataLength) != 0 ) {
      opserr << "WARNING invalid input";
      return 0;
    }
  }

  if ( strcmp(sData,"-Direct") == 0 ) {

    numData = 1;
    if (OPS_GetDoubleInput(&numData, &E) != 0) {
      opserr << "WARNING invalid input, want: uniaxialMaterial multiSurfaceKinematicHardening tag E \n";
      return 0;
    }

    // Determine number of surfaces from number of remaining arguments
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if ( numRemainingArgs < 3 || numRemainingArgs%3 != 0 ) {
      opserr << "WARNING invalid number of arguments\n";
      return 0;
    }

    // Initialize arrays
    numSurfaces = numRemainingArgs/3;
    startCenters = new double[numSurfaces];
    radii = new double[numSurfaces];
    hardeningModulii = new double[numSurfaces];

    // Parse input data
    numData = numRemainingArgs;
    dData = new double[numData];
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid input, want: uniaxialMaterial multiSurfaceKinematicHardening tag $E $numSurfaces <$center $radius $hardeningModulus> \n";
      return 0;
    }

    for (int i=0; i<numSurfaces; i++) {
      startCenters[i] = dData[3*i];
      radii[i] = dData[3*i+1];
      hardeningModulii[i] = dData[3*i+2];
    }

  } else if ( strcmp(sData,"-StressStrainSymmetric") == 0 ) {

    // Determine number of surfaces from number of remaining arguments
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if ( numRemainingArgs < 2 ) {
      opserr << "WARNING invalid number of arguments\n";
      return 0;
    }

    // Initialize arrays
    numSurfaces = numRemainingArgs/2;
    startCenters = new double[numSurfaces];
    radii = new double[numSurfaces];
    hardeningModulii = new double[numSurfaces];

    // Parse input data
    numData = numRemainingArgs;
    dData = new double[numData];
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid input, want: uniaxialMaterial multiSurfaceKinematicHardening tag $E $numSurfaces <$center $radius $hardeningModulus> \n";
      return 0;
    }

    double * strain = new double[numSurfaces];
    double * stress = new double[numSurfaces];
    double lastTangentModulus;
    for (int i=0; i<numSurfaces; i++) {
      stress[i] = dData[2*i];
      strain[i] = dData[2*i+1];
    }
    if (numRemainingArgs%2 == 1) {
      lastTangentModulus = dData[numData-1];
    } else {
      lastTangentModulus = 0.0;
    }

    E = stress[0]/strain[0];
    for (int i=0; i<numSurfaces; i++) {
      startCenters[i] = 0.0;
      radii[i] = stress[i];
      double temp;

      if ( i == numSurfaces-1 && lastTangentModulus == 0.0) {
        hardeningModulii[i] = 0.0;
      } else {
        if ( i == numSurfaces-1 ) {
          temp = 1/lastTangentModulus - 1/E;
        } else {
          temp = (strain[i+1]-strain[i])/(stress[i+1]-stress[i]) - 1/E;
        }
        for (int j=0; j<i; j++) {
          temp -= 1/hardeningModulii[j];
        }
        hardeningModulii[i] = 1/temp;
      }
    }

    delete[] strain;
    delete[] stress;

  } else {
    opserr << "WARNING invalid InputType: "<< sData <<"\n";
    return 0;
  }




  // Create material
  theMaterial = new multiSurfaceKinematicHardening(tag, E, numSurfaces, startCenters, radii, hardeningModulii, sigma0);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type multiSurfaceKinematicHardening\n";
    return 0;
  }

  delete[] dData;
  delete[] startCenters;
  delete[] radii;
  delete[] hardeningModulii;

  return theMaterial;
}




multiSurfaceKinematicHardening::multiSurfaceKinematicHardening(int tag, double iE, int iNumSurfaces, double iCenters[], double iRadii[], double iHardeningModulii[], double iSigma0)
:UniaxialMaterial(tag,MAT_TAG_multiSurfaceKinematicHardening),
 E(iE), numSurfaces(iNumSurfaces), initStress(iSigma0)
{
  // Check input data
  if ( E <= 0.0 ) {
    opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
    opserr << "The modulus of elasticity, E, cannot be negative!\n";
    exit(-1);
  }

  if ( numSurfaces < 1 ) {
    opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
    opserr << "The number of surfaces must be greater than zero! "<<numSurfaces<<"\n";
    exit(-1);
  }

  // Initialize Arrays
  trialCenters = new double[numSurfaces];
  committedCenters = new double[numSurfaces];
  startCenters = new double[numSurfaces];
  radii = new double[numSurfaces];
  hardeningModulii = new double[numSurfaces];

  // Sset and check surface data
  for (int i=0; i<numSurfaces; i++) {
    startCenters[i] = iCenters[i];
    radii[i] = iRadii[i];
    hardeningModulii[i] = iHardeningModulii[i];

    // Check radii
    if ( radii[i] <= 0.0 ) {
      opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
      opserr << "radii cannot be negative!\n";
      exit(-1);
    }

    // Check centers
    if ( i == 0 ) {
      if ( ((startCenters[i]-radii[i]) > 0) || ((startCenters[i]+radii[i]) < 0) ) {
        opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
        opserr << "first yield surfaces must contain zero stress!\n";
        exit(-1);
      }
    } else {
      if ( ((startCenters[i]-radii[i]) > (startCenters[i-1]-radii[i-1])) || ((startCenters[i]+radii[i]) < (startCenters[i-1]+radii[i-1])) ) {
        opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
        opserr << "yield surfaces cannot overlap!\n";
        exit(-1);
      }
    }

    // Check hardening modulii
    if ( i == numSurfaces-1 ) {
      if ( hardeningModulii[i] < 0.0 ) {
        opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
        opserr << "hardeningModulii cannot be negative!\n";
        exit(-1);
      }
    } else {
      if ( hardeningModulii[i] <= 0.0 ) {
        opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
        opserr << "hardeningModulii cannot be negative!\n";
        exit(-1);
      }
    }

    // Check Initial Stress
    if ( (initStress < startCenters[0]-radii[0]) || (initStress > startCenters[0]+radii[0])  ) {
      opserr << "multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()\n";
      opserr << "initStress is outside of the elastic range!\n";
      exit(-1);
    }
  }

  // Set Remaining Variables
	this->revertToStart();
}

multiSurfaceKinematicHardening::multiSurfaceKinematicHardening()
:UniaxialMaterial(0,MAT_TAG_multiSurfaceKinematicHardening),
 E(0.0), numSurfaces(0)
{
  // Initialize Arrays
  trialCenters = new double[numSurfaces];
  committedCenters = new double[numSurfaces];
  startCenters = new double[numSurfaces];
  radii = new double[numSurfaces];
  hardeningModulii = new double[numSurfaces];

	this->revertToStart();
}

multiSurfaceKinematicHardening::~multiSurfaceKinematicHardening()
{
  // Delete Arrays
  delete[] trialCenters;
  delete[] committedCenters;
  delete[] startCenters;
  delete[] radii;
  delete[] hardeningModulii;
}

int multiSurfaceKinematicHardening::setTrialStrain(double strain, double strainRate) {

	// Define trial strain and strain increment
	int activeSurface;
	double strain_incr, yieldStress;
  if ( initStress == 0) {
    trialStrain = strain;
  } else {
    trialStrain = strain + initStress/E;
  }
  strain_incr = trialStrain - committedStrain;

	if ( strain_incr > 0 ) {
	  // Moving in the tensile direction

	  // Determine committed active surface
	  activeSurface = numSurfaces;
	  for (int i=0; i<numSurfaces; i++) {
	    if ( committedStress < committedCenters[i]+radii[i] )  {
	      activeSurface = i;
	      break;
	    }
	  }

	  // Determine trial state, including active surface
    trialTangent = tangentModulus(activeSurface);
    trialStress = committedStress+strain_incr*trialTangent;

    if ( activeSurface < numSurfaces ) {
      yieldStress =  committedCenters[activeSurface]+radii[activeSurface];

      while ( trialStress > yieldStress ) {
        strain_incr -= (yieldStress-committedStress)/trialTangent;
        activeSurface++;
        trialTangent = tangentModulus(activeSurface);
        trialStress = yieldStress + strain_incr*trialTangent;
        if ( activeSurface < numSurfaces ) {
          yieldStress = committedCenters[activeSurface]+radii[activeSurface];
        } else {
          break;
        }
      }
    }

    // Update trial surface centers
    for (int i=0; i<activeSurface; i++) {
      trialCenters[i] = trialStress - radii[i];
    }

	} else if ( strain_incr < 0) {
    // Moving in the compressive direction

    // Determine committed active surface
    activeSurface = numSurfaces;
    for (int i=0; i<numSurfaces; i++) {
      if ( committedStress > committedCenters[i]-radii[i] )  {
        activeSurface = i;
        break;
      }
    }

    // Determine trial state, including active surface
    trialTangent = tangentModulus(activeSurface);
    trialStress = committedStress+strain_incr*trialTangent;

    if ( activeSurface < numSurfaces ) {
      yieldStress =  committedCenters[activeSurface]-radii[activeSurface];

      while ( trialStress < yieldStress ) {
        strain_incr -= (yieldStress-committedStress)/trialTangent;
        activeSurface++;
        trialTangent = tangentModulus(activeSurface);
        trialStress = yieldStress + strain_incr*trialTangent;
        if ( activeSurface < numSurfaces ) {
          yieldStress = committedCenters[activeSurface]-radii[activeSurface];
        } else {
          break;
        }
      }
    }

    // Update trial surface centers
    for (int i=0; i<activeSurface; i++) {
      trialCenters[i] = trialStress + radii[i];
    }
	} else {
	  // Not moving
	  trialStress = committedStress;
	  trialTangent = committedTangent;
	}

	return 0;
}

double multiSurfaceKinematicHardening::tangentModulus(int activeSurface) {

  // Check if perfectly plastic
  if (activeSurface == numSurfaces && hardeningModulii[numSurfaces-1] == 0.0 ) {
    return 0.0;
  }

  // Determine tangent modulus
  double temp = 1/E;
  for (int i=0; i<activeSurface; i++) {
    temp += 1/hardeningModulii[i];
  }
  double Et = 1/temp;
  return Et;
}

double multiSurfaceKinematicHardening::getStrain(void) {
  if (initStress == 0.0) {
    return trialStrain;
  } else {
    return trialStrain - initStress/E;
  }
}

double multiSurfaceKinematicHardening::getStress(void) {
  return trialStress;
}


double multiSurfaceKinematicHardening::getTangent(void) {
  return trialTangent;
}

double multiSurfaceKinematicHardening::getInitialTangent(void) {
  return E;
}

int multiSurfaceKinematicHardening::commitState(void) {
  committedStrain = trialStrain;
  committedStress = trialStress;
  committedTangent = trialTangent;

  for (int i=0; i<numSurfaces; i++) {
    committedCenters[i] = trialCenters[i];
  }

  return 0;
}

int multiSurfaceKinematicHardening::revertToLastCommit(void) {
  trialStrain =  committedStrain;
  trialStress = committedStress;
  trialTangent = committedTangent;

  for (int i=0; i<numSurfaces; i++) {
    trialCenters[i] = committedCenters[i];
  }

  return 0;
}


int multiSurfaceKinematicHardening::revertToStart(void) {
  trialTangent = E;
  trialStress = initStress;
  trialStrain = initStress/E;

  for (int i=0; i<numSurfaces; i++) {
    trialCenters[i] = startCenters[i];
  }

  this->commitState();
	return 0;
}


UniaxialMaterial * multiSurfaceKinematicHardening::getCopy(void) {
  multiSurfaceKinematicHardening *theCopy = new multiSurfaceKinematicHardening(this->getTag(), E, numSurfaces, startCenters, radii, hardeningModulii, initStress);

  theCopy->trialStrain = this->trialStrain;
  theCopy->committedStrain = this->committedStrain;
  theCopy->trialStress = this->trialStress;
  theCopy->committedStress = this->committedStress;
  theCopy->trialTangent = this->trialTangent;
  theCopy->committedTangent = this->committedTangent;

  for (int i=0; i<numSurfaces; i++) {
    theCopy->trialCenters[i] = this->trialCenters[i];
    theCopy->committedCenters[i] = this->committedCenters[i];
  }

  return theCopy;
}


int multiSurfaceKinematicHardening::sendSelf(int cTag, Channel &theChannel) {
//  int res = 0;
//  static Vector data(106);
//  data(0) = this->getTag();
//  data(1) = Fc_n;
//  data(2) = ec_n;
//
//
//  res = theChannel.sendVector(this->getDbTag(), cTag, data);
//  if (res < 0)
//    opserr << "multiSurfaceKinematicHardening::sendSelf() - failed to send data\n";
//
//  return res;
  return -1;
}

int multiSurfaceKinematicHardening::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
//  int res = 0;
//  static Vector data(6);
//  res = theChannel.recvVector(this->getDbTag(), cTag, data);
//  if (res < 0)
//    opserr << "multiSurfaceKinematicHardening::recvSelf() - failed to recv data\n";
//  else {
//    this->setTag(data(0));
//
//    Fc_n = data(1);
//    ec_n = data(2);
//
//  }
//
//  return res;
  return -1;
}

void multiSurfaceKinematicHardening::Print(OPS_Stream &s, int flag) {
	s<<"multiSurfaceKinematicHardening, tag: "<<this->getTag()<<"\n";
	s<<" E: "<<E<<"\n";
	s<<" Initial Stress: "<<initStress<<"\n";
	s<<" Number of Surfaces: "<<numSurfaces<<"\n";
	for (int i=0; i<numSurfaces; i++) {
    s<<"    #"<<i+1<<" Center = "<<startCenters[i]<<" Radius = "<<radii[i]<<" Hardening Modulus = "<<hardeningModulii[i]<<"\n";
  }
	return;
}

