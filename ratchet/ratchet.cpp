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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/ratchet/ratchet.cpp,v $

// Documentation: Ratchet Material
//
// uniaxialMaterial ratchet $tag $direction $E
//
// Required Input Parameters:
//   $tag           integer tag identifying material
//   $direction     direction of allowed movement
//   $E             modulus in direction opposite of allowed movement
//
// Optional Input:
//
// References:
//   1.
//

#include <elementAPI.h>
#include "ratchet.h"

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
  OPS_Error("ratchet unaxial material \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign\n", 1);
}

OPS_Export void *
OPS_ratchet()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numData;
  int sDataLength = 40;
  char *sData = new char[sDataLength];

  int tag;
  double E;
  double direction;
  int numSurfaces;

  numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ratchet tag \n" << endln;
    return 0;
  }


  if ( OPS_GetStringCopy(&sData) != 0 ) {
    opserr << "WARNING invalid input";
    return 0;
  }

  if ( strcmp(sData,"tension") == 0 ) {
    direction = 1.0;
  } else if ( strcmp(sData,"compression") == 0 ) {
    direction = -1.0;
  } else {
    opserr << "WARNING invalid direction: " << sData << "\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &E) != 0) {
    opserr << "WARNING invalid input, want: uniaxialMaterial ratchet tag direction E \n";
    return 0;
  }

  // Create material
  theMaterial = new ratchet(tag, direction, E);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ratchet\n";
    return 0;
  }

  return theMaterial;
}




ratchet::ratchet(int tag, double iDirection, double iE)
:UniaxialMaterial(tag,MAT_TAG_ratchet),
 direction(iDirection), E(iE)
{
  // Check input data

  // Set Remaining Variables
	this->revertToStart();
}

ratchet::ratchet()
:UniaxialMaterial(0,MAT_TAG_ratchet),
 direction(0.0), E(0.0)
{
	this->revertToStart();
}

ratchet::~ratchet()
{

}

int ratchet::setTrialStrain(double strain, double strainRate) {

	// Define trial strain and strain increment
	trialStrain = strain;

	double relativeStrain = direction*trialStrain - committedEp;

	if ( relativeStrain < 0.0 ) {
	  trialStress = relativeStrain*E;
	  trialTangent = E;
	} else {
    trialStress = 0.0;
    trialTangent = 0.0;
    ep = direction*trialStrain;
	}

	return 0;
}

double ratchet::getStrain(void) {
  return trialStrain;
}

double ratchet::getStress(void) {
  return trialStress;
}


double ratchet::getTangent(void) {
  return trialTangent;
}

double ratchet::getInitialTangent(void) {
  return E;
}

int ratchet::commitState(void) {
  committedStrain = trialStrain;
  committedStress = trialStress;
  committedTangent = trialTangent;

  committedEp = ep;

  return 0;
}

int ratchet::revertToLastCommit(void) {
  trialStrain =  committedStrain;
  trialStress = committedStress;
  trialTangent = committedTangent;

  ep = committedEp;

  return 0;
}


int ratchet::revertToStart(void) {
  trialStrain = 0.0;
  trialStress = 0.0;
  trialTangent = 0.0;

  ep = 0.0;

  this->commitState();
	return 0;
}


UniaxialMaterial * ratchet::getCopy(void) {
  ratchet *theCopy = new ratchet(this->getTag(), direction, E);

  theCopy->trialStrain = this->trialStrain;
  theCopy->committedStrain = this->committedStrain;
  theCopy->trialStress = this->trialStress;
  theCopy->committedStress = this->committedStress;
  theCopy->trialTangent = this->trialTangent;
  theCopy->committedTangent = this->committedTangent;

  theCopy->ep = this->ep;
  theCopy->committedEp = this->committedEp;

  return theCopy;
}


int ratchet::sendSelf(int cTag, Channel &theChannel) {
//  int res = 0;
//  static Vector data(106);
//  data(0) = this->getTag();
//  data(1) = Fc_n;
//  data(2) = ec_n;
//
//
//  res = theChannel.sendVector(this->getDbTag(), cTag, data);
//  if (res < 0)
//    opserr << "ratchet::sendSelf() - failed to send data\n";
//
//  return res;
  return -1;
}

int ratchet::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
//  int res = 0;
//  static Vector data(6);
//  res = theChannel.recvVector(this->getDbTag(), cTag, data);
//  if (res < 0)
//    opserr << "ratchet::recvSelf() - failed to recv data\n";
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

void ratchet::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s<<"ratchet, tag: "<<this->getTag()<<endln;
    s<<" Forward Direction: "<<direction<<endln;
    s<<" E: "<<E<<endln;
  } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ratchet\", ";
    s << "\"direction\": " << direction << ", ";
    s << "\"E\": " << E << "}";
  }
	return;
}

