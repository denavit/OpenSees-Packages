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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/mixedBeamColumn3d/mixedBeamColumn3d.cpp,v $

#include "mixedBeamColumn3d.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <iostream>
#include <fstream>
#include <Node.h>
#include <Message.h>

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <RegularizedHingeIntegration.h>

// Constants that define the dimensionality
#define  NDM   3                       // dimension of the problem (3d)
#define  NND   6                       // number of nodal dof's
#define  NEGD  12                      // number of element global dof's
#define  NDM_SECTION  3                // number of section dof's without torsion
#define  NDM_NATURAL  5                // number of element dof's in the basic system without torsion
#define  NDM_NATURAL_WITH_TORSION  6   // number of element dof's in the basic system with torsion

using namespace std;


Matrix mixedBeamColumn3d::theMatrix(NEGD,NEGD);
Vector mixedBeamColumn3d::theVector(NEGD);
double mixedBeamColumn3d::workArea[400];
Matrix mixedBeamColumn3d::transformNaturalCoords(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
Matrix mixedBeamColumn3d::transformNaturalCoordsT(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
int mixedBeamColumn3d::maxNumSections = 10;

Vector *mixedBeamColumn3d::sectionDefShapeFcn = 0;
Matrix *mixedBeamColumn3d::nldhat = 0;
Matrix *mixedBeamColumn3d::nldhatT = 0;
Matrix *mixedBeamColumn3d::nd1 = 0;
Matrix *mixedBeamColumn3d::nd2 = 0;
Matrix *mixedBeamColumn3d::nd1T = 0;
Matrix *mixedBeamColumn3d::nd2T = 0;



#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif


OPS_Export void
localInit()
{
  OPS_Error("mixedBeamColumn3d element \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign, Copyright 2010\n", 1);
}

// Documentation: Three Dimensional Mixed Beam Column Element
// element mixedBeamColumn3d $tag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>
//   <-integration $intType> <-doRayleigh $rFlag>
//
// Required Input Parameters:
//   $tag                   integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts           number of integration points along the element length
//   $secTag                identifier for previously-defined section object
//   $transfTag             identifier for previously-defined coordinate-transformation (CrdTransf) object
//
// Optional Input:
//   -mass $massDens
//       $massDens          element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)
//   -integration $intType
//       $intType           numerical integration type, options are Lobotto, Legendre, Radau, NewtonCotes, Trapezoidal (optional, default= Lobotto)
//   -doRayleigh $rFlag
//       $rFlag             optional, default = 1
//                              rFlag = 0 no rayleigh damping
//                              rFlag = 1 include rayleigh damping (default)
//   -geomLinear            perform analysis without internal geometric nonlinearity
//
//
// References:
//   1. Bulent N. Alemdar and Donald W. White, "Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis," Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, "Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns," Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.
//

OPS_Export void *
OPS_mixedBeamColumn3d()
{
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  int sDataLength = 40;
  char *sData  = new char[sDataLength];
  char *sData2 = new char[sDataLength];
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != NDM) {
     opserr << "ERROR: mixedBeamColumn3d: invalid number of dimensions\n";
   return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
     opserr << "ERROR: mixedBeamColumn3d: invalid number of degrees of freedom\n";
   return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "ERROR: mixedBeamColumn3d: too few arguments\n";
    return 0;
  }

  // Get required input data
  numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - mixedBeamColumn3d\n";
    return 0;
  }
  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int numIntgrPts = iData[3];
  int secTag = iData[4];
  int transfTag = iData[5];

  // Get the section
  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(secTag);
  if (theSection == 0) {
    opserr << "WARNING section with tag " << secTag << "not found for element " << eleTag << endln;
    return 0;
  }

  SectionForceDeformation **sections = new SectionForceDeformation *[numIntgrPts];
  for (int i = 0; i < numIntgrPts; i++) {
    sections[i] = theSection;
  }

  // Get the coordinate transformation
  CrdTransf *theTransf = OPS_GetCrdTransfPtr(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
    return 0;
  }

  // Set Default Values for Optional Input
  int doRayleigh = 1;
  double massDens = 0.0;
  bool geomLinear = false;
  BeamIntegration *beamIntegr = 0;

  // Loop through remaining arguments to get optional input
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
    if ( OPS_GetString(sData, sDataLength) != 0 ) {
      opserr << "WARNING invalid input";
      return 0;
    }

    if ( strcmp(sData,"-mass") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -mass $massDens \n";
        return 0;
      }
      massDens = dData[0];

    } else if ( strcmp(sData,"-integration") == 0 ) {
      if ( OPS_GetString(sData2, sDataLength) != 0 ) {
        opserr << "WARNING invalid input, want: -integration $intType";
        return 0;
      }

      if (strcmp(sData2,"Lobatto") == 0) {
        beamIntegr = new LobattoBeamIntegration();
      } else if (strcmp(sData2,"Legendre") == 0) {
        beamIntegr = new LegendreBeamIntegration();
      } else if (strcmp(sData2,"Radau") == 0) {
        beamIntegr = new RadauBeamIntegration();
      } else if (strcmp(sData2,"NewtonCotes") == 0) {
        beamIntegr = new NewtonCotesBeamIntegration();
      } else if (strcmp(sData2,"Trapezoidal") == 0) {
        beamIntegr = new TrapezoidalBeamIntegration();
      } else if (strcmp(sData2,"RegularizedLobatto") == 0 || strcmp(sData2,"RegLobatto") == 0) {
        numData = 4;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING invalid input, want: -integration RegularizedLobatto $lpI $lpJ $zetaI $zetaJ \n";
          return 0;
        }
        BeamIntegration *otherBeamInt = 0;
        otherBeamInt = new LobattoBeamIntegration();
        beamIntegr = new RegularizedHingeIntegration(*otherBeamInt, dData[0], dData[1], dData[2], dData[3]);
          if (otherBeamInt != 0) {
            delete otherBeamInt;
          }
      } else {
        opserr << "WARNING invalid integration type, element: " << eleTag;
        return 0;
      }

    } else if ( strcmp(sData,"-doRayleigh") == 0 ) {
        numData = 1;
        if (OPS_GetInt(&numData, &doRayleigh) != 0) {
          opserr << "WARNING: Invalid doRayleigh in element mixedBeamColumn3d " << eleTag;
          return 0;
        }

    } else if ( strcmp(sData,"-geomLinear") == 0 ) {
      geomLinear = true;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }


  // Set the beam integration object if not in options
  if (beamIntegr == 0) {
    beamIntegr = new LobattoBeamIntegration();
  }

  // now create the element and add it to the Domain
  Element *theElement = new mixedBeamColumn3d(eleTag, nodeI, nodeJ, numIntgrPts, sections, *beamIntegr, *theTransf, massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theElement;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
mixedBeamColumn3d::mixedBeamColumn3d (int tag, int nodeI, int nodeJ, int numSec,
                                      SectionForceDeformation **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf,
                                      double massDensPerUnitLength,
                                      int damp, bool geomLin):
  Element(tag,ELE_TAG_mixedBeamColumn3d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0), doRayleigh(damp), geomLinear(geomLin),
  rho(massDensPerUnitLength), deflength(0.0), lengthLastIteration(0.0), lengthLastStep(0.0), initialLength(0.0),
  initialFlag(0), initialFlagB(0), itr(0), cnvg(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  internalForceOpenSees(NDM_NATURAL_WITH_TORSION), committedInternalForceOpenSees(NDM_NATURAL_WITH_TORSION),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  c(NDM_NATURAL), commitedC(NDM_NATURAL),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  kvcommit(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0), sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;

  // get copy of the beam integration object
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr<<"Error: mixedBeamColumn3d::mixedBeamColumn3d: could not create copy of beam integration object" << endln;
    exit(-1);
  }

  // get copy of the transformation object
  crdTransf = coordTransf.getCopy3d();
  if (crdTransf == 0) {
    opserr << "Error: mixedBeamColumn3d::mixedBeamColumn3d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }


  //this->setSectionPointers(numSec,sec);
  if (numSec > maxNumSections) {
    opserr << "Error: mixedBeamColumn3d::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (sec == 0) {
    opserr << "Error: mixedBeamColumn3d::setSectionPointers -- invalid section pointer";
  }

  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: mixedBeamColumn3d::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {
    if (sec[i] == 0) {
      opserr << "Error: mixedBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (SectionForceDeformation*) sec[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: mixedBeamColumn3d::setSectionPointers -- could not create copy of section " << i << endln;
    }
  }

  // Mark's Variables
  sectionForceFibers = new Vector [numSections];
  commitedSectionForceFibers = new Vector [numSections];
  sectionDefFibers = new Vector [numSections];
  commitedSectionDefFibers = new Vector [numSections];
  sectionFlexibility = new Matrix [numSections];
  commitedSectionFlexibility = new Matrix [numSections];

  for (int i = 0; i < numSections; i++){
    // Mark's Variables
    sectionForceFibers[i] = Vector(NDM_SECTION);
    commitedSectionForceFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i] = Vector(NDM_SECTION);
    commitedSectionDefFibers[i] = Vector(NDM_SECTION);
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
  }

  for(int i = 0; i < numSec; i++){
    // Mark's Variables
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i].Zero();
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i].Zero();
    sectionFlexibility[i].Zero();
    commitedSectionFlexibility[i].Zero();
  }

  V.Zero();
  committedV.Zero();
  internalForceOpenSees.Zero();
  committedInternalForceOpenSees.Zero();
  naturalForce.Zero();
  commitedNaturalForce.Zero();
  lastNaturalDisp.Zero();
  commitedLastNaturalDisp.Zero();
  c.Zero();
  commitedC.Zero();
  Hinv.Zero();
  commitedHinv.Zero();
  GMH.Zero();
  commitedGMH.Zero();
  kv.Zero();
  kvcommit.Zero();

  if ( transformNaturalCoords(1,1) != 1 ) {
    // if transformNaturalCoords hasn't been set yet then set it
    // This is needed because the way the state determination algorithm was formulated
    // and OpenSees assume a different order of natural degrees of freedom
    // Formulation --> { e, theta_zi, theta_yi, theta_zj, theta_yj, twist }
    // OpenSees    --> { e, theta_zi, theta_zj, theta_yi, theta_yj, twist }
    transformNaturalCoords.Zero();
    transformNaturalCoords(0,0) = 1;
    transformNaturalCoords(1,1) = 1;
    transformNaturalCoords(2,3) = 1;
    transformNaturalCoords(3,2) = 1;
    transformNaturalCoords(4,4) = 1;
    transformNaturalCoords(5,5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0,0) = 1;
    transformNaturalCoordsT(1,1) = 1;
    transformNaturalCoordsT(3,2) = 1;
    transformNaturalCoordsT(2,3) = 1;
    transformNaturalCoordsT(4,4) = 1;
    transformNaturalCoordsT(5,5) = 1;
  }

  if (sectionDefShapeFcn == 0)
    sectionDefShapeFcn  = new Vector [maxNumSections];
  if (nldhat == 0)
    nldhat  = new Matrix [maxNumSections];
  if (nldhatT == 0)
    nldhatT  = new Matrix [maxNumSections];
  if (nd1 == 0)
    nd1  = new Matrix [maxNumSections];
  if (nd2 == 0)
    nd2  = new Matrix [maxNumSections];
  if (nd1T == 0)
    nd1T  = new Matrix [maxNumSections];
  if (nd2T == 0)
    nd2T  = new Matrix [maxNumSections];
  if (!sectionDefShapeFcn || !nldhat || !nldhatT || !nd1 || !nd2 || !nd1T || !nd2T ) {
    opserr << "mixedBeamColumn3d::mixedBeamColumn3d() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i<maxNumSections; i++ ){
    nldhatT[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING
mixedBeamColumn3d::mixedBeamColumn3d():
  Element(0,ELE_TAG_mixedBeamColumn3d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0), doRayleigh(0), geomLinear(false),
  rho(0.0), deflength(0.0), lengthLastIteration(0.0), lengthLastStep(0.0), initialLength(0.0),
  initialFlag(0), initialFlagB(0), itr(0), cnvg(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  internalForceOpenSees(NDM_NATURAL_WITH_TORSION), committedInternalForceOpenSees(NDM_NATURAL_WITH_TORSION),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  c(NDM_NATURAL), commitedC(NDM_NATURAL),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION), kvcommit(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0), sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  sectionForceFibers = new Vector [numSections];
  commitedSectionForceFibers = new Vector [numSections];
  sectionDefFibers = new Vector [numSections];
  commitedSectionDefFibers = new Vector [numSections];
  sectionFlexibility = new Matrix [numSections];
  commitedSectionFlexibility = new Matrix [numSections];


  for (int i = 0; i < numSections; i++){
    sectionForceFibers[i] = Vector(NDM_SECTION);
    commitedSectionForceFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i] = Vector(NDM_SECTION);
    commitedSectionDefFibers[i] = Vector(NDM_SECTION);
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
  }

  for(int i = 0; i < numSections; i++){
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i].Zero();
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i].Zero();
    sectionFlexibility[i].Zero();
    commitedSectionFlexibility[i].Zero();
  }

  V.Zero();
  committedV.Zero();
  internalForceOpenSees.Zero();
  committedInternalForceOpenSees.Zero();
  naturalForce.Zero();
  commitedNaturalForce.Zero();
  lastNaturalDisp.Zero();
  commitedLastNaturalDisp.Zero();
  c.Zero();
  commitedC.Zero();
  Hinv.Zero();
  commitedHinv.Zero();
  GMH.Zero();
  commitedGMH.Zero();
  kv.Zero();
  kvcommit.Zero();

  if ( transformNaturalCoords(1,1) != 1 ) {
    // if transformNaturalCoords hasn't been set yet then set it
    transformNaturalCoords.Zero();
    transformNaturalCoords(0,0) = 1;
    transformNaturalCoords(1,1) = 1;
    transformNaturalCoords(2,3) = -1;
    transformNaturalCoords(3,2) = 1;
    transformNaturalCoords(4,4) = -1;
    transformNaturalCoords(5,5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0,0) = 1;
    transformNaturalCoordsT(1,1) = 1;
    transformNaturalCoordsT(3,2) = -1;
    transformNaturalCoordsT(2,3) = 1;
    transformNaturalCoordsT(4,4) = -1;
    transformNaturalCoordsT(5,5) = 1;
  }

  if (sectionDefShapeFcn == 0)
    sectionDefShapeFcn  = new Vector [maxNumSections];
  if (nldhat == 0)
    nldhat  = new Matrix [maxNumSections];
  if (nldhatT == 0)
    nldhatT  = new Matrix [maxNumSections];
  if (nd1 == 0)
    nd1  = new Matrix [maxNumSections];
  if (nd2 == 0)
    nd2  = new Matrix [maxNumSections];
  if (nd1T == 0)
    nd1T  = new Matrix [maxNumSections];
  if (nd2T == 0)
    nd2T  = new Matrix [maxNumSections];
  if (!sectionDefShapeFcn || !nldhat || !nldhatT || !nd1 || !nd2 || !nd1T || !nd2T ) {
    opserr << "mixedBeamColumn3d::mixedBeamColumn3d() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i<maxNumSections; i++ ){
    nldhatT[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }

}

mixedBeamColumn3d::~mixedBeamColumn3d() {

  if (sections) {
    for (int i=0; i < numSections; i++) {
      if (sections[i]) {
        delete sections[i];
      }
    }
    delete [] sections;
  }

  if (crdTransf)
    delete crdTransf;

  if (beamIntegr != 0)
    delete beamIntegr;

  if (Ki != 0)
    delete Ki;

  if (sectionForceFibers != 0)
    delete [] sectionForceFibers;

  if (commitedSectionForceFibers != 0)
    delete [] commitedSectionForceFibers;

  if (sectionDefFibers != 0)
    delete [] sectionDefFibers;

  if (commitedSectionDefFibers != 0)
    delete [] commitedSectionDefFibers;

  if (sectionFlexibility != 0)
    delete [] sectionFlexibility;

  if (commitedSectionFlexibility != 0)
    delete [] commitedSectionFlexibility;
}

int mixedBeamColumn3d::getNumExternalNodes(void) const {
   return 2;
}

const ID & mixedBeamColumn3d::getExternalNodes(void) {
   return connectedExternalNodes;
}

Node ** mixedBeamColumn3d::getNodePtrs(void) {
   return theNodes;
}

int mixedBeamColumn3d::getNumDOF(void) {
   return NEGD;
}

void mixedBeamColumn3d::setDomain(Domain *theDomain) {

  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;

    opserr << "mixedBeamColumn3d::setDomain:  theDomain = 0 ";
    exit(0);
  }

  // get pointers to the nodes
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0) {
    opserr << "mixedBeamColumn3d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }

  if (theNodes[1] == 0) {
    opserr << "mixedBeamColumn3d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }

  // call the DomainComponent class method
  this->DomainComponent::setDomain(theDomain);

  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();

  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "mixedBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }

  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "mixedBeamColumn3d::setDomain(): Error initializing coordinate transformation";
    exit(0);
  }

  deflength = crdTransf->getInitialLength();
  initialLength = deflength;
  lengthLastIteration = deflength;
  lengthLastStep = deflength;

  // Check element length
  if (deflength == 0.0) {
    opserr << "mixedBeamColumn3d::setDomain(): Zero element length:" << this->getTag();
    exit(0);
  }

}

int mixedBeamColumn3d::commitState() {
  int err = 0; // error flag
  int i = 0; // integer for loops

  cnvg = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "mixedBeamColumn3d::commitState () - failed in base class";
    return err;
  }

  // commit the sections
  do {
    err = sections[i++]->commitState();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;

  // commit the element variables state
  lengthLastStep = deflength;
  committedV = V;
  committedInternalForceOpenSees = internalForceOpenSees;
  commitedNaturalForce = naturalForce;
  commitedLastNaturalDisp = lastNaturalDisp;
  commitedC = c;
  commitedHinv = Hinv;
  commitedGMH = GMH;
  kvcommit = kv;
  for( i = 0; i < numSections; i++){
    commitedSectionForceFibers[i] = sectionForceFibers[i];
    commitedSectionDefFibers[i] = sectionDefFibers[i];
    commitedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int mixedBeamColumn3d::revertToLastCommit() {
  int err;
  int i = 0;

  cnvg = 1;

  do {
    err = sections[i]->revertToLastCommit();
    i++;
  } while (err == 0 && i < numSections);


  if (err)
    return err;

  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  // revert the element state to last commit
  V = committedV;
  internalForceOpenSees = committedInternalForceOpenSees;
  naturalForce = commitedNaturalForce;
  lastNaturalDisp = commitedLastNaturalDisp;
  c = commitedC;
  Hinv = commitedHinv;
  GMH = commitedGMH;
  kv   = kvcommit;
  for( i = 0; i < numSections; i++){
    sectionForceFibers[i] = commitedSectionForceFibers[i];
    sectionDefFibers[i] = commitedSectionDefFibers[i];
    sectionFlexibility[i] = commitedSectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int mixedBeamColumn3d::revertToStart()
{
  // revert the sections state to start
  int err;
  int i = 0;

  do {
     err = sections[i++]->revertToStart();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start
  // @todo revertToStart element state



  // Reset iteration counter
  itr = 0;

  return err;
}

const Matrix & mixedBeamColumn3d::getInitialStiff(void) {

  // Check if it has been computed already
  if (Ki != 0)
    return *Ki;

  int i,j,k; // for loops
  double GJ;

  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Vector of zeros to use at initial natural displacements
  Vector myZeros(NDM_NATURAL);
  myZeros.Zero();

  // Set initial shape functions
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, myZeros, initialLength);
    nd1[i] = this->getNd1(i, myZeros, initialLength);
    nd2[i] = this->getNd2(i, 0, initialLength);

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Set initial and section flexibility and GJ
  for ( i = 0; i < numSections; i++ ){
    Matrix ks(NDM_SECTION,NDM_SECTION);
    getSectionTangent(i,2,ks,GJ);
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
  }

  // Compute the following matrices: G, G2, H, H12, H22, Md, Kg
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();
  for( i = 0; i < numSections; i++ ){
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    // Md is zero since deformations are zero
    Kg  = Kg  + initialLength * wt[i] * this->getKg(i, 0.0, initialLength);
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta

  // Compute the transposes of the following matrices: G2, GMH
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute the stiffness matrix without the torsion term
  Matrix K_temp_noT(NDM_NATURAL,NDM_NATURAL);
  K_temp_noT = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  //K_temp_noT = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta

  // Add in the torsional stiffness term
  Matrix K_temp_withT(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
  K_temp_withT.Zero();
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j < NDM_NATURAL; j++ ) {
      K_temp_withT(i,j) = K_temp_noT(i,j);
    }
  }

  K_temp_withT(5,5) =  GJ/initialLength; // Torsional Stiffness GJ/L

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(K_temp_withT));

  return *Ki;

}

const Matrix & mixedBeamColumn3d::getTangentStiff(void) {
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Matrix ktOpenSees = transformNaturalCoordsT*kv*transformNaturalCoords;
  return crdTransf->getGlobalStiffMatrix(ktOpenSees,internalForceOpenSees);
}

const Vector & mixedBeamColumn3d::getResistingForce(void) {
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0(NDM_NATURAL);
  p0.Zero();
  return crdTransf->getGlobalResistingForce(internalForceOpenSees, p0);
}

int mixedBeamColumn3d::update() {

  int i,j,k; // integers for loops

  itr = itr + 1; // says how many times update has been called since the last commit state

  // Current and Initial Length
  crdTransf->update();
  double currentLength = crdTransf->getDeformedLength();
  const double initialLength = crdTransf->getInitialLength();

  // Compute the natural displacements
  Vector naturalDispWithTorsion = crdTransf->getBasicTrialDisp();
  naturalDispWithTorsion = transformNaturalCoords*naturalDispWithTorsion;
    // convert to the arrangement of natural deformations that the element likes

  Vector naturalDisp(NDM_NATURAL);
  for ( i = 0; i < NDM_NATURAL; i++ ) {
    naturalDisp(i) = naturalDispWithTorsion(i); //all but the torsional component
  }
  double twist = naturalDispWithTorsion(5);

  Vector naturalIncrDeltaDisp(NDM_NATURAL);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;

  // Get the numerical integration scheme
  double xi[maxNumSections]; // location of sections or gauss points or integration points
  beamIntegr->getSectionLocations(numSections, initialLength, xi);
  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Define Variables
  double temp_x, temp_A, temp_B;
  Matrix temp_Md(NDM_NATURAL,NDM_NATURAL);
  Matrix K_temp(NDM_NATURAL,NDM_NATURAL);
  double GJ;
  Vector sectionForceShapeFcn[numSections];
  for ( i = 0; i < numSections; i++ ) {
    sectionForceShapeFcn[i] = Vector(NDM_SECTION);
  }

  // If update had never been called before, set variables to initial status
  // This should be in revertToStart
  if ( initialFlagB == 0 ) {

      // Compute initial H11, G, Md, H12, matrices and c vector
      Vector myZeros(NDM_NATURAL);
      myZeros.Zero();

        for ( i = 0; i < numSections; i++ ){
           nldhat[i] = this->getNld_hat(i, myZeros, lengthLastStep);
          sectionDefShapeFcn[i] = this->getd_hat(i, myZeros, lengthLastStep);
           nd1[i] = this->getNd1(i, myZeros, lengthLastStep);
           nd2[i] = this->getNd2(i, 0, lengthLastStep);

           for( j = 0; j < NDM_SECTION; j++ ){
              for( k = 0; k < NDM_NATURAL; k++ ){
                nd1T[i](k,j) = nd1[i](j,k);
                nd2T[i](k,j) = nd2[i](j,k);
                nldhatT[i](k,j) = nldhat[i](j,k);
              }
           }
        }

       // Set initial and committed section flexibility
        for ( i = 0; i < numSections; i++ ){
        Matrix ks(NDM_SECTION,NDM_SECTION);
        getSectionTangent(i,2,ks,GJ);
        invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
          commitedSectionFlexibility[i] = sectionFlexibility[i];
        }

        for ( i = 0; i < numSections; i++ ){
          sectionForceFibers[i].Zero();
          commitedSectionForceFibers[i].Zero();
          sectionDefFibers[i].Zero();
          commitedSectionDefFibers[i].Zero();
        }

        // Compute the H matrix and its inverse
        Matrix H(NDM_NATURAL,NDM_NATURAL);
        H.Zero();
        for( i = 0; i < numSections; i++ ){
           H = H + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
        }
        invertMatrix(NDM_NATURAL, H, Hinv);
        commitedHinv = Hinv;

        // Compute the H12 matrix
        Matrix H12(NDM_NATURAL,NDM_NATURAL);
        H12.Zero();
        for( i = 0; i < numSections; i++ ){
           H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
        }

        // Compute the G matrix and its transpose
        Matrix G(NDM_NATURAL,NDM_NATURAL);
        G.Zero();
        for( i = 0; i < numSections; i++ ){
           G = G + initialLength * wt[i] * nd1T[i] * nldhat[i];
        }

        Matrix Md(NDM_NATURAL,NDM_NATURAL);
        Md.Zero();	// should be zero at converged state

        GMH = G + Md - H12;
        //GMH = G; // Omit P-small delta

        commitedGMH = GMH;
        c.Zero();	// should be zero at initial state
       commitedC.Zero();
        naturalForce.Zero();
        commitedNaturalForce.Zero();

        initialFlagB = 1; // Don't need to come back here again
   }


  // Compute shape functions and their transposes
  for ( i = 0; i < numSections; i++ ){
    // Shape Functions
    nldhat[i] = this->getNld_hat(i, naturalDisp, currentLength);
    sectionDefShapeFcn[i] = this->getd_hat(i, naturalDisp, currentLength);
    nd1[i] = this->getNd1(i, naturalDisp, currentLength);
    nd2[i] = this->getNd2(i, internalForceOpenSees(0), currentLength);

    // Transpose of shape functions
    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
        nldhatT[i](k,j) = nldhat[i](j,k);
      }
    }
  }

  // Update natural force
  //naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Alemdar's Scheme
  naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V ); // Nukala's Scheme
  //naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Omit P-small delta

  // Update sections
  for ( i = 0; i < numSections; i++){
    // Compute section deformations
    sectionForceShapeFcn[i] = nd1[i] * naturalForce;
    sectionDefFibers[i] = sectionDefFibers[i] + sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );

    // Send section deformation to section object
    double twist = 0.0; // set torsional strain to zero
    setSectionDeformation(i,sectionDefFibers[i],twist);

    // Get section force vector
    double torsionalForce;
    getSectionStress(i,sectionForceFibers[i],torsionalForce);

    // Get section tangent matrix
    Matrix ks(NDM_SECTION,NDM_SECTION);
    getSectionTangent(i,1,ks,GJ);

    // Compute section flexibility matrix
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
  }

  // Compute the following matrices: V, c, V2, G, G2, H, H12, H22, Md, Kg
  Vector V2(NDM_NATURAL);
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

  V.Zero();
  c.Zero();
  V2.Zero();
  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();

  for( i = 0; i < numSections; i++ ){
    V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
    c   = c   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
    V2  = V2  + initialLength * wt[i] * nd2T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    Kg  = Kg  + initialLength * wt[i] * this->getKg(i, sectionForceFibers[i](0), currentLength);
      // sectionForceFibers[i](0) is the axial load, P

      temp_x = currentLength * xi[i];
      temp_A =  ( temp_x / currentLength - 2 * pow( temp_x / currentLength, 2 ) + pow( temp_x / currentLength, 3 ) ) * currentLength;
      temp_B =  ( - pow( temp_x / currentLength, 2 ) + pow( temp_x / currentLength, 3 ) ) * currentLength;
      temp_Md.Zero();
      temp_Md(0,1) = temp_A * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
      temp_Md(0,2) = temp_A * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
      temp_Md(0,3) = temp_B * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
      temp_Md(0,4) = temp_B * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
    Md  = Md  + initialLength * wt[i] * temp_Md;
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta

  // Compute the transposes of the following matrices: G, G2, GMH
  Matrix GT(NDM_NATURAL,NDM_NATURAL);
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      GT(i,j) = G(j,i);
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute new internal force
  Vector internalForce(NDM_NATURAL);
  internalForce.Zero();

  //internalForce = GT * naturalForce + V2 + GMHT * Hinv * c; // Alemdar's Scheme
  internalForce = GT * naturalForce + V2 + GMHT * Hinv * V; // Nukala's Scheme
  //internalForce = GT * naturalForce + GMHT * Hinv * c; // Omit P-small delta


  // Compute internal force for OpenSees ( i.e., add torsion and rearrange )
  for ( i = 0; i < NDM_NATURAL; i++ ) {
    internalForceOpenSees(i) = internalForce(i);
  }
  internalForceOpenSees(5) = twist*GJ/currentLength; // Add in torsional force
  internalForceOpenSees = transformNaturalCoordsT*internalForceOpenSees;

  // Compute the stiffness matrix without the torsion term
  K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  //K_temp = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta

  // Add in the torsional stiffness term
  kv.Zero();
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j< NDM_NATURAL; j++ ) {
      kv(i,j) = K_temp(i,j);
    }
  }

  kv(5,5) =  GJ/currentLength; // Torsional Stiffness GJ/L

  return 0;
}

const Matrix & mixedBeamColumn3d::getMass(void) {
  theMatrix.Zero();

  if (rho != 0.0) {
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
    theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*initialLength*rho;
  }

  return theMatrix;
}


const Vector & mixedBeamColumn3d::getResistingForceIncInertia() {

  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Add the inertial forces
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(2) += m*accel1(2);
    theVector(6) += m*accel2(0);
    theVector(7) += m*accel2(1);
    theVector(8) += m*accel2(2);
  }

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theVector += this->getRayleighDampingForces();
  }

  return theVector;
}


void mixedBeamColumn3d::Print(OPS_Stream &s, int flag) {

  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho;
    for (int i = 0; i < numSections; i++)
      s << "\nSection "<<i<<" :" << *sections[i];

  } else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
    double xi[maxNumSections]; // location of sections or gauss points or integration points
    beamIntegr->getSectionLocations(numSections, initialLength, xi);
    double wt[maxNumSections]; // weights of sections or gauss points of integration points
    beamIntegr->getSectionWeights(numSections, initialLength, wt);
    s << "\n section xi wt";
    for (int i = 0; i < numSections; i++)
      s << "\n"<<i<<" "<<xi[i]<<" "<<wt[i];

  } else {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
  }

}


OPS_Stream &operator<<(OPS_Stream &s, mixedBeamColumn3d &E) {
  E.Print(s);
  return s;
}


Response* mixedBeamColumn3d::setResponse(const char **argv, int argc,
                                         OPS_Stream &output) {

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ForceBeamColumn3d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"force") == 0 ||
      strcmp(argv[0],"globalForce") == 0 ||
      strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Pz_1");
    output.tag("ResponseType","Mx_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","Mx_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse = new ElementResponse(this, 1, theVector);

  } else if (strcmp(argv[0],"section") ==0) {
    if (argc > 2) {

      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number",sectionNum);
        output.attr("eta",xi[sectionNum-1]*L);

        theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);

        output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int mixedBeamColumn3d::getResponse(int responseID, Information &eleInfo) {
  if (responseID == 1) { // global forces
    return eleInfo.setVector(this->getResistingForce());

  } else {
    return -1;

  }
}

Vector mixedBeamColumn3d::getd_hat(int sec, const Vector &v, double L) {
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

  Vector D_hat(NDM_SECTION);
  D_hat.Zero();

  temp_x = L * xi[sec];
  temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
  temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
  temp_C = 1 / L;
  temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
  temp_E = - 4 / L + 6 * temp_x / ( L * L );
  temp_F =  - 2 / L + 6 * temp_x / ( L * L );

  D_hat(0) =  temp_C * v(0) +
              0.5 * ( temp_C * temp_C * v(0) ) * v(0) +
              0.5 * ( temp_A * temp_A * v(1) + temp_A * temp_B * v(3) ) * v(1) +
              0.5 * ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) ) * v(2) +
              0.5 * ( temp_A * temp_B * v(1) + temp_B * temp_B * v(3) ) * v(3) +
              0.5 * ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) ) * v(4);

  D_hat(1) =  temp_E * v(1) + temp_F * v(3);

  D_hat(2) =  temp_E * v(2) + temp_F * v(4);

  return D_hat;
}

Matrix mixedBeamColumn3d::getKg(int sec, double P, double L) {
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B;

  temp_x = L * xi[sec];

  Matrix kg(NDM_NATURAL,NDM_NATURAL);
  kg.Zero();

  temp_A = 1 - 4 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );
  temp_B = - 2 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );

  kg(0,0) = P / ( L * L );
  kg(1,1) = P * temp_A * temp_A;
  kg(1,3) = P * temp_A * temp_B;
  kg(2,2) = P * temp_A * temp_A;
  kg(2,4) = P * temp_A * temp_B;
  kg(3,1) = P * temp_A * temp_B;
  kg(3,3) = P * temp_B * temp_B;
  kg(4,2) = P * temp_A * temp_B;
  kg(4,4) = P * temp_B * temp_B;

  return kg;
}


Matrix mixedBeamColumn3d::getNld_hat(int sec, const Vector &v, double L) {
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

  temp_x = L * xi[sec];
  temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
  temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
  temp_C =  1 / L;
  temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
  temp_E = - 4 / L + 6 * temp_x / ( L * L );
  temp_F =  - 2 / L + 6 * temp_x / ( L * L );

  Matrix Nld_hat(NDM_SECTION,NDM_NATURAL);
  Nld_hat.Zero();

  Nld_hat(0,0) = temp_C + temp_C * temp_C * v(0);
  Nld_hat(0,1) =  temp_A * temp_A * v(1) + temp_A * temp_B * v(3);
  Nld_hat(0,2) =  temp_A * temp_A * v(2) + temp_A * temp_B * v(4);
  Nld_hat(0,3) =  temp_A * temp_B * v(1) + temp_B * temp_B * v(3);
  Nld_hat(0,4) =  temp_A * temp_B * v(2) + temp_B * temp_B * v(4);
  Nld_hat(1,1) = temp_E;
  Nld_hat(1,3) = temp_F;
  Nld_hat(2,2) = temp_E;
  Nld_hat(2,4) = temp_F;

  return Nld_hat;
}

Matrix mixedBeamColumn3d::getNd2(int sec, double P, double L) {
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B;

  temp_x = L * xi[sec];

  Matrix Nd2(NDM_SECTION,NDM_NATURAL);
  Nd2.Zero();

  temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
  temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

  Nd2(1,1) = P * temp_A;
  Nd2(1,3) = P * temp_B;
  Nd2(2,2) = P * temp_A;
  Nd2(2,4) = P * temp_B;

  return Nd2;
}

Matrix mixedBeamColumn3d::getNd1(int sec, const Vector &v, double L) {
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B;

  temp_x = L * xi[sec];

  temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[1]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[3];

  temp_B = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[4];

  Matrix Nd1(NDM_SECTION,NDM_NATURAL);
  Nd1.Zero();

  Nd1(0,0)   = 1.0;
  Nd1(1,0)   = temp_A;
  Nd1(1,1)   = - temp_x / L + 1.0;
  Nd1(1,3)   = temp_x / L;
  Nd1(2,0)   = temp_B;
  Nd1(2,2)   = - temp_x / L + 1.0;
  Nd1(2,4)   = temp_x / L;

  return Nd1;
}

void mixedBeamColumn3d::getSectionTangent(int sec,int type,Matrix &kSection,
                                          double &GJ) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize formulation friendly variables
  kSection.Zero();
  GJ = 0.0;

  // Get the stress resultant from section
  Matrix sectionTangent(order,order);
  if ( type == 1 ) {
    sectionTangent = sections[sec]->getSectionTangent();
  } else if ( type == 2 ) {
    sectionTangent = sections[sec]->getInitialTangent();
  } else {
    sectionTangent.Zero();
  }

  // Set Components of Section Tangent
  int i,j;
  for (i = 0; i < order; i++) {
    for (j = 0; j < order; j++) {
      switch(code(i)) {
        case SECTION_RESPONSE_P:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(0,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(0,1) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MY:
              kSection(0,2) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_MZ:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(1,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(1,1) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MY:
              kSection(1,2) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_MY:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(2,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(2,1) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MY:
              kSection(2,2) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_T:
          GJ = sectionTangent(i,i);
          break;
        default:
          break;
      }
    }
  }
}

void mixedBeamColumn3d::getSectionStress(int sec,Vector &fSection,
                                         double &torsion) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Get the stress resultant from section
  Vector stressResultant = sections[sec]->getStressResultant();

  // Initialize formulation friendly variables
  fSection.Zero();
  torsion = 0.0;

  // Set Components of Section Stress Resultant
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        fSection(0) = stressResultant(j);
        break;
      case SECTION_RESPONSE_MZ:
        fSection(1) = stressResultant(j);
        break;
      case SECTION_RESPONSE_MY:
        fSection(2) = stressResultant(j);
        break;
      case SECTION_RESPONSE_T:
        torsion = stressResultant(j);
        break;
      default:
        break;
    }
  }
}

void mixedBeamColumn3d::setSectionDeformation(int sec,Vector &defSection,
                                              double &twist) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize Section Deformation Vector
  Vector sectionDeformation(order);
  sectionDeformation.Zero();

  // Set Components of Section Deformations
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        sectionDeformation(j) = defSection(0);
        break;
      case SECTION_RESPONSE_MZ:
        sectionDeformation(j) = defSection(1);
        break;
      case SECTION_RESPONSE_MY:
        sectionDeformation(j) = defSection(2);
        break;
      case SECTION_RESPONSE_T:
        sectionDeformation(j) = twist;
        break;
      default:
        break;
    }
  }

  // Set the section deformations
  int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}


int mixedBeamColumn3d::sendSelf(int commitTag, Channel &theChannel) {
  // @todo write mixedBeamColumn3d::sendSelf
  opserr << "Error: mixedBeamColumn3d::sendSelf -- not yet implemented for mixedBeamColumn3d element";
  return -1;
}

int mixedBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker) {
  // @todo write mixedBeamColumn3d::recvSelf
  opserr << "Error: mixedBeamColumn3d::sendSelf -- not yet implemented for mixedBeamColumn3d element";
  return -1;
}
