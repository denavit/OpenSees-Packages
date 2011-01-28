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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/mixedBeamColumn2d/mixedBeamColumn2d.cpp,v $
                                                                        
#include "mixedBeamColumn2d.h"
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
#define  NDM   2                      // dimension of the problem (2d)
#define  NND   3                      // number of nodal dof's
#define  NEGD  6                      // number of element global dof's
#define  NDM_SECTION  2               // number of section dof's without torsio
#define  NDM_NATURAL  3               // number of element dof's in the basic system without torsion
#define  NDM_NATURAL_WITH_TORSION  3  // number of element dof's in the basic system with torsion

using namespace std;


Matrix mixedBeamColumn2d::theMatrix(NEGD,NEGD);
Vector mixedBeamColumn2d::theVector(NEGD);
double mixedBeamColumn2d::workArea[400];
Matrix mixedBeamColumn2d::transformNaturalCoords(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
Matrix mixedBeamColumn2d::transformNaturalCoordsT(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
int mixedBeamColumn2d::maxNumSections = 10;

Vector *mixedBeamColumn2d::sectionDefShapeFcn = 0;
Matrix *mixedBeamColumn2d::nldhat = 0;
Matrix *mixedBeamColumn2d::nldhatT = 0;
Matrix *mixedBeamColumn2d::nd1 = 0;
Matrix *mixedBeamColumn2d::nd2 = 0;
Matrix *mixedBeamColumn2d::nd1T = 0;
Matrix *mixedBeamColumn2d::nd2T = 0;



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
  OPS_Error("mixedBeamColumn2d element \nWritten by Mark D Denavit, University of Illinois at Urbana-Champaign, Copyright 2010\n", 1);
}

// Documentation: Two Dimensional Mixed Beam Column Element
// element mixedBeamColumn2d $tag $iNode $jNode $numIntgrPts $secTag $transfTag
//
// Required Input Parameters:
//   $tag					integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts 			number of integration points along the element length
//   $secTag 				identifier for previously-defined section object
//   $transfTag   			identifier for previously-defined coordinate-transformation (CrdTransf) object
//
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
//   1. Bulent N. Alemdar and Donald W. White, “Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis,” Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, “Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns,” Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.
//

OPS_Export void *
OPS_mixedBeamColumn2d()
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
     opserr << "ERROR: mixedBeamColumn2d: invalid number of dimensions\n";
	 return 0;	  
  }
  
  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
     opserr << "ERROR: mixedBeamColumn2d: invalid number of degrees of freedom\n";
	 return 0;	  
  }
  
  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 6) {
	  opserr << "ERROR: mixedBeamColumn2d: too few arguments\n";
	  return 0;	  
  }
  
  numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - mixedBeamColumn2d\n";
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
  for (int i = 0; i < numIntgrPts; i++)
    sections[i] = theSection;
  
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
		    if (otherBeamInt != 0)
		      delete otherBeamInt;
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
  if (beamIntegr == 0)
    beamIntegr = new LobattoBeamIntegration();
  
  
  // now create the element and add it to the Domain
  Element *theElement = new mixedBeamColumn2d(eleTag, nodeI, nodeJ, numIntgrPts, sections, *beamIntegr, *theTransf, massDens, doRayleigh, geomLinear);
  
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theElement;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
mixedBeamColumn2d::mixedBeamColumn2d (int tag, int nodeI, int nodeJ, int numSec, SectionForceDeformation **sec,
				      BeamIntegration &bi, CrdTransf &coordTransf, double massDensPerUnitLength, int damp, bool geomLin):
  Element(tag,ELE_TAG_mixedBeamColumn2d), 
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
  kv(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION), kvcommit(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
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
   if(beamIntegr == 0) {
     opserr<<"Error: mixedBeamColumn2d::mixedBeamColumn2d: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   crdTransf = coordTransf.getCopy2d();
   if (crdTransf == 0) {
      opserr << "Error: mixedBeamColumn2d::mixedBeamColumn2d: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }

   
   //this->setSectionPointers(numSec,sec);
   if (numSec > maxNumSections) {
     opserr << "Error: mixedBeamColumn2d::setSectionPointers -- max number of sections exceeded";
   }

   numSections = numSec;

   if (sec == 0) {
     opserr << "Error: mixedBeamColumn2d::setSectionPointers -- invalid section pointer";
   }

   sections = new SectionForceDeformation *[numSections];
   if (sections == 0) {
     opserr << "Error: mixedBeamColumn2d::setSectionPointers -- could not allocate section pointers";
   }

   for (int i = 0; i < numSections; i++) {

     if (sec[i] == 0) {
       opserr << "Error: mixedBeamColumn2d::setSectionPointers -- null section pointer " << i << endln;
     }

     sections[i] = (SectionForceDeformation*) sec[i]->getCopy();

     if (sections[i] == 0) {
       opserr << "Error: mixedBeamColumn2d::setSectionPointers -- could not create copy of section " << i << endln;
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
	   transformNaturalCoords.Zero();
	   transformNaturalCoords(0,0) = 1;
	   transformNaturalCoords(1,1) = 1;
	   transformNaturalCoords(2,2) = 1;
	   transformNaturalCoordsT.Zero();
	   transformNaturalCoordsT(0,0) = 1;
	   transformNaturalCoordsT(1,1) = 1;
	   transformNaturalCoordsT(2,2) = 1;
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
     opserr << "mixedBeamColumn2d::mixedBeamColumn2d() -- failed to allocate static section arrays";   
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

mixedBeamColumn2d::mixedBeamColumn2d():
  Element(0,ELE_TAG_mixedBeamColumn2d), 
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
	   transformNaturalCoords(2,2) = 1;
	   transformNaturalCoordsT.Zero();
	   transformNaturalCoordsT(0,0) = 1;
	   transformNaturalCoordsT(1,1) = 1;
	   transformNaturalCoordsT(2,2) = 1;
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
    opserr << "mixedBeamColumn2d::mixedBeamColumn2d() -- failed to allocate static section arrays";   
    exit(-1);
  }
  
  int i;
  for ( i=0; i<maxNumSections; i++ ){
	  nldhatT[i] = Matrix(NDM_NATURAL,NDM_SECTION);
	  nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
	  nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }
  
}

mixedBeamColumn2d::~mixedBeamColumn2d()
{

   if (sections) {
      for (int i=0; i < numSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if (crdTransf)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;
   
   if (Ki != 0)
     delete Ki;
   
   if(sectionForceFibers != 0)
     delete [] sectionForceFibers;
   
   if(commitedSectionForceFibers != 0)
     delete [] commitedSectionForceFibers;

   if(sectionDefFibers != 0)
     delete [] sectionDefFibers;
   
   if(commitedSectionDefFibers != 0)
     delete [] commitedSectionDefFibers;
   
   if(sectionFlexibility != 0)
     delete [] sectionFlexibility;
   
   if(commitedSectionFlexibility != 0)
     delete [] commitedSectionFlexibility;
}


int
mixedBeamColumn2d::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
mixedBeamColumn2d::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
mixedBeamColumn2d::getNodePtrs(void)
{
   return theNodes;
}

int
mixedBeamColumn2d::getNumDOF(void)
{
   return NEGD;
}


void
mixedBeamColumn2d::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;

      opserr << "mixedBeamColumn2d::setDomain:  theDomain = 0 ";
      exit(0);
   }

   // get pointers to the nodes

   int Nd1 = connectedExternalNodes(0);
   int Nd2 = connectedExternalNodes(1);

   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);

   if (theNodes[0] == 0)
   {
      opserr << "mixedBeamColumn2d::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "mixedBeamColumn2d::setDomain: Nd2: ";
      opserr << Nd2 << "does not exist in model\n";
      exit(0);
   }

   // call the DomainComponent class method
   this->DomainComponent::setDomain(theDomain);

   // ensure connected nodes have correct number of dof's
   int dofNode1 = theNodes[0]->getNumberDOF();
   int dofNode2 = theNodes[1]->getNumberDOF();

   if ((dofNode1 != NND) || (dofNode2 != NND))
   {
      opserr << "mixedBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "mixedBeamColumn2d::setDomain(): Error initializing coordinate transformation";
      exit(0);
   }

   deflength = crdTransf->getInitialLength();
   initialLength = deflength;
   lengthLastIteration = deflength;
   lengthLastStep = deflength;
   
   // Check element length
   if (deflength == 0.0)
   {
      opserr << "mixedBeamColumn2d::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

}

int
mixedBeamColumn2d::commitState()
{
   int err = 0; // error flag
   int i = 0; 
   int j = 0; // integers for loops

   cnvg = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "mixedBeamColumn2d::commitState () - failed in base class";
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
   itr = 0;
   
   return err;
}


int mixedBeamColumn2d::revertToLastCommit()
{
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
   itr = 0;

   return err;
}


int mixedBeamColumn2d::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;

   do {
       err = sections[i++]->revertToStart();

   }while (err == 0 && i < numSections);

   if (err)
      return err;

   // revert the transformation to start
   if ((err = crdTransf->revertToStart()) != 0)
      return err;

   // revert the element state to start
   // @todo revertToStart element state
//   kvcommit.Zero();
//   lengthLastIteration = crdTransf->getInitialLength();
//   lengthLastStep = crdTransf->getInitialLength();
//
//   kv.Zero();
//
//   // THis needs work
//   
//   initialFlag = 0;
//   // this->update();
   
   
   return err;
}

const Matrix &
mixedBeamColumn2d::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;
  
  int i,j,k; // for loops 
  
  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);
  
  // Vector of zeros to use at inital natural dispalcements
  Vector myZeros(NDM_NATURAL);
  myZeros.Zero();
 
  // Set inital shape functions 
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
	getSectionTangent(i,2,ks);
	invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);

	//Matrix sectionTangent = sections[i]->getInitialTangent();
  	//invertMatrix(NDM_SECTION,sectionTangent,sectionFlexibility[i]);
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
   	Kg  = Kg  + initialLength * wt[i] * this->getKg(i, 0.0, initialLength);
  	// Md is zero since deformations are zero
  }
	
  
  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
	
  // Compute the transposes of the following matricies: G2, GMH
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
   	for( j = 0; j < NDM_NATURAL; j++ ){
    		G2T(i,j) = G2(j,i);
    		GMHT(i,j) = GMH(j,i);
  	}
  }  
  
  // Compute the stiffness matrix without the torsion term
  Matrix K_temp(NDM_NATURAL,NDM_NATURAL);
  K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  //K_temp = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta
  
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(K_temp));
  
  return *Ki;

}

const Matrix &
mixedBeamColumn2d::getTangentStiff(void){
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Matrix ktOpenSees = transformNaturalCoordsT*kv*transformNaturalCoords;
  return crdTransf->getGlobalStiffMatrix(ktOpenSees,internalForceOpenSees);
}

const Vector &
mixedBeamColumn2d::getResistingForce(void){
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0(NDM_NATURAL);
  p0.Zero();
  return crdTransf->getGlobalResistingForce(internalForceOpenSees, p0);	 
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int mixedBeamColumn2d::update()
{
//  if( cnvg == 1 ){ // I think this prevents update from running right after a revertToLastCommit
//    cnvg = 0;
//    return 0;
//  }
  
  int i,j,k,m; // integers for loops

  itr = itr + 1; // says how many times update has been called since the last commmit state 

  //mixedBeamColumn2dSectionStiffness<<"\n Iteration: "<<itr<<" Element: "<<this->getTag()<<endl;

  crdTransf->update();
  double currentLength = crdTransf->getDeformedLength();
  const double initialLength = crdTransf->getInitialLength();
  
  if( initialFlag == 2 ) { // i don't think this is ever called
	  this->revertToLastCommit();
  }

  Vector naturalDisp = crdTransf->getBasicTrialDisp();   
  naturalDisp = transformNaturalCoords*naturalDisp; 
    
  Vector naturalIncrDeltaDisp(NDM_NATURAL);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;
  
  double xi[maxNumSections]; // location of sections or gauss points or intergration points 
  beamIntegr->getSectionLocations(numSections, initialLength, xi);

  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  	double temp_x, temp_A, temp_B;
  	Matrix temp_Md(NDM_NATURAL,NDM_NATURAL);
  	Matrix K_temp(NDM_NATURAL,NDM_NATURAL);

  	Vector sectionForceShapeFcn[numSections];
  	Vector sectionDefWithTorsion(NDM_SECTION);
  	Vector sectionForceFibersWithTorsion(NDM_SECTION);
  	
	for ( i = 0; i < numSections; i++ )
		sectionForceShapeFcn[i] = Vector(NDM_SECTION);
		
  	if ( initialFlagB == 0 ) {
  		
  		// Compute initial H11, G, Md, H12, matricies and c vector
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
  	  	  	getSectionTangent(i,2,ks);
  	  	  	invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);

  	  		//Matrix sectionTangent = sections[i]->getInitialTangent();
  	  		//invertMatrix(NDM_SECTION,sectionTangent,sectionFlexibility[i]);

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
     	nldhat[i] = this->getNld_hat(i, naturalDisp, currentLength);
    	sectionDefShapeFcn[i] = this->getd_hat(i, naturalDisp, currentLength);
     	nd1[i] = this->getNd1(i, naturalDisp, currentLength);
     	nd2[i] = this->getNd2(i, internalForceOpenSees(0), currentLength);

     	for( j = 0; j < NDM_SECTION; j++ ){
        	for( k = 0; k < NDM_NATURAL; k++ ){
        		nd1T[i](k,j) = nd1[i](j,k);
        		nd2T[i](k,j) = nd2[i](j,k);
        		nldhatT[i](k,j) = nldhat[i](j,k);
        	}
     	}
  	}  	

  	
	//naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Alemdar's Scheme
	naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V ); // Nukala's Scheme
	//naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Omit P-small delta

	for ( i = 0; i < numSections; i++){
		// Compute section deformations
		sectionForceShapeFcn[i] = nd1[i] * naturalForce;
		sectionDefFibers[i] = sectionDefFibers[i] + sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );

		// Send section deformation to section object
		setSectionDeformation(i,sectionDefFibers[i]);

		// Get section force vector
		getSectionStress(i,sectionForceFibers[i]);

		// Get section tangent matrix
		Matrix ks(NDM_SECTION,NDM_SECTION);
		getSectionTangent(i,1,ks);

		// Compute section flexibility matrix
		invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
  			
	}

  	// Compute the following matricies: G, G2, V, c, V2, H, H12, H22, Md, Kg
  	Matrix G(NDM_NATURAL,NDM_NATURAL);
  	Matrix G2(NDM_NATURAL,NDM_NATURAL);
  	Matrix H(NDM_NATURAL,NDM_NATURAL);
  	Matrix H12(NDM_NATURAL,NDM_NATURAL);
  	Matrix H22(NDM_NATURAL,NDM_NATURAL);
  	Matrix Md(NDM_NATURAL,NDM_NATURAL);
  	Matrix Kg(NDM_NATURAL,NDM_NATURAL);
  	Vector V2(NDM_NATURAL);
  	
  	G.Zero();
  	G2.Zero();
  	V.Zero();
  	c.Zero();
  	V2.Zero();
  	H.Zero();
  	H12.Zero();
  	H22.Zero();
  	Md.Zero();
  	Kg.Zero();
  	for( i = 0; i < numSections; i++ ){
       	G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
       	G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
       	V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
       	c   = c   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
       	V2  = V2  + initialLength * wt[i] * nd2T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
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
       		temp_Md(0,2) = temp_B * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
     	Md  = Md  + initialLength * wt[i] * temp_Md;
  	}
  	
 	// Compute the inverse of the H matrix
  	invertMatrix(NDM_NATURAL, H, Hinv);

  	// Compute the GMH matrix ( G + Md - H12 ) and its transpose
  	GMH = G + Md - H12;
  	//GMH = G; // Omit P-small delta
  	
  	// Compute the transposes of the following matricies: G, G2, GMH
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
  	
 	
  	// Define the internal force
  	Vector internalForce(NDM_NATURAL);
  	internalForce.Zero();
  		
	// Compute new internal force
	//internalForce = GT * naturalForce + V2 + GMHT * Hinv * c; // Alemdar's Scheme
	internalForce = GT * naturalForce + V2 + GMHT * Hinv * V; // Nukala's Scheme
	//internalForce = GT * naturalForce + GMHT * Hinv * c; // Omit P-small delta
  	
    
    // Compute internal force for OpenSees ( i.e., add torsion and rearrange )
    for( i = 0; i < NDM_NATURAL; i++ )
    	internalForceOpenSees(i) = internalForce(i);
    	
    internalForceOpenSees = transformNaturalCoordsT*internalForceOpenSees;
  	
    // Compute the stiffness matrix without the torsion term
    K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
    //K_temp = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta
    
    // Add in the torsional stiffness term
    kv.Zero();
    for( i = 0; i < NDM_NATURAL; i++ )
    	for( j = 0; j< NDM_NATURAL; j++ )
    		kv(i,j) = K_temp(i,j);
    
  	return 0;
}

const Matrix &
mixedBeamColumn2d::getMass(void)
{
  theMatrix.Zero();

  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) =
	theMatrix(3,3) = theMatrix(4,4) = 0.5*initialLength*rho;

  return theMatrix;
}


const Vector &
mixedBeamColumn2d::getResistingForceIncInertia()
{
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


void
mixedBeamColumn2d::Print(OPS_Stream &s, int flag)
{
   if (flag == 1) {
      s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn2d ";
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
      s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn2d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, mixedBeamColumn2d &E)
{
    E.Print(s);
    return s;
}


Response*
mixedBeamColumn2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

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
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, theVector);

    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","V_1");
      output.tag("ResponseType","M_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","V_2");
      output.tag("ResponseType","M_2");

      theResponse = new ElementResponse(this, 2, theVector);

    } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","M_1");
        output.tag("ResponseType","M_2");

        theResponse = new ElementResponse(this, 3, Vector(3));

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


int 
mixedBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) { // global forces
	  return eleInfo.setVector(this->getResistingForce());
  } else if (responseID == 2) { // local forces
	    theVector(3) =  internalForceOpenSees(0);
	    theVector(0) = -internalForceOpenSees(0);
	    theVector(2) = internalForceOpenSees(1);
	    theVector(5) = internalForceOpenSees(2);
	    double V = (internalForceOpenSees(1)+internalForceOpenSees(2))/crdTransf->getInitialLength();
	    theVector(1) =  V;
	    theVector(4) = -V;
	    return eleInfo.setVector(theVector);
  } else if (responseID == 4) { // basic forces
	  return eleInfo.setVector(internalForceOpenSees);
  } else {
    return -1;
  
  }
}

Vector
mixedBeamColumn2d::getd_hat(int sec, const Vector &v, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

   Vector D_hat(NDM_SECTION);
   D_hat.Zero();

   temp_x = L * xi[sec];
   temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_C = 1 / L;
   //temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );

   D_hat(0) =  temp_C * v(0) + 0.5 * temp_A * temp_A * v(1) * v(1) +
               temp_A * temp_B * v(1) * v(2) + 0.5 * temp_B * temp_B * v(2) * v(2);  

   D_hat(1) =  temp_E * v(1) + temp_F * v(2);

   return D_hat;
}

Matrix
mixedBeamColumn2d::getKg(int sec, double P, double L){
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
   kg(2,2) = P * temp_B * temp_B;

   kg(1,2) = P * temp_A * temp_B;
   kg(2,1) = P * temp_A * temp_B;

   return kg;
}


Matrix
mixedBeamColumn2d::getNld_hat(int sec, const Vector &v, double L){
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
   Nld_hat(0,1) = temp_A * temp_A * v(1) + temp_A * temp_B * v(2);
   Nld_hat(0,2) = temp_A * temp_B * v(1) + temp_B * temp_B * v(2);
   Nld_hat(1,1) = temp_E;
   Nld_hat(1,2) = temp_F;

   return Nld_hat;
}

Matrix
mixedBeamColumn2d::getNd2(int sec, double P, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   Matrix Nd2(NDM_SECTION,NDM_NATURAL);
   Nd2.Zero();

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
   temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

   Nd2(1,1) = P * temp_A;
   Nd2(1,2) = P * temp_B;

   return Nd2;
}


Matrix
mixedBeamColumn2d::getNd1(int sec, const Vector &v, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A;

   temp_x = L * xi[sec];

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[1]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2];

   Matrix Nd1(NDM_SECTION,NDM_NATURAL);
   Nd1.Zero();

   Nd1(0,0)   = 1.0;
   Nd1(1,0)   = temp_A;
   Nd1(1,1)   = - temp_x / L + 1.0;
   Nd1(1,2)   = temp_x / L;

   return Nd1;
}
	

void
mixedBeamColumn2d::getSectionTangent(int sec,int type,Matrix &kSection) {
    int order = sections[sec]->getOrder();
    const ID &code = sections[sec]->getType();

    // Initialize formulation friendly variables
    kSection.Zero();

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
			  default:
				break;
			  }
			  break;
		  default:
			  break;
		  }
		}
    }
}

void
mixedBeamColumn2d::getSectionStress(int sec,Vector &fSection) {
    int order = sections[sec]->getOrder();
    const ID &code = sections[sec]->getType();

    // Get the stress resultant from section
    Vector stressResultant = sections[sec]->getStressResultant();

    // Initialize formulation friendly variables
    fSection.Zero();

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
      default:
        break;
      }
    }
}

void
mixedBeamColumn2d::setSectionDeformation(int sec,Vector &defSection) {
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
      default:
        break;
      }
    }

    // Set the section deformations
    int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}


int
mixedBeamColumn2d::sendSelf(int commitTag, Channel &theChannel){  
  // @todo write mixedBeamColumn2d::sendSelf
  opserr << "Error: mixedBeamColumn2d::sendSelf -- not yet implemented for mixedBeamColumn2d element";
  return -1;
}
    
int
mixedBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker){
  // @todo write mixedBeamColumn2d::recvSelf
  opserr << "Error: mixedBeamColumn2d::sendSelf -- not yet implemented for mixedBeamColumn2d element";
  return -1;
}
