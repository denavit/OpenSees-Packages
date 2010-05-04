
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
                                                                        
// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign 
//
// Description: This file contains the implementation for the mixedBeamColumn3d class.
//
// What: "@(#) mixedBeamColumn3d.C, revA"


// we specify what header files we need
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


#define  NDM   3         // dimension of the problem (3d)
#define  NL    2         // size of uniform load vector
#define  NND   6         // number of nodal dof's
#define  NEGD  12        // number of element global dof's
#define  NUM_SECTION_DISPLACEMENTS_WITH_TORSION  4
#define  NUM_NATURAL_DISPLACEMENTS_WITH_TORSION  6 // number of element dof's in the basic system
#define  NUM_SECTION_DISPLACEMENTS  3
#define  NUM_NATRUAL_DISPLACEMENTS  5

#define  NUKALA_NN_ALGORITHM 1  // code for Nukala's NN algorithm
#define  NUKALA_LL_ALGORITHM 2  // code for Nukala's LL algorithm


using namespace std;


Matrix mixedBeamColumn3d::theMatrix(NEGD,NEGD);
Vector mixedBeamColumn3d::theVector(NEGD);
double mixedBeamColumn3d::workArea[400];
Matrix mixedBeamColumn3d::transformNaturalCoords(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION);
Matrix mixedBeamColumn3d::transformNaturalCoordsT(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION);
int mixedBeamColumn3d::maxNumSections = 10;

Vector *mixedBeamColumn3d::sectionDefShapeFcn = 0;
Matrix *mixedBeamColumn3d::ks = 0;
Matrix *mixedBeamColumn3d::fs = 0;
Matrix *mixedBeamColumn3d::ksa = 0;
Matrix *mixedBeamColumn3d::fsa = 0;
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
  OPS_Error("mixedBeamColumn3d element \nWritten by Mark D Denavit, University of Illinois at Urbana-Champaign, Copyright 2010\n", 1);
}

OPS_Export void *
OPS_mixedBeamColumn3d()
{
  // get the id and end nodes 
  int iData[6];
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
  
  numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - mixedBeamColumn3d\n";
    return 0;
  }

  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int numIntgrPts = iData[3];

  // Get the section
  int secTag = iData[4];
  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(secTag);
  
  if (theSection == 0) {
    opserr << "WARNING section with tag " << secTag << "not found for element " << eleTag << endln;
    return 0;
  }

  FiberSectionGJ **sections = new FiberSectionGJ *[numIntgrPts];
  for (int i = 0; i < numIntgrPts; i++)
    sections[i] = (FiberSectionGJ*) theSection;
  
  
  // Get the coordinate transformation
  int transfTag = iData[5];
  CrdTransf3d *theTransf3d = OPS_GetGetCrdTransf3d(transfTag);

  if (theTransf3d == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
    return 0;
  }

  
  LobattoBeamIntegration beamIntegr;
  
  double massDens = 0.0;
  
  // now create the element and add it to the Domain
  Element *theElement = new mixedBeamColumn3d(eleTag, nodeI, nodeJ, numIntgrPts, sections, beamIntegr, *theTransf3d, massDens);
  
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theElement;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
mixedBeamColumn3d::mixedBeamColumn3d (int tag, int nodeI, int nodeJ, int numSec, FiberSectionGJ **sec,
				      BeamIntegration &bi, CrdTransf3d &coordTransf, double massDensPerUnitLength):
  Element(tag,ELE_TAG_mixedBeamColumn3d), 
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(massDensPerUnitLength), deflength(0.0), lengthLastIteration(0.0), lengthLastStep(0.0), initialLength(0.0),
  initialFlag(0), initialFlagB(0), itr(0), cnvg(0),
  V(NUM_NATRUAL_DISPLACEMENTS), committedV(NUM_NATRUAL_DISPLACEMENTS), 
  internalForceOpenSees(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), committedInternalForceOpenSees(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), 
  naturalForce(NUM_NATRUAL_DISPLACEMENTS), commitedNaturalForce(NUM_NATRUAL_DISPLACEMENTS), 
  lastNaturalDisp(NUM_NATRUAL_DISPLACEMENTS), commitedLastNaturalDisp(NUM_NATRUAL_DISPLACEMENTS), 
  c(NUM_NATRUAL_DISPLACEMENTS), commitedC(NUM_NATRUAL_DISPLACEMENTS), 
  Hinv(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), commitedHinv(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), 
  GMH(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), commitedGMH(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), 
  kv(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), kvcommit(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), 
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
     opserr<<"Error: mixedBeamColumn3d::mixedBeamColumn3d: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   crdTransf = coordTransf.getCopy();
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

   sections = new FiberSectionGJ *[numSections];
   if (sections == 0) {
     opserr << "Error: mixedBeamColumn3d::setSectionPointers -- could not allocate section pointers";
   }

   for (int i = 0; i < numSections; i++) {

     if (sec[i] == 0) {
       opserr << "Error: mixedBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
     }

     sections[i] = (FiberSectionGJ*) sec[i]->getCopy();

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
         sectionForceFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
         commitedSectionForceFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
         sectionDefFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
         commitedSectionDefFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
         sectionFlexibility[i] = Matrix(NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
         commitedSectionFlexibility[i] = Matrix(NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
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
   if (ks == 0) 
	   ks  = new Matrix [maxNumSections];
   if (fs == 0) 
	   fs  = new Matrix [maxNumSections];
   if (ksa == 0) 
	   ksa  = new Matrix [maxNumSections];
   if (fsa == 0) 
	   fsa  = new Matrix [maxNumSections];
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
   if (!sectionDefShapeFcn || !ks || !fs || !ksa || !fsa || !nldhat || !nldhatT || !nd1 || !nd2 || !nd1T || !nd2T ) {
     opserr << "mixedBeamColumn3d::mixedBeamColumn3d() -- failed to allocate static section arrays";   
     exit(-1);
   }
   
   int i;
   for ( i=0; i<maxNumSections; i++ ){
 	  nldhatT[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
 	  nd1T[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
 	  nd2T[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
   }
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

mixedBeamColumn3d::mixedBeamColumn3d():
  Element(0,ELE_TAG_mixedBeamColumn3d), 
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(0.0), deflength(0.0), lengthLastIteration(0.0), lengthLastStep(0.0), initialLength(0.0),
  initialFlag(0), initialFlagB(0), itr(0), cnvg(0),
  V(NUM_NATRUAL_DISPLACEMENTS), committedV(NUM_NATRUAL_DISPLACEMENTS), 
  internalForceOpenSees(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), committedInternalForceOpenSees(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), 
  naturalForce(NUM_NATRUAL_DISPLACEMENTS), commitedNaturalForce(NUM_NATRUAL_DISPLACEMENTS), 
  lastNaturalDisp(NUM_NATRUAL_DISPLACEMENTS), commitedLastNaturalDisp(NUM_NATRUAL_DISPLACEMENTS), 
  c(NUM_NATRUAL_DISPLACEMENTS), commitedC(NUM_NATRUAL_DISPLACEMENTS), 
  Hinv(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), commitedHinv(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), 
  GMH(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), commitedGMH(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS), 
  kv(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), kvcommit(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION), 
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
        sectionForceFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
        commitedSectionForceFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
        sectionDefFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
        commitedSectionDefFibers[i] = Vector(NUM_SECTION_DISPLACEMENTS);
        sectionFlexibility[i] = Matrix(NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
        commitedSectionFlexibility[i] = Matrix(NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
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
  if (ks == 0) 
	   ks  = new Matrix [maxNumSections];
  if (fs == 0) 
	   fs  = new Matrix [maxNumSections];
  if (ksa == 0) 
	   ksa  = new Matrix [maxNumSections];
  if (fsa == 0) 
	   fsa  = new Matrix [maxNumSections];
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
  if (!sectionDefShapeFcn || !ks || !fs || !ksa || !fsa || !nldhat || !nldhatT || !nd1 || !nd2 || !nd1T || !nd2T ) {
    opserr << "mixedBeamColumn3d::mixedBeamColumn3d() -- failed to allocate static section arrays";   
    exit(-1);
  }
  
  int i;
  for ( i=0; i<maxNumSections; i++ ){
	  nldhatT[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
	  nd1T[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
	  nd2T[i] = Matrix(NUM_NATRUAL_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
  }
  
}

mixedBeamColumn3d::~mixedBeamColumn3d()
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
mixedBeamColumn3d::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
mixedBeamColumn3d::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
mixedBeamColumn3d::getNodePtrs(void)
{
   return theNodes;
}

int
mixedBeamColumn3d::getNumDOF(void)
{
   return NEGD;
}


void
mixedBeamColumn3d::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
   if (theDomain == 0)
   {
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

   if (theNodes[0] == 0)
   {
      opserr << "mixedBeamColumn3d::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "mixedBeamColumn3d::setDomain: Nd2: ";
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
      opserr << "mixedBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "mixedBeamColumn3d::setDomain(): Error initializing coordinate transformation";
      exit(0);
   }

   deflength = crdTransf->getInitialLength();
   initialLength = deflength;
   lengthLastIteration = deflength;
   lengthLastStep = deflength;
   
   // Check element length
   if (deflength == 0.0)
   {
      opserr << "mixedBeamColumn3d::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

}

int
mixedBeamColumn3d::commitState()
{
   int err = 0; // error flag
   int i = 0; 
   int j = 0; // integers for loops

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
   itr = 0;
   
   return err;
}


int mixedBeamColumn3d::revertToLastCommit()
{
	//ofstream nnOutput;
	//nnOutput.open("nnOutput.dat",ios::app);
	//nnOutput<<"\nrevertToLastCommit()\n";
	
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


int mixedBeamColumn3d::revertToStart()
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
mixedBeamColumn3d::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;
  
  int i,j,k; // for loops 
  double GJ;
  
  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);
  
  // Vector of zeros to use at inital natural dispalcements
  Vector myZeros(NUM_NATRUAL_DISPLACEMENTS);
  myZeros.Zero();
 
  // Set inital shape functions 
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, myZeros, initialLength);
    nd1[i] = this->getNd1(i, myZeros, initialLength);
    nd2[i] = this->getNd2(i, 0, initialLength);

    for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ ){
      for( k = 0; k < NUM_NATRUAL_DISPLACEMENTS; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
   	}
  }  	
  
  // Set inital and section flexibility and GJ
  for ( i = 0; i < numSections; i++ ){
	Matrix sectionTangentWithTorsion = sections[i]->getInitialTangent();
	
	if ( i == 0 )
		GJ = sectionTangentWithTorsion(NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);
	
  	Matrix sectionFlexibilityWithTorsion;
  	invertMatrix(NUM_SECTION_DISPLACEMENTS_WITH_TORSION,sectionTangentWithTorsion,sectionFlexibilityWithTorsion);
		
  	// fs and ks are the section tangent and flexibility matricies without the GJ term (size = NUM_SECTION_DISPLACEMENTS)
  	for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ )
  	  for( k = 0; k < NUM_SECTION_DISPLACEMENTS; k++ )
  	    sectionFlexibility[i](j,k) = sectionFlexibilityWithTorsion(j,k);
  }
  
  // Compute the following matricies: G, G2, H, H12, H22, Md, Kg
  Matrix G(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix G2(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix H(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix H12(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix H22(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix Md(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix Kg(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
	
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
   	
//    	temp_x = currentLength * xi[i];
//    	temp_A =  ( temp_x / currentLength - 2 * pow( temp_x / currentLength, 2 ) + pow( temp_x / currentLength, 3 ) ) * currentLength;
//    	temp_B =  ( - pow( temp_x / currentLength, 2 ) + pow( temp_x / currentLength, 3 ) ) * currentLength;
//    	temp_Md.Zero();
//    	temp_Md(0,1) = temp_A * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
//    	temp_Md(0,2) = temp_A * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
//    	temp_Md(0,3) = temp_B * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
//    	temp_Md(0,4) = temp_B * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
//   	Md  = Md  + initialLength * wt[i] * temp_Md;
   	// Md is zero since deformations are zero
  }
	
  
  // Compute the inverse of the H matrix
  invertMatrix(NUM_NATRUAL_DISPLACEMENTS, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
	
  // Compute the transposes of the following matricies: G2, GMH
  Matrix G2T(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  Matrix GMHT(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  for( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ ){
   	for( j = 0; j < NUM_NATRUAL_DISPLACEMENTS; j++ ){
    		G2T(i,j) = G2(j,i);
    		GMHT(i,j) = GMH(j,i);
  	}
  }  
  
  // Compute the stiffness matrix without the torsion term
  Matrix K_temp_noT(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  K_temp_noT = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  //K_temp_noT = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta
  
  // Add in the torsional stiffness term
  Matrix K_temp_withT(NUM_NATURAL_DISPLACEMENTS_WITH_TORSION,NUM_NATURAL_DISPLACEMENTS_WITH_TORSION);
  K_temp_withT.Zero();
  for( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ )
  	for( j = 0; j < NUM_NATRUAL_DISPLACEMENTS; j++ )
  		K_temp_withT(i,j) = K_temp_noT(i,j);
  
  K_temp_withT(5,5) =  GJ/initialLength; // Torsional Stiffness GJ/L  
    
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(K_temp_withT));
  
  return *Ki;

}

const Matrix &
mixedBeamColumn3d::getTangentStiff(void){
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Matrix ktOpenSees = transformNaturalCoordsT*kv*transformNaturalCoords;
  return crdTransf->getGlobalStiffMatrix(ktOpenSees,internalForceOpenSees);
}

const Vector &
mixedBeamColumn3d::getResistingForce(void){
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0(NUM_NATRUAL_DISPLACEMENTS);
  p0.Zero();
  return crdTransf->getGlobalResistingForce(internalForceOpenSees, p0);	 
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int mixedBeamColumn3d::update()
{
//  if( cnvg == 1 ){ // I think this prevents update from running right after a revertToLastCommit
//    cnvg = 0;
//    return 0;
//  }
  
  int i,j,k,m; // integers for loops

  itr = itr + 1; // says how many times update has been called since the last commmit state 

  //mixedBeamColumn3dSectionStiffness<<"\n Iteration: "<<itr<<" Element: "<<this->getTag()<<endl;

  crdTransf->update();
  double currentLength = crdTransf->getDeformedLength();
  const double initialLength = crdTransf->getInitialLength();
  
  if( initialFlag == 2 ) { // i don't think this is ever called
	  this->revertToLastCommit();
  }

  Vector naturalDispWithTorsion = crdTransf->getBasicTrialDisp();   
  naturalDispWithTorsion = transformNaturalCoords*naturalDispWithTorsion; 
  	// convert to the arrangement of natural deformations that the element likes
    
  Vector naturalDisp(NUM_NATRUAL_DISPLACEMENTS);
  for ( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ )
	  naturalDisp(i) = naturalDispWithTorsion(i); //remove the torsional component
  
  double twist = naturalDispWithTorsion(5);
  
  Vector naturalIncrDeltaDisp(NUM_NATRUAL_DISPLACEMENTS);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;
  
  double xi[maxNumSections]; // location of sections or gauss points or intergration points 
  beamIntegr->getSectionLocations(numSections, initialLength, xi);

  double wt[maxNumSections]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  	double temp_x, temp_A, temp_B;
  	Matrix temp_Md(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix K_temp(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	double GJ;
  	Vector sectionForceShapeFcn[numSections];
  	Vector sectionDefWithTorsion(NUM_SECTION_DISPLACEMENTS_WITH_TORSION);
  	Vector sectionForceFibersWithTorsion(NUM_SECTION_DISPLACEMENTS_WITH_TORSION);
  	
	for ( i = 0; i < numSections; i++ )
		sectionForceShapeFcn[i] = Vector(NUM_SECTION_DISPLACEMENTS);
		
  	int algorithmType = NUKALA_LL_ALGORITHM;


  	if ( initialFlagB == 0 ) {
  		
  		// Compute initial H11, G, Md, H12, matricies and c vector
  		Vector myZeros(NUM_NATRUAL_DISPLACEMENTS);
  		myZeros.Zero();
  		
  	  	for ( i = 0; i < numSections; i++ ){
  	     	nldhat[i] = this->getNld_hat(i, myZeros, lengthLastStep);
  	    	sectionDefShapeFcn[i] = this->getd_hat(i, myZeros, lengthLastStep);
  	     	nd1[i] = this->getNd1(i, myZeros, lengthLastStep);
  	     	nd2[i] = this->getNd2(i, 0, lengthLastStep);

  	     	for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ ){
  	        	for( k = 0; k < NUM_NATRUAL_DISPLACEMENTS; k++ ){
  	        		nd1T[i](k,j) = nd1[i](j,k);
  	        		nd2T[i](k,j) = nd2[i](j,k);
  	        		nldhatT[i](k,j) = nldhat[i](j,k);
  	        	}
  	     	}
  	  	}  	
  	  	
     	// Set inital and commited section flexibility
  	  	for ( i = 0; i < numSections; i++ ){
   	  		Matrix sectionTangentWithTorsion = sections[i]->getInitialTangent();
  	  		Matrix sectionFlexibilityWithTorsion;
  	  		invertMatrix(NUM_SECTION_DISPLACEMENTS_WITH_TORSION,sectionTangentWithTorsion,sectionFlexibilityWithTorsion);
     		
  	  		// fs and ks are the section tangent and flexibility matricies without the GJ term (size = NUM_SECTION_DISPLACEMENTS)
  	  		for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ )
  	  			for( k = 0; k < NUM_SECTION_DISPLACEMENTS; k++ )
  	  				sectionFlexibility[i](j,k) = sectionFlexibilityWithTorsion(j,k);

  	  		commitedSectionFlexibility[i] = sectionFlexibility[i];
  	  	}
  	  	
  	  	for ( i = 0; i < numSections; i++ ){
  	  		sectionForceFibers[i].Zero();
  	  		commitedSectionForceFibers[i].Zero();
  	  		sectionDefFibers[i].Zero();
  	  		commitedSectionDefFibers[i].Zero();
  	  	}
  	  	
  	  	// Compute the H matrix and its inverse
  	  	Matrix H(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	  	H.Zero();
  	  	for( i = 0; i < numSections; i++ ){
  	     	H = H + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
  	  	}
  	  	invertMatrix(NUM_NATRUAL_DISPLACEMENTS, H, Hinv);
  	  	commitedHinv = Hinv;
  	  	
  	  	// Compute the H12 matrix
  	  	Matrix H12(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	  	H12.Zero();
  	  	for( i = 0; i < numSections; i++ ){
  	     	H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
  	  	}

  	  	// Compute the G matrix and its transpose
  	  	Matrix G(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	  	G.Zero();
  	  	for( i = 0; i < numSections; i++ ){
       		G = G + initialLength * wt[i] * nd1T[i] * nldhat[i];
  	  	}
  	  	
  	  	Matrix Md(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
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

     	for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ ){
        	for( k = 0; k < NUM_NATRUAL_DISPLACEMENTS; k++ ){
        		nd1T[i](k,j) = nd1[i](j,k);
        		nd2T[i](k,j) = nd2[i](j,k);
        		nldhatT[i](k,j) = nldhat[i](j,k);
        	}
     	}
  	}  	

  	
  	if ( algorithmType == NUKALA_LL_ALGORITHM ) {

  	  	//naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Alemdar's Scheme
  	 	naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V ); // Nukala's Scheme
  	  	//naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Omit P-small delta
  		
  		for ( i = 0; i < numSections; i++){
  			
  			sectionForceShapeFcn[i] = nd1[i] * naturalForce;
  			sectionDefFibers[i] = sectionDefFibers[i] + sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );
  			
  			for( m = 0; m < NUM_SECTION_DISPLACEMENTS; m++ )
  				sectionDefWithTorsion(m) = sectionDefFibers[i](m);
  			
  			sectionDefWithTorsion(3) = 0.0; // set torsional strain to zero
     		     			
  			// sectionDeformations is the total strain which is sent to the section
  			int res = sections[i]->setTrialSectionDeformation(sectionDefWithTorsion);
  			
  			// sectionForceFibersWithTorsion is the total stress
  			sectionForceFibersWithTorsion = sections[i]->getStressResultant();
  			
  			// Get section tangent and flexibility matricies with the GJ term (size = NUM_SECTION_DISPLACEMENTS_WITH_TORSION)
  			ksa[i] = sections[i]->getSectionTangent();
  			invertMatrix(NUM_SECTION_DISPLACEMENTS_WITH_TORSION,ksa[i],fsa[i]);
  			GJ = ksa[i](3,3);
  			
  			// fs and ks are the section tangent and flexibility matricies without the GJ term (size = NUM_SECTION_DISPLACEMENTS)
  			for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ )
  				for( k = 0; k< NUM_SECTION_DISPLACEMENTS; k++ )
  					sectionFlexibility[i](j,k) = fsa[i](j,k);
  			
  			for( m = 0; m < NUM_SECTION_DISPLACEMENTS; m++ )
  				sectionForceFibers[i](m) = sectionForceFibersWithTorsion(m);
  			
  		}
  		
  	} else if ( algorithmType == NUKALA_NN_ALGORITHM ) {
  		
  		Vector sectionForceResidual(NUM_SECTION_DISPLACEMENTS_WITH_TORSION);
  		int iSectionIter, iElementIter;
  		int maxSectionIter = 10;
  		int maxElementIter = 10;
  		double tolSectionForce = 0.001;
  		double tolV = 0.001;
  		
  		// compute V
  		V.Zero();
  		for( i = 0; i < numSections; i++ )
  			V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
  			
  		naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + c ); // Alemdar's Scheme
  		
  		// begin loop for element convregence
  		for ( iElementIter = 0; iElementIter < maxElementIter; iElementIter++ ) {
  		
  			// loop over integration points
  			for ( i = 0; i < numSections; i++ ) {
  				
  				sectionForceShapeFcn[i] = nd1[i] * naturalForce;
  			
  				// begin loop for section convergence
  				for ( iSectionIter = 0; iSectionIter < maxSectionIter; iSectionIter++ ) {
  					
  					// compute new sectionDefFibers
  					sectionDefFibers[i] = sectionDefFibers[i] +  sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );
  		
  					// send to section
  					for( m = 0; m < NUM_SECTION_DISPLACEMENTS; m++ )
  						sectionDefWithTorsion(m) = sectionDefFibers[i](m);
  					sectionDefWithTorsion(3) = 0.0; // set torsional strain to zero
	     					
  					int res = sections[i]->setTrialSectionDeformation(sectionDefWithTorsion);
  					
  					sectionForceFibersWithTorsion = sections[i]->getStressResultant();
  					for( m = 0; m < NUM_SECTION_DISPLACEMENTS; m++ )
  						sectionForceFibers[i](m) = sectionForceFibersWithTorsion(m);
  					
  					ksa[i] = sections[i]->getSectionTangent();
  					GJ = ksa[i](NUM_SECTION_DISPLACEMENTS,NUM_SECTION_DISPLACEMENTS);			
  					invertMatrix(NUM_SECTION_DISPLACEMENTS_WITH_TORSION,ksa[i],fsa[i]);
  					// fs and ks are the section tangent and flexibility matricies without the GJ term (size = NUM_SECTION_DISPLACEMENTS)
  					for( j = 0; j < NUM_SECTION_DISPLACEMENTS; j++ )
  						for( k = 0; k< NUM_SECTION_DISPLACEMENTS; k++ )
  							sectionFlexibility[i](j,k) = fsa[i](j,k);

  					// check if sectionForceShapeFcn ~ sectionForceFibers, break if so
  					sectionForceResidual = sectionForceShapeFcn[i] - sectionForceFibers[i];
  					
  					if ( sectionForceResidual.Norm() <= tolSectionForce ) 
  						break;  					
  					
  				}// end loop for section convergence
  			
  			} // end loop over integration points
  		
  			// compute V
  	  		V.Zero();
  	  		for( i = 0; i < numSections; i++ )
  	  			V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] );
  	  			//V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
  			
  	  		// check if V ~ 0, break if so
  	  		if ( V.Norm() <= tolV ) 
  	  			break;
  	  		
  			// compute H
  	  		Matrix H(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  			H.Zero();
  			for( i = 0; i < numSections; i++ )
  				H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
  			invertMatrix(NUM_NATRUAL_DISPLACEMENTS, H, Hinv);
  			
  			// compute new naturalForce
  			naturalForce = naturalForce + Hinv * V;
  			
  		}// end loop for element convergence
  		
  	}

  	  	
  	// Compute the following matricies: G, G2, V, c, V2, H, H12, H22, Md, Kg
  	Matrix G(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix G2(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix H(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix H12(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix H22(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix Md(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix Kg(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Vector V2(NUM_NATRUAL_DISPLACEMENTS);
  	
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
       		temp_Md(0,2) = temp_A * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
       		temp_Md(0,3) = temp_B * ( sectionDefShapeFcn[i](1) - sectionDefFibers[i](1) );
       		temp_Md(0,4) = temp_B * ( sectionDefShapeFcn[i](2) - sectionDefFibers[i](2) );
     	Md  = Md  + initialLength * wt[i] * temp_Md;
  	}
  	
 	// Compute the inverse of the H matrix
  	invertMatrix(NUM_NATRUAL_DISPLACEMENTS, H, Hinv);

  	// Compute the GMH matrix ( G + Md - H12 ) and its transpose
  	GMH = G + Md - H12;
  	//GMH = G; // Omit P-small delta
  	
  	// Compute the transposes of the following matricies: G, G2, GMH
  	Matrix GT(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix G2T(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	Matrix GMHT(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
  	for( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ ){
     	for( j = 0; j < NUM_NATRUAL_DISPLACEMENTS; j++ ){
       		GT(i,j) = G(j,i);
      		G2T(i,j) = G2(j,i);
      		GMHT(i,j) = GMH(j,i);
    	}
  	}
  	
 	
  	// Define the internal force
  	Vector internalForce(NUM_NATRUAL_DISPLACEMENTS);
  	internalForce.Zero();
  	if ( algorithmType == NUKALA_LL_ALGORITHM ) {
  		
  		// Compute new internal force
  		//internalForce = GT * naturalForce + V2 + GMHT * Hinv * c; // Alemdar's Scheme
  		internalForce = GT * naturalForce + V2 + GMHT * Hinv * V; // Nukala's Scheme
  		//internalForce = GT * naturalForce + GMHT * Hinv * c; // Omit P-small delta
  	
  	} else if ( algorithmType == NUKALA_NN_ALGORITHM ) {
  	
  		internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  		
  	}
    
    // Compute internal force for OpenSees ( i.e., add torsion and rearrange )
    for( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ )
    	internalForceOpenSees(i) = internalForce(i);
    	
    internalForceOpenSees(5) = twist*GJ/currentLength; // Add in torsional force
    internalForceOpenSees = transformNaturalCoordsT*internalForceOpenSees;
  	
    // Compute the stiffness matrix without the torsion term
    K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
    //K_temp = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta
    
    // Add in the torsional stiffness term
    kv.Zero();
    for( i = 0; i < NUM_NATRUAL_DISPLACEMENTS; i++ )
    	for( j = 0; j< NUM_NATRUAL_DISPLACEMENTS; j++ )
    		kv(i,j) = K_temp(i,j);
    
  	kv(5,5) =  GJ/currentLength; // Torsional Stiffness GJ/L

  	return 0;
}

const Matrix &
mixedBeamColumn3d::getMass(void)
{
  theMatrix.Zero();

  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
	theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*initialLength*rho;

  return theMatrix;
}


const Vector &
mixedBeamColumn3d::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

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

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	theVector += this->getRayleighDampingForces();

  } else {
	  
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	theVector += this->getRayleighDampingForces();
  }

  return theVector;    
}


void
mixedBeamColumn3d::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, mixedBeamColumn3d &E)
{
    E.Print(s);
    return s;
}


Response*
mixedBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
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

    }  else if (strcmp(argv[0],"section") ==0) { 
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
mixedBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) { // global forces
	  return eleInfo.setVector(this->getResistingForce());
	  
  } else {
    return -1;
  
  }
}

Vector
mixedBeamColumn3d::getd_hat(int sec, const Vector &v, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

   Vector D_hat(NUM_SECTION_DISPLACEMENTS);
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

Matrix
mixedBeamColumn3d::getKg(int sec, double P, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   Matrix kg(NUM_NATRUAL_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
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


Matrix
mixedBeamColumn3d::getNld_hat(int sec, const Vector &v, double L){
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
 
   Matrix Nld_hat(NUM_SECTION_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
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

Matrix
mixedBeamColumn3d::getNd2(int sec, double P, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   Matrix Nd2(NUM_SECTION_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
   Nd2.Zero();

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
   temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

   Nd2(1,1) = P * temp_A;
   Nd2(1,3) = P * temp_B;
   Nd2(2,2) = P * temp_A;
   Nd2(2,4) = P * temp_B;

   return Nd2;
}


Matrix
mixedBeamColumn3d::getNd1(int sec, const Vector &v, double L){
   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[1]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[3];

   temp_B = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[4];

   Matrix Nd1(NUM_SECTION_DISPLACEMENTS,NUM_NATRUAL_DISPLACEMENTS);
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
	
int
mixedBeamColumn3d::sendSelf(int commitTag, Channel &theChannel){  
  // @todo write mixedBeamColumn3d::sendSelf
  opserr << "Error: mixedBeamColumn3d::sendSelf -- not yet implemented for mixedBeamColumn3d element";
  return -1;
}
    
int
mixedBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker){
  // @todo write mixedBeamColumn3d::recvSelf
  opserr << "Error: mixedBeamColumn3d::sendSelf -- not yet implemented for mixedBeamColumn3d element";
  return -1;
}
