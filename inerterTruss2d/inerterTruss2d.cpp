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

#include "inerterTruss2d.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix inerterTruss2d::inerterTruss2dM4(4,4);
Matrix inerterTruss2d::inerterTruss2dM6(6,6);
Vector inerterTruss2d::inerterTruss2dV4(4);
Vector inerterTruss2d::inerterTruss2dV6(6);


#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

// Documentation
// Truss Element with Inertance
//
// element inerterTruss2d $tag $iNode $jNode $k $c $b <-rho $rho> <-consistentMass $cFlag> <-doRayleigh $rFlag> <-geomNonlinear $gFlag>
//
// Required Input Parameters:
//   $tag                   integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $k                     stiffness (units of force per length)
//   $c                     damping   (units of force per velocity)
//   $b                     inertance (units of force per acceleration)
//
// Optional Input:
//   $rho              element mass density (units of mass per unit length) default = 0.0
//   $cFlag                 consistent mass flag
//                              cFlag = 0, use lumped mass matrix (default)
//                              cFlag = 1, use consistent mass matrix
//   $rFlag                 Rayleigh damping flag
//                              rFlag = 0, no Rayleigh damping (default)
//                              rFlag = 1, include Rayleigh damping
//   $gFlag                 geometric nonlinearity flag
//                              gFlag = 0, use geometric linear formulation (default)
//                              gFlag = 1, use geometric nonlinear formulation
//
// Element Notes:
//   To avoid conflict with mass matrix used to develop forces in the ground motion,
//   inertance is defined with a restoring force. Accordingly, a nonlinear solver is
//   necessary when using this element, even for linear problems.
//

OPS_Export void localInit() {
    OPS_Error("inerterTruss2d element\n  "
       "Developed by Mark D. Denavit, Abdollah Javidialesaadi, and Nicholas E. Wierschem\n  "
       "University of Tennessee, Knoxville\n", 1);
}

OPS_Export void * OPS_inerterTruss2d() {

    Element *theTruss = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingArgs == 0) { // parallel processing
        theTruss = new inerterTruss2d();
        return theTruss;
    }

    if (numRemainingArgs < 6) {
        opserr << "ERROR - inerterTruss2d not enough args provided, want: element inerterTruss2d tag? iNode? jNode? k? c? b?\n";
    }

    // get the id and end nodes
    int iData[3];
    double dData[3];
    char *sData  = new char[40];
    int numData;
    int ndm = OPS_GetNDM();

    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data\n";
        return 0;
    }
    
    int eleTag = iData[0];
    
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element area for element" << eleTag << endln;
        return 0;
    }    

    // Set Default Values for Optional Input
    double rho = 0.0;
    int cFlag = 0;
    int rFlag = 0;
    int gFlag = 0;

    // Loop through remaining arguments to get optional input
    while ( OPS_GetNumRemainingInputArgs() > 0 ) {
        if ( OPS_GetStringCopy(&sData) != 0 ) {
            opserr << "WARNING invalid input";
            return 0;
        }

        if ( strcmp(sData,"-rho") == 0 ) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &rho) != 0) {
                opserr << "WARNING invalid input, want: -rho $rho \n";
                return 0;
            }

        } else if ( (strcmp(sData,"-consistentMass") == 0) || (strcmp(sData,"-cMass") == 0) ) {
            numData = 1;
            if (OPS_GetInt(&numData, &cFlag) != 0) {
                opserr << "WARNING: Invalid input, want: -consistentMass $cFlag in element inerterTruss2d " << eleTag;
                return 0;
            }

        } else if ( strcmp(sData,"-doRayleigh") == 0 ) {
            numData = 1;
            if (OPS_GetInt(&numData, &rFlag) != 0) {
                opserr << "WARNING: Invalid input, want: -doRayleigh $rFlag in element inerterTruss2d " << eleTag;
                return 0;
            }

        } else if ( strcmp(sData,"-geomNonlinear") == 0 ) {
            numData = 1;
            if (OPS_GetInt(&numData, &gFlag) != 0) {
                opserr << "WARNING: Invalid input, want: -geomNonlinear $gFlag in element inerterTruss2d " << eleTag;
                return 0;
            }

        } else {
            opserr << "WARNING unknown option " << sData << "\n";
        }
    }    

    
    // now create the truss and add it to the Domain
    theTruss = new inerterTruss2d(eleTag, ndm, iData[1], iData[2], dData[0], dData[1], dData[2], 
        rho, cFlag, rFlag, gFlag);

    if (theTruss == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }

    return theTruss;
}

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the inerterTruss2d end nodes.

inerterTruss2d::inerterTruss2d(int tag, int dim,
           int Nd1, int Nd2, double k1, double c1, double b1,
           double rho1, int cFlag1, int rFlag1, int gFlag1):
    Element(tag, 0), externalNodes(2), 
    numDIM(dim), numDOF(0),
    theMatrix(0), theVector(0),
    k(k1), c(c1), b(b1), Lo(0.0), Lc(0.0), so(0.0), sc(0.0), co(0.0), cc(0.0),
    rho(rho1), cFlag(cFlag1), rFlag(rFlag1), gFlag(gFlag1)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (externalNodes.Size() != 2) {
      opserr << "FATAL inerterTruss2d::inerterTruss2d - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
inerterTruss2d::inerterTruss2d():
    Element(0, 0), externalNodes(2), 
    numDIM(0), numDOF(0),
    theMatrix(0), theVector(0),
    k(0.0), c(0.0), b(0.0), Lo(0.0), Lc(0.0), so(0.0), sc(0.0), co(0.0), cc(0.0),
    rho(0.0), cFlag(0), rFlag(0), gFlag(0)
{
    // ensure the connectedExternalNode ID is of correct size 
    if (externalNodes.Size() != 2) {
        opserr << "FATAL inerterTruss2d::inerterTruss2d - failed to create an ID of size 2\n";
        exit(-1);
    }

    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
inerterTruss2d::~inerterTruss2d() {
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
}


int
inerterTruss2d::getNumExternalNodes(void) const {
    return 2;
}

const ID &
inerterTruss2d::getExternalNodes(void) {
    return externalNodes;
}

Node **
inerterTruss2d::getNodePtrs(void) {
  return theNodes;
}

int
inerterTruss2d::getNumDOF(void) {
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the inerterTruss2d element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
inerterTruss2d::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        Lo = 0;
        Lc = 0;
        return;
    }

    // first set the node pointers
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
        if (theNodes[0] == 0)
            opserr <<"inerterTruss2d::setDomain() - inerterTruss2d" << this->getTag() << " node " << Nd1 <<
                "does not exist in the model\n";
        else
            opserr <<"inerterTruss2d::setDomain() - inerterTruss2d" << this->getTag() << " node " << Nd2 <<
                "does not exist in the model\n";

        // fill this in so don't segment fault later
        numDOF = 4;
        theMatrix = &inerterTruss2dM4;
        theVector = &inerterTruss2dV4;	

        return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
        opserr <<"WARNING inerterTruss2d::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
            "have differing dof at ends for inerterTruss2d " << this->getTag() << endln;

        // fill this in so don't segment fault later
        numDOF = 4;    
        theMatrix = &inerterTruss2dM4;
        theVector = &inerterTruss2dV4;	
	
        return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 2 && dofNd1 == 2) {
        numDOF = 4;
        theMatrix = &inerterTruss2dM4;
        theVector = &inerterTruss2dV4;	
    } else if (numDIM == 2 && dofNd1 == 3) {
        numDOF = 6;	
        theMatrix = &inerterTruss2dM6;
        theVector = &inerterTruss2dV6;		
    } else {
        opserr <<"WARNING inerterTruss2d::setDomain cannot handle " << numDIM << " dofs at nodes in " << 
            dofNd1  << " problem\n";

        numDOF = 4;
        theMatrix = &inerterTruss2dM4;
        theVector = &inerterTruss2dV4;	
        return;
    }


    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    // Set undeformed and initial length
    double Lx = end2Crd(0)-end1Crd(0);
    double Ly = end2Crd(1)-end1Crd(1);
    Lo = sqrt(Lx*Lx + Ly*Ly);

    if (Lo == 0.0) {
        opserr <<"WARNING inerterTruss2d::setDomain() - inerterTruss2d " << this->getTag() << " has zero length\n";
        return;
    }

    co = Lx/Lo;
    so = Ly/Lo;

    // Set current to initial
    Lc = Lo;
    cc = co;
    sc = so;

    return;
}   	 


int
inerterTruss2d::commitState() {
    return this->Element::commitState();
}

int
inerterTruss2d::revertToLastCommit() {
    return 0;
}

int
inerterTruss2d::revertToStart() {
    return 0;
}

int
inerterTruss2d::update(void) {
    // Nodal displacements, velocities, and accelerations	
    const Vector &end1Disp  = theNodes[0]->getTrialDisp();
    const Vector &end2Disp  = theNodes[1]->getTrialDisp();    
    const Vector &end1Vel   = theNodes[0]->getTrialVel();
    const Vector &end2Vel   = theNodes[1]->getTrialVel();	
    const Vector &end1Accel = theNodes[0]->getTrialAccel();
    const Vector &end2Accel = theNodes[1]->getTrialAccel();


    if (gFlag == 0) { // geometric linear     
        dc =   (end2Disp(0)-end1Disp(0))*co +   (end2Disp(1)-end1Disp(1))*so;
        vc =     (end2Vel(0)-end1Vel(0))*co +     (end2Vel(1)-end1Vel(1))*so;
        ac = (end2Accel(0)-end1Accel(0))*co + (end2Accel(1)-end1Accel(1))*so; 
 
    } else { // geometric nonlinear 
        const Vector &end1Crd   = theNodes[0]->getCrds();
        const Vector &end2Crd   = theNodes[1]->getCrds();
        double Lx = (end2Crd(0)+end2Disp(0))-(end1Crd(0)+end1Disp(0));
        double Ly = (end2Crd(1)+end2Disp(1))-(end1Crd(1)+end1Disp(1));
        Lc = sqrt(Lx*Lx+Ly*Ly);
        cc = Lx/Lc;
        sc = Ly/Lc;
        dc = Lc-Lo;
        vc =     (end2Vel(0)-end1Vel(0))*cc +     (end2Vel(1)-end1Vel(1))*sc;
        ac = (end2Accel(0)-end1Accel(0))*cc + (end2Accel(1)-end1Accel(1))*sc; 
    }

    return 0;
}


const Matrix &
inerterTruss2d::getTangentStiff(void) {
    Matrix &K = *theMatrix;
    K.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return K;
    }

    // transformation matrix
    static Matrix T(4,4);
    T.Zero();
    double cosq, sinq;
    if (gFlag == 0) { 
        cosq = co; // geometric linear  
        sinq = so;
    } else { 
        cosq = cc; // geometric nonlinear  
        sinq = sc;
    }
    T(0,0) =  cosq;
    T(1,0) = -sinq;
    T(0,1) =  sinq;
    T(1,1) =  cosq;
    T(2,2) =  cosq;
    T(3,2) = -sinq;
    T(2,3) =  sinq;
    T(3,3) =  cosq;

    static Matrix kb(4,4);
    kb.Zero();
    kb(0,0) =  k;
    kb(2,0) = -k;
    kb(0,2) = -k;
    kb(2,2) =  k;
    if (gFlag == 1) { 
        double kg1 = (k*dc + c*vc + b*ac)/Lc;
        kb(1,1) =  kg1;
        kb(3,1) = -kg1;
        kb(1,3) = -kg1;
        kb(3,3) =  kg1;
    }

    static Matrix ke(4,4);
    ke.addMatrixTripleProduct(0.0, T, kb, 1.0);

    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i,j)                 += ke(i,j);
            K(i+numDOF2,j)         += ke(i+numDIM,j);
            K(i,j+numDOF2)         += ke(i,j+numDIM);
            K(i+numDOF2,j+numDOF2) += ke(i+numDIM,j+numDIM);
        }
    }

    // opserr << K(0,0) << " " << K(0,1) << " " << K(0,2) << " " << K(0,3) << " " << K(0,4) << " " << K(0,5) << "\n";
    // opserr << K(1,0) << " " << K(1,1) << " " << K(1,2) << " " << K(1,3) << " " << K(1,4) << " " << K(1,5) << "\n";
    // opserr << K(2,0) << " " << K(2,1) << " " << K(2,2) << " " << K(2,3) << " " << K(2,4) << " " << K(2,5) << "\n";
    // opserr << K(3,0) << " " << K(3,1) << " " << K(3,2) << " " << K(3,3) << " " << K(3,4) << " " << K(3,5) << "\n";
    // opserr << K(4,0) << " " << K(4,1) << " " << K(4,2) << " " << K(4,3) << " " << K(4,4) << " " << K(4,5) << "\n";
    // opserr << K(5,0) << " " << K(5,1) << " " << K(5,2) << " " << K(5,3) << " " << K(5,4) << " " << K(5,5) << "\n\n";

    return K;
}


const Matrix &
inerterTruss2d::getInitialStiff(void) {
    Matrix &K = *theMatrix;
    K.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return K;
    }
    
    // transformation matrix
    static Matrix T(4,4);
    T.Zero();
    T(0,0) =  co;
    T(1,0) =  so;
    T(0,1) = -so;
    T(1,1) =  co;
    T(2,2) =  co;
    T(3,2) =  so;
    T(2,3) = -so;
    T(3,3) =  co;

    static Matrix kb(4,4);
    kb.Zero();
    kb(0,0) =  k;
    kb(2,0) = -k;
    kb(0,2) = -k;
    kb(2,2) =  k;

    static Matrix ke(4,4);
    ke.addMatrixTripleProduct(0.0, T, kb, 1.0);

    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i,j)                 += ke(i,j);
            K(i+numDOF2,j)         += ke(i+numDIM,j);
            K(i,j+numDOF2)         += ke(i,j+numDIM);
            K(i+numDOF2,j+numDOF2) += ke(i+numDIM,j+numDIM);
        }
    }

    return K;
}

const Matrix &
inerterTruss2d::getDamp(void) {
    Matrix &C = *theMatrix;
    C.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return C;
    }

    // get Rayleigh damping
    if (rFlag == 1)
        C = this->Element::getDamp();

    // add in defined damping
    if (c != 0.0) {
        // transformation matrix
        static Matrix T(4,4);
        T.Zero();
        double cosq, sinq;
        if (gFlag == 0) { 
            cosq = co; // geometric linear  
            sinq = so;
        } else { 
            cosq = cc; // geometric nonlinear  
            sinq = sc;
        }
        T(0,0) =  cosq;
        T(1,0) = -sinq;
        T(0,1) =  sinq;
        T(1,1) =  cosq;
        T(2,2) =  cosq;
        T(3,2) = -sinq;
        T(2,3) =  sinq;
        T(3,3) =  cosq;

        static Matrix cb(4,4);
        cb.Zero();
        cb(0,0) =  c;
        cb(0,2) = -c;
        cb(2,0) = -c;
        cb(2,2) =  c;

        static Matrix ce(4,4);
        ce.Zero();
        ce.addMatrixTripleProduct(0.0, T, cb, 1.0);

        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            for (int j = 0; j < numDIM; j++) {
                C(i,j)                 += ce(i,j);
                C(i+numDOF2,j)         += ce(i+numDIM,j);
                C(i,j+numDOF2)         += ce(i,j+numDIM);
                C(i+numDOF2,j+numDOF2) += ce(i+numDIM,j+numDIM);
            }
        }
    }

    //opserr << C(0,0) << " " << C(0,1) << " " << C(0,2) << " " << C(0,3) << "\n";
    //opserr << C(1,0) << " " << C(1,1) << " " << C(1,2) << " " << C(1,3) << "\n";
    //opserr << C(2,0) << " " << C(2,1) << " " << C(2,2) << " " << C(2,3) << "\n";
    //opserr << C(3,0) << " " << C(3,1) << " " << C(3,2) << " " << C(3,3) << "\n";

    return C;
}


const Matrix &
inerterTruss2d::getMass(void) {
    // zero the matrix
    Matrix &M = *theMatrix;
    M.Zero();
  
    // check for quick return
    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return M;
    }
  
    if (rho != 0.0) {
        if (cFlag == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*Lo;
            int numDOF2 = numDOF/2;
            for (int i = 0; i < numDIM; i++) {
                M(i,i) = m;
                M(i+numDOF2,i+numDOF2) = m;
            }
        } else {
            // consistent mass matrix
            double m = rho*Lo/6.0;
            int numDOF2 = numDOF/2;
            for (int i = 0; i < numDIM; i++) {
                M(i,i) = 2.0*m;
                M(i,i+numDOF2) = m;
                M(i+numDOF2,i) = m;
                M(i+numDOF2,i+numDOF2) = 2.0*m;
            }
        }
    }    


    if (b != 0.0) {
        // transformation matrix
        static Matrix T(4,4);
        T.Zero();
        double cosq, sinq;
        if (gFlag == 0) { 
            cosq = co; // geometric linear  
            sinq = so;
        } else { 
            cosq = cc; // geometric nonlinear  
            sinq = sc;
        }
        T(0,0) =  cosq;
        T(1,0) = -sinq;
        T(0,1) =  sinq;
        T(1,1) =  cosq;
        T(2,2) =  cosq;
        T(3,2) = -sinq;
        T(2,3) =  sinq;
        T(3,3) =  cosq;

        static Matrix bb(4,4);
        bb.Zero();
        bb(0,0) =  b;
        bb(0,2) = -b;
        bb(2,0) = -b;
        bb(2,2) =  b;

        static Matrix be(4,4);
        be.addMatrixTripleProduct(0.0, T, bb, 1.0);

        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            for (int j = 0; j < numDIM; j++) {
                M(i,j)                 += be(i,j);
                M(i+numDOF2,j)         += be(i+numDIM,j);
                M(i,j+numDOF2)         += be(i,j+numDIM);
                M(i+numDOF2,j+numDOF2) += be(i+numDIM,j+numDIM);
            }
        }
    }

    return M;
}

/*
void 
inerterTruss2d::zeroLoad(void) {
    // nothing to do here
}

int 
inerterTruss2d::addLoad(ElementalLoad *theLoad, double loadFactor) {  
    opserr << "Error: inerterTruss2d::addLoad -- not yet implemented";
    return -1;
}

int 
inerterTruss2d::addInertiaLoadToUnbalance(const Vector &accel) {
    opserr << "Error: inerterTruss2d::addInertiaLoadToUnbalance -- not yet implemented";
    return -1;
}
*/

const Vector &
inerterTruss2d::getResistingForce() {	
    Vector &P = *theVector;
    P.Zero();

    if (Lo == 0.0) { // - problem in setDomain() no further warnings
        return P;
    }
    
    // R = Ku - Pext
    // Ku = F * transformation
    double force = k*dc + c*vc + b*ac;

    double cosq, sinq;
    if (gFlag == 0) { 
        cosq = co; // geometric linear  
        sinq = so;
    } else { 
        cosq = cc; // geometric nonlinear  
        sinq = sc;
    }

    int numDOF2 = numDOF/2;
    P(0)         = -force*cosq;
    P(1)         = -force*sinq;
    P(numDOF2+0) =  force*cosq;
    P(numDOF2+1) =  force*sinq;


    // opserr << P(0) << "\n";
    // opserr << P(1) << "\n";
    // opserr << P(2) << "\n";
    // opserr << P(3) << "\n";
    // opserr << P(4) << "\n";
    // opserr << P(5) << "\n\n";

    return P;
}


const Vector &
inerterTruss2d::getResistingForceIncInertia() {	

    Vector &P = *theVector;
    P = this->getResistingForce();

    // now include the mass portion
    if (Lo != 0.0 && rho != 0.0) {
    
        // add inertia forces from element mass
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();	
        
        int numDOF2 = numDOF/2;
        
        if (cFlag == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*Lo;
            for (int i = 0; i < numDIM; i++) {
                P(i) += m*accel1(i);
                P(i+numDOF2) += m*accel2(i);
            }
        } else  {
            // consistent mass matrix
            double m = rho*Lo/6.0;
            for (int i=0; i<numDIM; i++) {
                P(i) += 2.0*m*accel1(i) + m*accel2(i);
                P(i+numDOF2) += m*accel1(i) + 2.0*m*accel2(i);
            }
        }
    }    
    
    // add the damping forces if rayleigh damping
    if (rFlag == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
        theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);

        //        P += this->getRayleighDampingForces();
      
    return P;
}

int
inerterTruss2d::sendSelf(int commitTag, Channel &theChannel) {
    opserr << "Error: inerterTruss2d::sendSelf -- not yet implemented\n";
    return -1;  
}

int
inerterTruss2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    opserr << "Error: inerterTruss2d::recvSelf -- not yet implemented\n";
    return -1;
}
  
int
inerterTruss2d::displaySelf(Renderer &theViewer, int displayMode, float fact, 
		   const char **displayModes, int numModes) {
    opserr << "Error: inerterTruss2d::displaySelf -- not yet implemented\n";
    return -1;
}

void
inerterTruss2d::Print(OPS_Stream &s, int flag) {
    s << "Element: " << this->getTag() << "\n";
    s << " type: inerterTruss2d\n";
    s << " iNode: " << externalNodes(0) << "\n";
    s << " jNode: " << externalNodes(1) << "\n";
    s << " k: " << k << "\n";
    s << " c: " << c << "\n";
    s << " b: " << b << "\n";
}

Response*
inerterTruss2d::setResponse(const char **argv, int argc, OPS_Stream &output) {
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","inerterTruss2d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",externalNodes[0]);
    output.attr("node2",externalNodes[1]);

    //
    // we compare argv[0] for known response types for the inerterTruss2d
    //

    if (strcmp(argv[0],"force") == 0 || 
        strcmp(argv[0],"forces") == 0 || 
        strcmp(argv[0],"globalForce") == 0 || 
        strcmp(argv[0],"globalForces") == 0) {
            
        char outputData[10];
        int numDOFperNode = numDOF/2;
        for (int i=0; i<numDOFperNode; i++) {
            sprintf(outputData,"P1_%d", i+1);
            output.tag("ResponseType", outputData);
        }
        for (int j=0; j<numDOFperNode; j++) {
            sprintf(outputData,"P2_%d", j+1);
            output.tag("ResponseType", outputData);
        }
        theResponse =  new ElementResponse(this, 1, Vector(numDOF));

    } else if (strcmp(argv[0],"axialForce") == 0 || 
               strcmp(argv[0],"basicForce") == 0 || 
               strcmp(argv[0],"localForce") == 0 || 
               strcmp(argv[0],"basicForces") == 0) {
        
        output.tag("ResponseType", "N");
        theResponse =  new ElementResponse(this, 2, Vector(1));

    } else if (strcmp(argv[0],"defo") == 0 || 
               strcmp(argv[0],"deformation") == 0 ||
               strcmp(argv[0],"deformations") == 0 || 
               strcmp(argv[0],"basicDefo") == 0 ||
               strcmp(argv[0],"basicDeformation") == 0 || 
               strcmp(argv[0],"basicDeformations") == 0) {

        output.tag("ResponseType", "U");
        theResponse = new ElementResponse(this, 3, Vector(1));

    } 

    output.endTag();

    return theResponse;
}

int 
inerterTruss2d::getResponse(int responseID, Information &eleInfo) {
    static Vector fVec(1);
    static Matrix kVec(1,1);

    switch (responseID) {
        case 1:
            return eleInfo.setVector(this->getResistingForce());

        case 2:
            fVec(0) = k*dc + c*vc + b*ac;
            return eleInfo.setVector(fVec);

        case 3:
            fVec(0) = dc;
            return eleInfo.setVector(fVec);
      
        default:
            return 0;
    }
}
