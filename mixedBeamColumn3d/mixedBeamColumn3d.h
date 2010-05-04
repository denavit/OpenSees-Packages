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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/mixedBeamColumn3d/mixedBeamColumn3d.h,v $
                                                                        
#ifndef mixedBeamColumn3d_h
#define mixedBeamColumn3d_h

// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign 
//
// Description: This file contains the interface for the mixedBeamColumn3d class.
// It defines the class interface and the class attributes. 
//
// What: "@(#) mixedBeamColumn3d.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

//#include <BeamIntegration.h>
//#include <SectionForceDeformation.h>
#include <FiberSectionGJ.h> // @todo Make the element work for e 
#include <CrdTransf3d.h>

class Node;
class Channel;
class Response;
class BeamIntegration;
class SectionForceDeformation;
//class FiberSectionGJ
//class CrdTransf3d

#define ELE_TAG_mixedBeamColumn3d 30765

class mixedBeamColumn3d : public Element
{
  public:
	  // constructors 
    mixedBeamColumn3d (int tag, int nodeI, int nodeJ,
		    int numSections, FiberSectionGJ *sectionPtrs[], BeamIntegration &bi,
		    CrdTransf3d &coordTransf, double massDensPerUnitLength);
    mixedBeamColumn3d ();

    // destructor
    ~mixedBeamColumn3d();

    // public methods to obtain inforrmation about dof & connectivity
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);
    
    // public methods for output    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
    friend OPS_Stream &operator<<(OPS_Stream &s, mixedBeamColumn3d &E);
    
    Response* setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInfo);    
    
  protected:
  
  private:
	// private member functions - only available to objects of the class
	Matrix getNld_hat(int sec, const Vector &v, double L);
	Vector getd_hat(int sec, const Vector &v, double L);
	Matrix getNd1(int sec, const Vector &v, double L);
	Matrix getNd2(int sec, double P, double L);
	Matrix getKg(int sec, double P, double L);  
	  
    // private attributes - a copy for each object of the class
	ID connectedExternalNodes; // tags of the end nodes
    Node *theNodes[2];   // pointers to the nodes
	BeamIntegration *beamIntegr;
    int numSections;
    FiberSectionGJ **sections;          // array of pointers to sections
    CrdTransf3d *crdTransf;        // pointer to coordinate tranformation object

    double rho;                    // mass density per unit length
    double deflength;
    double lengthLastIteration;		// the deformed length of the element in the last iteration
    double lengthLastStep;		  	// the deformed length of the element at the end of the last step
    double initialLength;

    int initialFlag;            // indicates if the element has been initialized
    int initialFlagB;  // indicates if the inital local matricies need to be computed
    int itr;
    int cnvg;
        
    Vector V;
    Vector committedV;
    Vector internalForceOpenSees;
    Vector committedInternalForceOpenSees;
    Vector naturalForce;
    Vector commitedNaturalForce;
    Vector lastNaturalDisp;
    Vector commitedLastNaturalDisp;
    Vector c;
    Vector commitedC;
    Matrix Hinv;
    Matrix commitedHinv;
    Matrix GMH;
    Matrix commitedGMH;
    Matrix kv;                     // stiffness matrix in the basic system
    Matrix kvcommit;               // commited stiffness matrix in the basic system
    
    Matrix *Ki;
    
    Vector *sectionForceFibers;
    Vector *commitedSectionForceFibers;
    Vector *sectionDefFibers;
    Vector *commitedSectionDefFibers;
    Matrix *sectionFlexibility;
    Matrix *commitedSectionFlexibility;

    // static data - single copy for all objects of the class
    static int maxNumSections;
    static Matrix theMatrix;
    static Vector theVector;
    static double workArea[];
    static Matrix transformNaturalCoords; 
    static Matrix transformNaturalCoordsT;
    	// matrix to transform the natural coordinates from what the coordinate transformation uses and what the element uses 	
    
    // These variable are always recomputed, so there is no need to store them for each instance of the element
    static Vector *sectionDefShapeFcn;
    static Matrix *ks;                    // array of section stiffness matrices
    static Matrix *fs;
    static Matrix *ksa;
    static Matrix *fsa;
    static Matrix *nldhat;
    static Matrix *nldhatT;
    static Matrix *nd1; 
    static Matrix *nd2;
    static Matrix *nd1T;
    static Matrix *nd2T; 
};

#endif

