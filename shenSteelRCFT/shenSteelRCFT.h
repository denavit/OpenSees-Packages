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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/shenSteelRCFT/shenSteelRCFT.h,v $
                                                                        
#ifndef shenSteelRCFT_h
#define shenSteelRCFT_h

// Written: fmk 
//
// Description: This file contains the class definition for 
// shenSteelRCFT. shenSteelRCFT provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) shenSteelRCFT.h, revA"

#include <UniaxialMaterial.h>

#define MAT_TAG_shenSteelRCFT 58640

class shenSteelRCFT : public UniaxialMaterial
{
  public:
	shenSteelRCFT(int tag, double i1, double i2, double i3, double i4, double i5,
	  		double i6, double i7, double i8, double i9, double i10, 
	  		double i11, double i12, double i13, double i14, double i15, 
	  		double i16, double i17, double i18);
	shenSteelRCFT();

    ~shenSteelRCFT();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
	  // Input Variables
	  double fy, fu, Ee; // steel yield strength, steel ultimate strength, steel elastic modulus
	  double R; // tube slenderness
	  double epo; // initial plastic strain
	  double localBucklingStrain; // strain at initiation of local buckling
	  double Ksft; // softening slope of steel 
	  double frs; // local buckling residual ratio
	  double Rbso, Epoi, alfa, a, bb, c, w, ksi, e, fE; //parameters from the Japanese material model

	  // Computed Material Properties
	  double Rlso; //initial size of the loading surface
  
	  // State Variables
	  double trialTangent, committedTangent; //current tangent modulus of the fiber (both elastic and plastic)
	  double trialStress, committedStress; //current stress in the fiber
	  double trialStrain, committedStrain; //current total strain in the fiber
	  
	  // Flags to determine the general state of the material
	  int plasticityStatus, commitedPlasticityStatus; // flag for the status of the plasticity: 0 = Elastic, 1 = Tensile Plastic, 2 = Compressive Plastic
	  int localBucklingStatus, commitedLocalBucklingStatus; //flag if local buckling is occurring: 0=NOT occurring, 1=IS occurring and in descending branch, 2=IS occurring and in constant branch
	  int localBucklingHistory, commitedLocalBucklingHistory; //flag for local buckling history: 0=not local buckling yet, 1=some lb, still in the descending branch, 2=some lb, in the constant branch
	  int loadingDirection, committedLoadingDirection; //flag for direction of loading: 0=not set yet 1=positive 2=negative 
	  int lastYieldedIn, commitedLastYieldedIn; //flag for the the type of loading the last (or current) yielding was in 0 = no yielding yet, 1 = tension, 2 = compression
	 
	  // Updated at "end" of each step
	  double Tls_p, Cls_p, Tls_n, Cls_n; //current positive and negative values of the loading surface
	  double Tbs_p, Cbs_p, Tbs_n, Cbs_n; //current positive and negative values of the bounding surface
	  double Tmem_p, Cmem_p, Tmem_n, Cmem_n; //current positive and negative values of the memory surface
	  double Tvbs_p, Cvbs_p, Tvbs_n, Cvbs_n; //current positive and negative values of the virtual bounding surface
	  double ep, Cep; //plastic strain
	  double epmin, Cepmin, epmax, Cepmax; //minimum and maximum plastic strains experienced by the fiber
	  double W, CW; //accumulated plastic work
	  double Cmax_strs, Tmax_strs; //maximum stress attained by the material (for size of memory surface)
	  
	  // Updated upon transition from elastic to plastic
	  double delta_in, Cdelta_in; //value of delta at the initial yield state in the current loading path
	  
	  // Updated upon unloading from plasticity
	  double delta_yForT, Cdelta_yForT; //distance between bounding line and virtual bounding line used in tension
	  double delta_yForC, Cdelta_yForC; //distance between bounding line and virtual bounding line used in compression
	  double localBucklingReferenceStrain, committedLocalBucklingReferenceStrain;
	  double localBucklingCyclicReduction, committedLocalBucklingCyclicReduction;
	  
	  // local buckling state variables
	  double localBucklingConstantResidualStrain, committedLocalBucklingConstantResidualStrain; //strain at initiation of constant residual strength
  	  double localBucklingReferenceWork, committedLocalBucklingReferenceWork; // the accumulated plastic work a the initiation of the constant residual strength branch
  	  double localBucklingBaseConstantResidualStress, committedLocalBucklingBaseConstantResidualStress; // base constant residual local buckling stress (before cyclic reduction)
	  double localBucklingConstantResidualStress, committedLocalBucklingConstantResidualStress; //constant residual local buckling stress
	  double localBucklingBoundingStress, committedLocalBucklingBoundingStress; //local buckling bounding stress
};


#endif



