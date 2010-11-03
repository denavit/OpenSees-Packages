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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/shenSteelCCFT/shenSteelCCFT.h,v $
                                                                        
#ifndef shenSteelCCFT_h
#define shenSteelCCFT_h

// Written: fmk 
//
// Description: This file contains the class definition for 
// shenSteelCCFT. shenSteelCCFT provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) shenSteelCCFT.h, revA"

#include <UniaxialMaterial.h>

#define MAT_TAG_shenSteelCCFT 58639

class shenSteelCCFT : public UniaxialMaterial
{
  public:
	shenSteelCCFT(int tag, double i1, double i2, double i3, double i4, double i5, 
	  		double i6, double i7, double i8, double i9, double i10, 
	  		double i11, double i12, double i13, double i14, double i15, 
	  		double i16, double i17, double i18, double i19);  
	shenSteelCCFT();    

    ~shenSteelCCFT();

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
	  void commitStatevar(void); //function to commit all state variables en masse
	  void backtocommitStatevar(void); //function to revert all state variables to last commit en masse
	  
	  // Input Variables
	  double fy, fu, Ee; // steel yield strength, steel ultimate strength, steel elastic modulus
	  double R; // tube slenderness
	  double epo; // inital plastic strain
	  double alphaHoop; // ratio of lateral hoop stress to yield stress used to determine biaxial reductions
	  double localBucklingStrain; // strain at initiation of local buckling
	  double Ksft; // softening slope of steel 
	  double Rcrit; // critical tube slenderness for determining residual strength
	  double Rbso, Epoi, alfa, a, bb, c, w, ksi, e, fE; //parameters from the Japanese material model


	  // Material Properties
	  double Rlso; //initial size of the loading surface
	  
	  
	  // State Variables
	  
	  // Flags to determine the general state of the material
	  int elastic, Celastic; //flag for type of loading: 1=elastic 0=plastic
	  int plastic, Cplastic; //flag for type of loading: 1=plastic 0=elastic
	  int lb, Clb; //flag if local buckling is occuring: 1=IS occuring 0=NOT occuring
	  int elb, Celb; //flag for local buckling history: 0=not local buckling yet, 1=some lb, still in the decending branch, 2=some lb, in the constant branch
	  int loadingDirection, committedLoadingDirection; //flag for direction of loading: 0=not set yet 1=positive 2=negative 
	  int Trule, Crule; //flag for type of loading: 1=elastic 2=plastic
	  int lastYieldedIn, commitedLastYieldedIn; //flag for the the type of loading the last (or current) yielding was in 0 = no yielding yet, 1 = tension, 2 = compression
	 
	  // Updated at "beginning" of each step
	  double Ep, CEp; //current plastic stiffness
	  double Epo, CEpo; //slope of the current bounding line
	  double delta, Cdelta; //distance between bounding line and loading point
	  double h, Ch; //shape parameter for determination of Ep
	  double trialTangent, committedTangent; //current tangent modulus of the fiber (both elastic and plastic)
	  double trialStress, committedStress; //current stress in the fiber
	  double trialStrain, committedStrain; //current total strain in the fiber
	  
	  // Updated at "end" of each step
	  double Rls, CRls; //current size of the loading surface
	  double Tls_p, Cls_p, Tls_n, Cls_n; //current positive and negative values of the loading surface
	  double Tbs_p, Cbs_p, Tbs_n, Cbs_n; //current positive and negative values of the bounding surface
	  double Tmem_p, Cmem_p, Tmem_n, Cmem_n; //current positive and negative values of the memory surface
	  double Tvbs_p, Cvbs_p, Tvbs_n, Cvbs_n; //current positive and negative values of the virtual bounding surface
	  double ep, Cep; //plastic strain
	  double epmin, Cepmin, epmax, Cepmax; //minimum and maximum plastic strains experienced by the fiber
	  double ebar_p, Cebar_p; //range of plastic strain experienced by the fiber (epmax-epmin)
	  double W, CW; //accumulated plastic work
	  double Cmax_strs, Tmax_strs; //maximum stress attained by the material (for size of memory surface)
	  
	  // Updated upon reversal
	  double delta_y, Cdelta_y; //distance between bounding line and virtual bounding line
	  
	  // Updated upon transition from elastic to plastic
	  double delta_in, Cdelta_in; //value of delta at the inital yield state in the current loading path
	  
	  // Updated upon unloading from tensile plasticity
	  double localBucklingReferenceStrain, committedLocalBucklingReferenceStrain;
	  
	  // Updated upon entering local buckling
	  double localBucklingConstantResidualStrain, committedLocalBucklingConstantResidualStrain; //strain at initiation of constant residual strength
	  double localBucklingConstantResidualStress, committedLocalBucklingConstantResidualStress; //constant residual local buckling stress

	  // Updated each step during local buckling
	  double localBucklingBoundingStress, committedLocalBucklingBoundingStress; //local buckling bounding stress
	  
};


#endif



