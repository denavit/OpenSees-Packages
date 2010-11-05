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
// $Date: 2010/05/04 17:14:46 $
// $Source: /cvsroot/openseescomp/CompositePackages/shenSteel01/shenSteel01.h,v $

/* ****************************************************************** **
**                                                                    **
** This file was developed by:                                        **
**   Mark D. Denavit, University of Illinois at Urbana Champaign      **
**   Jerome F. Hajjar, Northeastern University                        **
**                                                                    **
********************************************************************* */

#ifndef shenSteel01_h
#define shenSteel01_h

#include <UniaxialMaterial.h>

#define MAT_TAG_shenSteel01 58641

class shenSteel01 : public UniaxialMaterial
{
  public:
	shenSteel01(int tag, double iA1, double iA2, double iA3,
			double iB1, double iB2, double iB3, double iB4, double iB5, double iB6, double iB7, double iB8, double iB9, double iB10,
	  		bool iC1, double iC2, double iC3, double iC4,
	  		double iD1, double iD2, double iD3,
	  		bool iE1, double iE2, double iE3, double iE4, int iE5,
	  		bool iF1, double iF2, double iF3,
	  		bool iG1, double iG2, double iG3,
	  		bool iH1, double iH2, double iH3);
	shenSteel01(int tag);
    ~shenSteel01();

    const char *getClassType(void) const {return "shenSteel01";};

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
	  double Rbso, Epoi, alfa, a, bb, c, w, ksi, e, fE; //parameters from the Japanese material model
	  bool modelYieldPlateau; // flag for modeling the yield plateau (true = yield plateau modeled)
	  double M, Est, est; // yield plateau parameters
	  double initStress; // initial stress (for residual stress)
	  double alphaLat; // biaxial stress ratio
	  double ep0; // initial plastic strain
	  bool modelLocalBuckling; // flag for modeling local buckling (true = local buckling modeled)
	  double localBucklingStrain, Ksft, alphaFulb; // local buckling parameters
	  int refFulb; // local buckling parameter
	  bool modelDegradeEp;  // flag for modeling degradation in Ep (true = it is modeled)
	  double degradeEpRate, degradeEpLimit; // degradation parameters
	  bool modelDegradeKappa;  // flag for modeling degradation in kappa (true = it is modeled)
	  double degradeKappaRate, degradeKappaLimit; // degradation parameter
	  bool modelDegradeFulb;  // flag for modeling degradation in Fulb (true = it is modeled)
	  double degradeFulbRate, degradeFulbLimit; // degradation parameter

	  // Computed Material Properties
	  double Rlso; //initial size of the loading surface

	  // State Variables
	  double trialTangent, committedTangent; //current tangent modulus of the fiber (both elastic and plastic)
	  double trialStress, committedStress; //current stress in the fiber
	  double trialStrain, committedStrain; //current total strain in the fiber
	  
	  // Flags to determine the general state of the material
	  int plasticityStatus, commitedPlasticityStatus; // flag for the status of the plasticity: 0 = Elastic, 1 = Tensile Plastic, 2 = Compressive Plastic
	  int yieldPlateauStatus, commitedYieldPlateauStatus; // flag for the status of the yield plateau: 0 = Vanished (not in post plateau status), 1 = Not Vanished, 2 = Vanished and in post plateau status
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



