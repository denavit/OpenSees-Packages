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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/multiSurfaceKinematicHardening/multiSurfaceKinematicHardening.h,v $

#ifndef multiSurfaceKinematicHardening_h
#define multiSurfaceKinematicHardening_h

// Written: Mark D. Denavit
//
// What: "@(#) multiSurfaceKinematicHardening.h, revA"

#include <UniaxialMaterial.h>

#define MAT_TAG_multiSurfaceKinematicHardening 58646

class multiSurfaceKinematicHardening : public UniaxialMaterial
{
  public:
    multiSurfaceKinematicHardening(int tag, double iE, int iNumSurfaces, double iCenters[], double iRadii[], double iHardeningModulii[], double iSigma0);
    multiSurfaceKinematicHardening();

    ~multiSurfaceKinematicHardening();

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
    double tangentModulus(int activeSurface);

    int numSurfaces; // Number of Surfaced
    double E; // Elastic Modulus
    double initStress; // initial stress (for residual stress)

    double trialStrain, committedStrain;
    double trialStress, committedStress;
    double trialTangent, committedTangent;

    double *trialCenters;
    double *committedCenters;
    double *startCenters;
    double *radii;
    double *hardeningModulii;
};


#endif



