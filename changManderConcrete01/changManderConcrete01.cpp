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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/changManderConcrete01/changManderConcrete01.cpp,v $

// Documentation: Chang and Mander Concrete Model
// uniaxialMaterial changManderConcrete01 $tag $fcc $ecc $Ec $rn_pre $rn_post $ft $et $rp $xp_cr <options>
//
// Required Input Parameters:
//   $tag           integer tag identifying material
//   $fcc           peak compressive stress (input as negative value)
//   $ecc           strain at peak compressive stress (input as negative value)
//   $Ec            initial modulus of elasticity
//   $rn_pre        value of factor r to be used in the compressive backbone curve before the peak stress
//   $rn_post       value of factor r to be used in the compressive backbone curve after the peak stress
//   $ft            peak tensile stress
//   $et            strain at peak tensile stress
//   $rp            value of factor r to be used in the tensile backbone curve
//   $xp_cr         critical strain ratio for the tensile backbone curve
//
// Optional Input:
//   -spall $xn_cr
//     $xn_cr       critical strain ratio for the tensile backbone curve
//
// References:
//   1. G. A Chang and J. B Mander, Seismic Energy Based Fatigue Damage Analysis of Bridge Columns:
//      Part I - Evaluation of Seismic Capacity (National Center for Earthquake Engineering Research,
//      State University of New York at Buffalo, Department of Civil Engineering, 1994).
//   2. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.

#include <elementAPI.h>
#include "changManderConcrete01.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#define SMALL_STRAIN 1.0e-20
//#define SMALL_STRESS 1.0e-16
#define SMALL_TANGENT 1.0e-7
#define SMALL_STRESS 0.0
//#define SMALL_TANGENT 0.0
#define SMALL_NUMBER 1.0e-20

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
  OPS_Error("changManderConcrete01 unaxial material \nWritten by Mark D. Denavit, University of Illinois at Urbana-Champaign\n", 1);
}

OPS_Export void *
OPS_changManderConcrete01()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[9];
  double xn_cr = 0.0;
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial changManderConcrete01 tag \n" << endln;
    return 0;
  }

  numData = 9;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input, want: uniaxialMaterial changManderConcrete01 tag fcc ecc Ec rn_pre rn_post ft et rp xp_cr \n";
    return 0;
  }

  // Check the data
  if ( dData[0] >= 0 ) {
	  opserr << "WARNING fcc should be input as a negative value \n";
	  return 0;
  }
  if ( dData[1] >= 0 ) {
	  opserr << "WARNING ecc should be input as a negative value \n";
	  return 0;
  }
  if ( dData[2] <= 0 ) {
	  opserr << "WARNING Ec should be a positive value \n";
	  return 0;
  }
  if ( dData[3] <= 0 ) {
	  opserr << "WARNING rn_pre should be a positive value \n";
	  return 0;
  }
  if ( dData[4] < 0 ) {
	  opserr << "WARNING rn_post should be a positive value (or zero if no spalling) \n";
	  return 0;
  }
  if ( dData[5] < 0 ) {
	  opserr << "WARNING ft should be input as a positive value \n";
	  return 0;
  }
  if ( dData[6] < 0 ) {
	  opserr << "WARNING et should be input as a positive value \n";
	  return 0;
  }
  if ( dData[7] <= 0 ) {
	  opserr << "WARNING rp should be a positive value \n";
	  return 0;
  }
  if ( dData[8] <= 1 ) {
	  opserr << "WARNING xp_cr should be greater than 1.0 \n";
	  return 0;
  }

  int sDataLength = 20;
  char *sData = new char[sDataLength];

  // Loop through remaining arguments
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
	  if ( OPS_GetString(sData, sDataLength) != 0 ) {
		  opserr << "WARNING invalid input \n";
		  return 0;
	  }

	  if ( strcmp(sData,"-spall") == 0 ) {
		  numData = 1;
		  if (OPS_GetDoubleInput(&numData, &xn_cr) != 0) {
		    opserr << "WARNING invalid input for -spall option, want: -spall xn_cr \n";
		    return 0;
		  }
		  if ( xn_cr <= 1.0 ) {
			  opserr << "WARNING xn_cr should be greater than 1.0 \n";
			  return 0;
		  }
      if ( dData[4] == 0.0 ) {
        opserr << "WARNING If modeling spalling, rn_post cannot be zero \n";
        return 0;
      }

	  } else {
		  opserr << "WARNING unknown option " << sData << "\n";
	  }
  }


  theMaterial = new changManderConcrete01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], xn_cr);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type changManderConcrete01\n";
    return 0;
  }

  return theMaterial;
}




changManderConcrete01::changManderConcrete01(int tag, double i1, double i2, double i3, double i4, double i5, double i6, double i7, double i8, double i9, double i10)
:UniaxialMaterial(tag,MAT_TAG_changManderConcrete01),
 Fc_n(i1), ec_n(i2), Ec(i3), r_n_pre(i4), r_n_post(i5),
 Fc_p(i6), ec_p(i7), r_p(i8), x_p_cr(i9), x_n_cr(i10)
{
	double ycr, zcr, n;

	n = fabs(Ec * ec_p / Fc_p);
	tsaiEquation(x_p_cr, r_p, n, ycr, zcr);
	ecrk_p = x_p_cr*ec_p - (Fc_p * ycr)/(Ec * zcr);

	if ( x_n_cr >= 1.0 ) {
		modelSpalling = true;
		n = fabs(Ec * ec_n / Fc_n);
		tsaiEquation(x_n_cr, r_n_post, n, ycr, zcr);
		espall = x_n_cr*ec_n - (Fc_n * ycr)/(Ec * zcr);
	} else {
		modelSpalling = false;
		espall = 0.0;
	}

	this->revertToStart();
}

changManderConcrete01::changManderConcrete01()
:UniaxialMaterial(0,MAT_TAG_changManderConcrete01),
Fc_n(0.0), ec_n(0.0), Ec(0.0), r_n_pre(0.0), r_n_post(0.0),
Fc_p(0.0), ec_p(0.0), r_p(0.0), x_p_cr(0.0), x_n_cr(0.0)
{
	ecrk_p = 0.0;
	modelSpalling = false;
	espall = 0.0;

	this->revertToStart();
}

changManderConcrete01::~changManderConcrete01()
{
  // does nothing
}

int
changManderConcrete01::setTrialStrain(double strain, double strainRate)
{
	// Return state variables to last commited values
	backToCommitStateVar();

	// Define trial strain and strain increment
	double strain_incr;
	Tstrain = strain; // For use with other elements
	strain_incr = Tstrain - Cstrain; // Works for both

	if(modelSpalling == true && isSpalled == true){
		Trule = 5;
	} else if(isCracked == false) { // Pre-Cracking
		switch ( Crule ) {
		case 0:
			// It should only arrive here before loading has started
			if( fabs(Tstrain) < SMALL_STRAIN ) {
				Trule = 0;
			} else if ( modelSpalling == true && Tstrain <= espall ){
				isSpalled = true;
				Trule = 5;
			} else if( Tstrain <= 0.0 ) {
				// Compressive Strain
				Trule = 1;
			} else if ( Tstrain <= ecrk_p + eo ) {
				// Tensile Strain less than cracking
				Trule = 2;
			} else {
				// Tensile strain greater than cracking
				isCracked = true;
				Trule = 6;
			}
			break;

		case 1:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive loading
				if ( modelSpalling == true && Tstrain <= espall ) {
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Unloading
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope(er1,fr1,scratch);
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;

				// Start of shifting the tensile branch
				double eo_old = eo;

				// if xu_p < xu_n force a new reversal from tensile branch and reset eo
				if ( fabs((er2-eo)/ec_p) < fabs(er1/ec_n) ) {
					er2 = er1 * ec_p / ec_n;
					eo = 0.0;
					eo_old = 0.0;
					// compute reversal for new er2
					double fr2, scratch;
					positiveEnvelope(er2,fr2,scratch);
					Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
					Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
					delta_f_p = 0.15*fr2;
					delta_e_p = 0.22*(er2-eo);
					epl_p = er2 - fr2 / Esec_p;
					fnew_p = fr2 - delta_f_p;
					if ( fnew_p == 0 ) {
						isCracked = true;
					}
					if ( fabs(er2-epl_p) <= SMALL_NUMBER ) {
						Enew_p = Ec;
					} else {
						Enew_p = fnew_p / (er2 - epl_p );
					}
					//trialConcData<<" Enew_p is calculated here  Enew_p: "<<Enew_p<<" fnew_p: "<<fnew_p<<" er2: "<<er2<<" epl_p: "<<epl_p<<endl;
					ere_p = er2 + delta_e_p;
				}

				// compute new eo
				double fr2;
				positiveEnvelope(er2,fr2,scratch);
				double Deo = 2*fr2/(Esec_p+Epl_p);
				eo = epl_n - (er2-eo-Deo);

				// shift all necessary strains (some may not be necessary)
				er2 += (eo-eo_old);
				er4 += (eo-eo_old);
				er6 += (eo-eo_old);
				er9 += (eo-eo_old);
				ea += (eo-eo_old);
				e82target += (eo-eo_old);
				ere_p += (eo-eo_old);
				ere_str_p += (eo-eo_old);
				epl_p += (eo-eo_old);
				// End of shifting the tensile branch


				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 2:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				if ( Tstrain <= ecrk_p + eo ) {
					// Before cracking
					Trule = 2;
				} else {
					// After Cracking
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Unloading
				er2 = Cstrain;
				double fr2, scratch;
				positiveEnvelope(er2,fr2,scratch);
				Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
				Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
				delta_f_p = 0.15*fr2;
				delta_e_p = 0.22*(er2-eo);
				epl_p = er2 - fr2 / Esec_p;
				fnew_p = fr2 - delta_f_p;
				if ( fabs(er2-epl_p) <= SMALL_NUMBER ) {
					Enew_p = Ec;
				} else {
					Enew_p = fnew_p / (er2 - epl_p );
				}
				//trialConcData<<" Enew_p is calculated here  Enew_p: "<<Enew_p<<" fnew_p: "<<fnew_p<<" er2: "<<er2<<" epl_p: "<<epl_p<<endl;
				ere_p = er2 + delta_e_p;
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					f4target = Cstress;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 3:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					// target remains as it is
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from Rule 3
				er3 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				double fr3 = Cstress;
				fnew_str_n = fr1 - delta_f_n * ( er1 - er3 )/( er1 - epl_n );
				Enew_str_n = ( fnew_str_n - fr3 ) / ( er1 - er3 );
				ere_str_n = er1 + delta_e_n * ( er1 - er3 ) / ( er1 - epl_n );

				if ( Tstrain >= er1 ) {
					Trule = 16;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 16 after partial unloading
				} else if ( Tstrain >= ere_str_n )  {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}

			break;

		case 4:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive unloading
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					// target remains as it is
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal
				er4 = Cstrain;
				double fr2, scratch;
				positiveEnvelope( er2, fr2, scratch );
				double fr4 = Tstress;
				fnew_str_p = fr2 - delta_f_p * ( er2 - er4 )/( er2 - epl_p );
				Enew_str_p = ( fnew_str_p - fr4 ) / ( er2 - er4 );
				ere_str_p = er2 + delta_e_p * ( er2 - er4 ) / ( er2 - epl_p );

				if ( Tstrain <= er2 ) {
					Trule = 17;
					f81target = fnew_str_p;
					E81target = Enew_str_p;
					e82target = ere_str_p; // Because arriving on Rule 17 after partial unloading
				} else if ( Tstrain <= ere_str_p ) {
					Trule = 8;
					f81target = fnew_str_p;
					E81target = Enew_str_p;
					e82target = ere_str_p; // Because arriving on Rule 8 after partial unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 7:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive
				if ( Tstrain >= e72target ) {
					Trule = 7;
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from Second Branch
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 8:
			if ( strain_incr >= SMALL_STRAIN ) {
				// Continued tensile
				if ( Tstrain <= e82target ) {
					Trule = 8;
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from Second Branch of Rule 8
				er2 = Cstrain;
				double fr2, scratch;
				positiveEnvelope(er2,fr2,scratch);
				Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
				Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
				delta_f_p = 0.15*fr2;
				delta_e_p = 0.22*er2;
				epl_p = er2 - fr2 / Esec_p;
				fnew_p = fr2 - delta_f_p;
				if ( fabs(er2-epl_p) <= SMALL_NUMBER ) {
					Enew_p = Ec;
				} else {
					Enew_p = fnew_p / (er2 - epl_p );
				}
				//trialConcData<<" Enew_p is calculated here  Enew_p: "<<Enew_p<<" fnew_p: "<<fnew_p<<" er2: "<<er2<<" epl_p: "<<epl_p<<endl;
				ere_p = er2 + delta_e_p;
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					f4target = Cstress;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain <= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 9:
			if ( strain_incr >= -SMALL_STRAIN ) {
				// Continued tensile loading
				if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from Rule 9
				er9 = Cstrain;
				eb = er1 - ((er9 - epl_n)/(er2 - epl_n))*(er1-epl_p);

				if ( Tstrain >= eb ) {
					Trule = 11;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 10:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive
				if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal
				er10 = Cstrain;
				ea = epl_n + ((er1 - er10)/(er1 - epl_p))*(er2-epl_n);

				if ( Tstrain <= ea ) {
					Trule = 12;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain < ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 11:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive
				if ( Tstrain >= eb ) {
					Trule = 11;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal
				if ( Tstrain <= er9 ) {
					Trule = 11;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain < ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 12:
			if ( strain_incr >= -SMALL_STRAIN ) {
				// Continued tensile loading
				if ( Tstrain < ea ) {
					Trule = 12;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from Rule 12
				if ( Tstrain >= er10 ) {
					Trule = 12;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 16:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er1 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain > e72target ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from First Branch of Rule 7, which is renamed Rule 16
				if ( Tstrain <= er3 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain <= epl_n ) {
					Trule = 3;
					// target remains as it is
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					isCracked = true;
					Trule = 6;
				}
			}
			break;

		case 17:
			if ( strain_incr >= SMALL_STRAIN ) {
				// Continued tensile
				if ( Tstrain <= er2 ) {
					Trule = 17;
					// Targets stay as is
				} else if ( Tstrain <= e82target ) {
					Trule = 8;
					f81target = fnew_str_p;
					E81target = Enew_str_p;
					e82target = ere_str_p; // Because arriving on Rule 8 after partial unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else {
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from First Branch
				if ( Tstrain >= er4 ) {
					Trule = 17;
					// Targets stay as is
				} else if ( Tstrain >= epl_p ) {
					Trule = 4;
					// target remains as it is
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		default:
			// It should not arrive here
			break;

	} // End of switch for Crule (for pre-cracking region)

	} else { // Post-Cracking
		//trialConcData<<"I got to the post cracking area"<<endl;

		switch ( Crule ) {
		case 1:
			//trialConcData<<"I got to the rule 1 area"<<endl;
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive loading
				if ( modelSpalling == true && Tstrain <= espall ) {
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else { // Unloading
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope(er1,fr1,scratch);
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				//trialConcData<<"I got to the unloading from rule 1 area"<<endl;
				//trialConcData<<" er1: "<<er1<<" fr1: "<<fr1<<" Esec_n: "<<Esec_n<<" Epl_n: "<<Epl_n<<" delta_f_n: "<<delta_f_n<<endl;
				//trialConcData<<" delta_e_n: "<<delta_e_n<<" epl_n: "<<epl_n<<" fnew_n: "<<fnew_n<<" Enew_n: "<<Enew_n<<" ere_n: "<<ere_n<<endl;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else {
					Trule = 6;
				}
			}
			break;

		case 3:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued unloading
				if ( Tstrain <= epl_n ) {
					Trule = 3;
				} else {
					Trule = 6;
				}
			} else { // Reloading (reversal)
				er3 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				double fr3 = Cstress;
				fnew_str_n = fr1 - delta_f_n * ( er1 - er3 )/( er1 - epl_n );
				Enew_str_n = ( fnew_str_n - fr3 ) / ( er1 - er3 );
				ere_str_n = er1 + delta_e_n * ( er1 - er3 ) / ( er1 - epl_n );

				if ( Tstrain >= er1 ) {
					Trule = 16;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else if ( Tstrain >= ere_str_n ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 6:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				Trule = 6;
			} else { // reversal from 6
				er6 = Cstrain;
				if ( Tstrain >= er1 ) {
					Trule = 13;
				} else if (Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 7:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > e72target ) {
					Trule = 7;
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from Second Branch
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else {
					Trule = 6;
				}
			}
			break;

		case 13:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain >= er1 ) {
					Trule = 13;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
					//trialConcData << "I just got to Rule 7 from Rule 13, f71target" << f71target << "E71target" << E71target<< "e72target" << e72target << endl;
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				er13 = Cstrain;
				eb = er13 - Cstress/Esec_n;
				if ( Tstrain < eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			}
			break;

		case 14:
			if ( strain_incr >= -SMALL_STRAIN ) {
				if ( Tstrain < eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			} else { // Reversal from Rule 14
				er14 = Cstrain;
				if ( Tstrain >= er13 ) {
					Trule = 15;
				} else if ( Tstrain > er1 ) {
					Trule = 13;
				} else if ( Tstrain > ere_n ) {
					Trule = 7;
					// will be the second branch
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			}
			break;

		case 15:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er13) {
					Trule = 15;
				} else if ( Tstrain > er1 ) {
					Trule = 13;
				} else if ( Tstrain > ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				if ( Tstrain <= er14 ) {
					Trule = 15; //stay in Rule 15
				} else if ( Tstrain <= eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			}
			break;

		case 16:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er1 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain > e72target ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else if ( modelSpalling == true && Tstrain <= espall ){
					isSpalled = true;
					Trule = 5;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from First Branch of Rule 7, which is renamed Rule 16
				if ( Tstrain <= er3 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain <= epl_n ) {
					Trule = 3;
				} else {
					Trule = 6;
				}
			}
			break;

		default:
			// It should not arrive here
			break;

	} // End of switch for Crule
	} // End of if statment for isCracked

	// Now that we have found the new rule and updated any state variables, we find the new stress and tangent
	switch ( Trule ) {
		case 0:
			// Loading hasn't started yet
			Tstress = 0.0;
			Ttangent = Ec;
			break;
		case 1:
			negativeEnvelope(Tstrain,Tstress,Ttangent);
			break;
		case 2:
			positiveEnvelope(Tstrain,Tstress,Ttangent);
			break;
		case 3:
			transitionCurve(Tstrain, Tstress, Ttangent, er1, f3target, Ec, epl_n, 0.0, Epl_n, Trule);
			break;
		case 4:
			transitionCurve(Tstrain, Tstress, Ttangent, er2, f4target, Ec, epl_p, 0.0, Epl_p, Trule);
			break;
		case 5:
			Tstress = SMALL_STRESS;
			Ttangent = SMALL_TANGENT;
			break;
		case 6:
			Tstress = SMALL_STRESS;
			Ttangent = SMALL_TANGENT;
			break;
		case 7:
			{
				// Second Branch of Rule 7 remains Rule 7
				// First Branch has been renamed Rule 16
				double ff, Ef;
				negativeEnvelope(e72target,ff,Ef); // Calculate target stress and tangent
				transitionCurve(Tstrain, Tstress, Ttangent, er1, f71target, E71target, e72target, ff, Ef, Trule);
			}
			break;
		case 8:
			{
				// Second Branch of RUle 8 remains Rule 8
				// First Branch has been renamed Rule 17
				double ff, Ef;
				positiveEnvelope(e82target,ff,Ef); // Calculate target stress and tangent
				transitionCurve(Tstrain, Tstress, Ttangent, er2, f81target, E81target, e82target, ff, Ef, Trule);
			}
			break;
		case 9:
			transitionCurve(Tstrain, Tstress, Ttangent, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p, Trule);
			break;
		case 10:
			transitionCurve(Tstrain, Tstress, Ttangent, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n, Trule);
			break;
		case 11:
			{
				double fr9, fb, Eb, scratch;
				//Find the inital stress (using Rule 9)
				transitionCurve(er9, fr9, scratch, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p, Trule);
				//Find the final stress and tangent (using Rule 10)
				transitionCurve(eb, fb, Eb, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n, Trule);
				//Set Tstress and Ttangent (Using Rule 11)
				transitionCurve(Tstrain, Tstress, Ttangent, er9, fr9, Ec, eb, fb, Eb, Trule);
			}
			break;
		case 12:
			{
				double fr10, fa, Ea, scratch;
				//Find the inital stress (using Rule 10)
				transitionCurve(er10, fr10, scratch, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n, Trule);
				//Find the final stress and tangent (using Rule 9)
				transitionCurve(ea, fa, Ea, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p, Trule);
				//Set Tstress and Ttangent (Using Rule 12)
				transitionCurve(Tstrain, Tstress, Ttangent, er10, fr10, Ec, ea, fa, Ea, Trule);
			}
			break;
		case 13:
			transitionCurve(Tstrain, Tstress, Ttangent, er6, 0.0, 0.0, er1, fnew_n, Enew_n, Trule);
			break;
		case 14:
			{
				double fr13, scratch;
				//Find the inital stress (Using Rule 13)
				transitionCurve(er13, fr13, scratch, er6, 0.0, 0.0, er1, fnew_n, Enew_n, Trule);
				//Set Tstress and Ttangent (Using Rule 14)
				transitionCurve(Tstrain, Tstress, Ttangent, er13, fr13, Ec, eb, 0.0, 0.0, Trule);
			}
			break;
		case 15:
			{
				double fr13, fr14, Er13, scratch;
				//Find fr13 and Er13 (Using Rule 13)
				transitionCurve(er13, fr13, Er13, er6, 0.0, 0.0, er1, fnew_n, Enew_n, Trule);
				//Find fr14 (Using Rule 14)
				transitionCurve(er14, fr14, scratch, er13, fr13, Ec, eb, 0.0, 0.0, Trule);
				//Set Tstress and Ttangent (Using Rule 15)
				transitionCurve(Tstrain, Tstress, Ttangent, er14, fr14, Ec, er13, fr13, Er13, Trule);
			}
			break;
		case 16:
			{
				// First Branch of Rule 7 (per Chang and Mander) is renamed Rule 16
				double fr3, scratch;
				transitionCurve(er3, fr3, scratch, er1, f3target, Ec, epl_n, 0.0, Epl_n, Trule);
					// starting stress for Rule 16 (first branch) comes from Rule 3
				transitionCurve(Tstrain, Tstress, Ttangent, er3, fr3, Ec, er1, f71target, E71target, Trule);
			}
			break;
		case 17:
			{
				// First Branch of Rule 8 (per Chang and Mander) is renamed Rule 17
				double fr4, scratch;
				transitionCurve(er4, fr4, scratch, er2, f4target, Ec, epl_p, 0.0, Epl_p, Trule);
					// starting stress for Rule 8 (first branch) comes from Rule 4
				transitionCurve(Tstrain, Tstress, Ttangent, er4, fr4, Ec, er2, f81target, E81target, Trule);
			}
			break;
		default:
			// It should not arrive here
			break;
	}

	return 0;
}

double
changManderConcrete01::getStrain(void)
{
  return Tstrain;
}

double
changManderConcrete01::getStress(void)
{
  return Tstress;
}


double
changManderConcrete01::getTangent(void)
{
  return Ttangent;
}

double
changManderConcrete01::getInitialTangent(void)
{
  return Ec;
}

int
changManderConcrete01::commitState(void)
{
	this->commitStateVar();
	return 0;
}


int
changManderConcrete01::revertToLastCommit(void)
{
	this->backToCommitStateVar();
	return 0;
}


int
changManderConcrete01::revertToStart(void)
{
	isCracked = false;	CisCracked = false;
	isSpalled = false;	CisSpalled = false;
	Trule = 0; Crule = 0;
	//isCracked = true;	CisCracked = true; // For Testing in Post-Cracking
	//Trule = 0; Crule = 1; // For Testing in Post-Cracking

	Tstrain = 0.0; Cstrain = 0.0;
	Tstress = 0.0; Cstress = 0.0;
	Ttangent = 0.0; Ctangent = 0.0;
	eo = 0.0; Ceo = 0.0;
	er1 = 0.0; er3 = 0.0; er6 = 0.0; er13 = 0.0; er14 = 0.0;
	Cer1 = 0.0; Cer3 = 0.0; Cer6 = 0.0; Cer13 = 0.0; Cer14 = 0.0;
	eb = 0.0; Ceb = 0.0;
	Esec_n = 0.0; CEsec_n = 0.0;
	Epl_n = 0.0; CEpl_n = 0.0;
	delta_f_n = 0.0; Cdelta_f_n = 0.0;
	delta_e_n = 0.0; Cdelta_e_n = 0.0;
	epl_n = 0.0; Cepl_n = 0.0;
	fnew_n = 0.0; Cfnew_n = 0.0;
	Enew_n = 0.0; CEnew_n = 0.0;
	ere_n = 0.0; Cere_n = 0.0;
	fnew_str_n = 0.0; Cfnew_str_n = 0.0;
	Enew_str_n = 0.0; CEnew_str_n = 0.0;
	ere_str_n = 0.0; Cere_str_n = 0.0;
	f71target = 0.0; Cf71target = 0.0;
	E71target = 0.0; CE71target = 0.0;
	e72target = 0.0; Ce72target = 0.0;
	f3target = 0.0; Cf3target = 0.0;
	f81target = 0.0; Cf81target = 0.0;
	E81target = 0.0; CE81target = 0.0;
	e82target = 0.0; Ce82target = 0.0;
	f4target = 0.0; Cf4target = 0.0;
	delta_e_p = 0.0; Cdelta_e_p = 0.0;
	delta_f_p = 0.0; Cdelta_f_p = 0.0;
	fnew_p = 0.0; Cfnew_p = 0.0;
	Enew_p = 0.0; CEnew_p = 0.0;
	ere_p = 0.0; Cere_p = 0.0;
	fnew_str_p = 0.0; Cfnew_str_p = 0.0;
	Enew_str_p = 0.0; CEnew_str_p = 0.0;
	ere_str_p = 0.0; Cere_str_p = 0.0;
	epl_p = 0.0; Cepl_p = 0.0;
	Epl_p = 0.0; CEpl_p = 0.0;
	er2 = 0.0; er4 = 0.0; er9 = 0.0; er10 = 0.0;
	Cer2 = 0.0; Cer4 = 0.0; Cer9 = 0.0; Cer10 = 0.0;
	Esec_p = 0.0; CEsec_p = 0.0;
	ea = 0.0; Cea = 0.0;

	return 0;
}


UniaxialMaterial *
changManderConcrete01::getCopy(void)
{
  changManderConcrete01 *theCopy = new changManderConcrete01(this->getTag(), Fc_n, ec_n, Ec, r_n_pre, r_n_post, Fc_p, ec_p, r_p, x_p_cr, x_n_cr);

  theCopy->isCracked = this->isCracked;
  theCopy->CisCracked = this->CisCracked;
  theCopy->isSpalled = this->isSpalled;
  theCopy->CisSpalled = this->CisSpalled;
  theCopy->Trule = this->Trule;
  theCopy->Crule = this->Crule;
  theCopy->Tstrain = this->Tstrain;
  theCopy->Cstrain = this->Cstrain;
  theCopy->Tstress = this->Tstress;
  theCopy->Cstress = this->Cstress;
  theCopy->Ttangent = this->Ttangent;
  theCopy->Ctangent = this->Ctangent;
  theCopy->eo = this->eo;
  theCopy->Ceo = this->Ceo;
  theCopy->er1 = this->er1;
  theCopy->er3 = this->er3;
  theCopy->er6 = this->er6;
  theCopy->er13 = this->er13;
  theCopy->er14 = this->er14;
  theCopy->Cer1 = this->Cer1;
  theCopy->Cer3 = this->Cer3;
  theCopy->Cer6 = this->Cer6;
  theCopy->Cer13 = this->Cer13;
  theCopy->Cer14 = this->Cer14;
  theCopy->eb = this->eb;
  theCopy->Ceb = this->Ceb;
  theCopy->Esec_n = this->Esec_n;
  theCopy->CEsec_n = this->CEsec_n;
  theCopy->Epl_n = this->Epl_n;
  theCopy->CEpl_n = this->CEpl_n;
  theCopy->delta_f_n = this->delta_f_n;
  theCopy->Cdelta_f_n = this->delta_f_n;
  theCopy->delta_e_n = this->delta_f_n;
  theCopy->Cdelta_e_n = this->Cdelta_e_n;
  theCopy->epl_n = this->epl_n;
  theCopy->Cepl_n = this->Cepl_n;
  theCopy->fnew_n = this->fnew_n;
  theCopy->Cfnew_n = this->Cfnew_n;
  theCopy->Enew_n = this->Enew_n;
  theCopy->CEnew_n = this->CEnew_n;
  theCopy->ere_n = this->ere_n;
  theCopy->Cere_n = this->Cere_n;
  theCopy->fnew_str_n = this->fnew_str_n;
  theCopy->Cfnew_str_n = this->Cfnew_str_n;
  theCopy->Enew_str_n = this->Enew_str_n;
  theCopy->CEnew_str_n = this->CEnew_str_n;
  theCopy->ere_str_n = this->ere_str_n;
  theCopy->Cere_str_n = this->Cere_str_n;
  theCopy->f71target = this->f71target;
  theCopy->Cf71target = this->Cf71target;
  theCopy->E71target = this->E71target;
  theCopy->CE71target = this->CE71target;
  theCopy->e72target = this->e72target;
  theCopy->Ce72target = this->Ce72target;
  theCopy->f3target = this->f3target;
  theCopy->Cf3target = this->Cf3target;
  theCopy->f81target = this->f81target;
  theCopy->Cf81target = this->Cf81target;
  theCopy->E81target = this->E81target;
  theCopy->CE81target = this->CE81target;
  theCopy->e82target = this->e82target;
  theCopy->Ce82target = this->Ce82target;
  theCopy->f4target = this->f4target;
  theCopy->Cf4target = this->Cf4target;
  theCopy->delta_e_p = this->delta_e_p;
  theCopy->Cdelta_e_p = this->Cdelta_e_p;
  theCopy->delta_f_p = this->delta_f_p;
  theCopy->Cdelta_f_p = this->Cdelta_f_p;
  theCopy->fnew_p = this->fnew_p;
  theCopy->Cfnew_p = this->Cfnew_p;
  theCopy->Enew_p = this->Enew_p;
  theCopy->CEnew_p = this->CEnew_p;
  theCopy->ere_p = this->ere_p;
  theCopy->Cere_p = this->Cere_p;
  theCopy->fnew_str_p = this->fnew_str_p;
  theCopy->Cfnew_str_p = this->Cfnew_str_p;
  theCopy->Enew_str_p = this->Enew_str_p;
  theCopy->CEnew_str_p = this->CEnew_str_p;
  theCopy->ere_str_p = this->ere_str_p;
  theCopy->Cere_str_p = this->Cere_str_p;
  theCopy->epl_p = this->epl_p;
  theCopy->Cepl_p = this->Cepl_p;
  theCopy->Epl_p = this->Epl_p;
  theCopy->CEpl_p = this->CEpl_p;
  theCopy->er2 = this->er2;
  theCopy->er4 = this->er4;
  theCopy->er9 = this->er9;
  theCopy->er10 = this->er10;
  theCopy->Cer2 = this->Cer2;
  theCopy->Cer4 = this->Cer4;
  theCopy->Cer9 = this->Cer9;
  theCopy->Cer10 = this->Cer10;
  theCopy->Esec_p = this->Esec_p;
  theCopy->CEsec_p = this->CEsec_p;
  theCopy->ea = this->ea;
  theCopy->Cea = this->Cea;

  return theCopy;
}


int
changManderConcrete01::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(106);
  data(0) = this->getTag();
  data(1) = Fc_n;
  data(2) = ec_n;
  data(3) = Ec;
  data(4) = r_n_pre;
  data(5) = r_n_post;
  data(6) = Fc_p;
  data(7) = ec_p;
  data(8) = r_p;
  data(9) = x_p_cr;
  data(10) = 0.0;
  data(11) = ecrk_p;

  if (isCracked)
    data(12) = 1;
  else
	data(12) = 0;

  if (CisCracked)
    data(13) = 1;
  else
	data(13) = 0;

  data(14) = Trule;
  data(15) = Crule;
  data(16) = Tstrain;
  data(17) = Cstrain;
  data(18) = Tstress;
  data(19) = Cstress;
  data(20) = Ttangent;
  data(21) = Ctangent;
  data(22) = eo;
  data(23) = Ceo;
  data(24) = er1;
  data(25) = er3;
  data(26) = er6;
  data(27) = er13;
  data(28) = er14;
  data(29) = Cer1;
  data(30) = Cer3;
  data(31) = Cer6;
  data(32) = Cer13;
  data(33) = Cer14;
  data(34) = eb;
  data(35) = Ceb;
  data(36) = Esec_n;
  data(37) = CEsec_n;
  data(38) = Epl_n;
  data(39) = CEpl_n;
  data(40) = delta_f_n;
  data(41) = Cdelta_f_n;
  data(42) = delta_e_n;
  data(43) = Cdelta_e_n;
  data(44) = epl_n;
  data(45) = Cepl_n;
  data(46) = fnew_n;
  data(47) = Cfnew_n;
  data(48) = Enew_n;
  data(49) = CEnew_n;
  data(50) = ere_n;
  data(51) = Cere_n;
  data(52) = fnew_str_n;
  data(53) = Cfnew_str_n;
  data(54) = Enew_str_n;
  data(55) = CEnew_str_n;
  data(56) = ere_str_n;
  data(57) = Cere_str_n;
  data(58) = f71target;
  data(59) = Cf71target;
  data(60) = E71target;
  data(61) = CE71target;
  data(62) = e72target;
  data(63) = Ce72target;
  data(64) = f3target;
  data(65) = Cf3target;
  data(66) = f81target;
  data(67) = Cf81target;
  data(68) = E81target;
  data(69) = CE81target;
  data(70) = e82target;
  data(71) = Ce82target;
  data(72) = f4target;
  data(73) = Cf4target;
  data(74) = delta_e_p;
  data(75) = Cdelta_e_p;
  data(76) = delta_f_p;
  data(77) = Cdelta_f_p;
  data(78) = fnew_p;
  data(79) = Cfnew_p;
  data(80) = Enew_p;
  data(81) = CEnew_p;
  data(82) = ere_p;
  data(83) = Cere_p;
  data(84) = fnew_str_p;
  data(85) = Cfnew_str_p;
  data(86) = Enew_str_p;
  data(87) = CEnew_str_p;
  data(88) = ere_str_p;
  data(89) = Cere_str_p;
  data(90) = epl_p;
  data(91) = Cepl_p;
  data(92) = Epl_p;
  data(93) = CEpl_p;
  data(94) = er2;
  data(95) = er4;
  data(96) = er9;
  data(97) = er10;
  data(98) = Cer2;
  data(99) = Cer4;
  data(100) = Cer9;
  data(101) = Cer10;
  data(102) = Esec_p;
  data(103) = CEsec_p;
  data(104) = ea;
  data(105) = Cea;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "changManderConcrete01::sendSelf() - failed to send data\n";

  return res;
}

int
changManderConcrete01::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "changManderConcrete01::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));

    Fc_n = data(1);
    ec_n = data(2);
    Ec = data(3);
    r_n_pre = data(4);
    r_n_post = data(5);
    Fc_p = data(6);
    ec_p = data(7);
    r_p = data(8);
    x_p_cr = data(9);

    ecrk_p = data(11);

	if (data(12) == 1)
		isCracked = true;
	else
		isCracked = false;

	if (data(13) == 1)
		CisCracked = true;
	else
		CisCracked = false;

	Trule = data(14);
	Crule = data(15);
	Tstrain = data(16);
	Cstrain = data(17);
	Tstress = data(18);
	Cstress = data(19);
	Ttangent = data(20);
	Ctangent = data(21);
	eo = data(22);
	Ceo = data(23);
	er1 = data(24);
	er3 = data(25);
	er6 = data(26);
	er13 = data(27);
	er14 = data(28);
	Cer1 = data(29);
	Cer3 = data(30);
	Cer6 = data(31);
	Cer13 = data(32);
	Cer14 = data(33);
	eb = data(34);
	Ceb = data(35);
	Esec_n = data(36);
	CEsec_n = data(37);
	Epl_n = data(38);
	CEpl_n = data(39);
	delta_f_n = data(40);
	Cdelta_f_n = data(41);
	delta_e_n = data(42);
	Cdelta_e_n = data(43);
	epl_n = data(44);
	Cepl_n = data(45);
	fnew_n = data(46);
	Cfnew_n = data(47);
	Enew_n = data(48);
	CEnew_n = data(49);
	ere_n = data(50);
	Cere_n = data(51);
	fnew_str_n = data(52);
	Cfnew_str_n = data(53);
	Enew_str_n = data(54);
	CEnew_str_n = data(55);
	ere_str_n = data(56);
	Cere_str_n = data(57);
	f71target = data(58);
	Cf71target = data(59);
	E71target = data(60);
	CE71target = data(61);
	e72target = data(62);
	Ce72target = data(63);
	f3target = data(64);
	Cf3target = data(65);
	f81target = data(66);
	Cf81target = data(67);
	E81target = data(68);
	CE81target = data(69);
	e82target = data(70);
	Ce82target = data(71);
	f4target = data(72);
	Cf4target = data(73);
	delta_e_p = data(74);
	Cdelta_e_p = data(75);
	delta_f_p = data(76);
	Cdelta_f_p = data(77);
	fnew_p = data(78);
	Cfnew_p = data(79);
	Enew_p = data(80);
	CEnew_p = data(81);
	ere_p = data(82);
	Cere_p = data(83);
	fnew_str_p = data(84);
	Cfnew_str_p = data(85);
	Enew_str_p = data(86);
	CEnew_str_p = data(87);
	ere_str_p = data(88);
	Cere_str_p = data(89);
	epl_p = data(90);
	Cepl_p = data(91);
	Epl_p = data(92);
	CEpl_p = data(93);
	er2 = data(94);
	er4 = data(95);
	er9 = data(96);
	er10 = data(97);
	Cer2 = data(98);
	Cer4 = data(99);
	Cer9 = data(100);
	Cer10 = data(101);
	Esec_p = data(102);
	CEsec_p = data(103);
	ea = data(104);
	Cea = data(105);
  }

  return res;
}

void
changManderConcrete01::Print(OPS_Stream &s, int flag)
{
	s<<"changManderConcrete01, tag: "<<this->getTag()<<endln;
	s<<" Ec:      "<<Ec<<endln;
	s<<" Fcn:     "<<Fc_n<<endln;
	s<<" ecn:     "<<ec_n<<endln;
	s<<" rn_pre:  "<<r_n_pre<<endln;
	s<<" rn_post: "<<r_n_post<<endln;
	s<<" xcrn:    "<<x_n_cr<<endln;
	s<<" Fcp:     "<<Fc_p<<endln;
	s<<" ecp:     "<<ec_p<<endln;
	s<<" rp:      "<<r_p<<endln;
	s<<" xcrp:    "<<x_p_cr<<endln;
	return;
}


void
changManderConcrete01::negativeEnvelope(double strain, double &stress, double &tangent)
{
	double x_n, n_n, y, z;
	x_n = strain/ec_n;
	n_n = fabs(Ec * ec_n / Fc_n);

	if ( modelSpalling ) {
	  if ( strain > 0.0 ) {
	    stress = 0.0;
	    tangent = 0.0;
	  } else if ( strain > ec_n ) {
	    tsaiEquation(x_n, r_n_pre, n_n, y, z);
	    stress = Fc_n * y;
	    tangent = Ec * z;
	  } else if ( strain > ec_n*x_n_cr ) {
	    tsaiEquation(x_n, r_n_post, n_n, y, z);
	    stress = Fc_n * y;
	    tangent = Ec * z;
	  } else if ( strain > espall ) {
	    // linear softening branch between critical and spalled
	    double y_cr, z_cr;
	    tsaiEquation(x_n_cr, r_n_post, n_n, y_cr, z_cr);
	    stress = Fc_n * (y_cr + n_n * z_cr * (x_n - x_n_cr));
	    tangent = Ec * z_cr;
	  } else {
      stress = 0.0;
      tangent = 0.0;
	  }
	} else {
    if ( strain > 0.0 ) {
      stress = 0.0;
      tangent = 0.0;
    } else if ( strain > ec_n ) {
      tsaiEquation(x_n, r_n_pre, n_n, y, z);
      stress = Fc_n * y;
      tangent = Ec * z;
    } else {
      if ( r_n_post == 0 ) {
        stress = Fc_n;
        tangent = 0.0;
      } else {
        tsaiEquation(x_n, r_n_post, n_n, y, z);
        stress = Fc_n * y;
        tangent = Ec * z;
      }
    }
	}
	return;
}

void
changManderConcrete01::positiveEnvelope(double strain, double &stress, double &tangent)
{
	double n_p, x_p;
	n_p = fabs(Ec * ec_p / Fc_p);
	x_p = (strain - eo)/ec_p;

	if( strain - eo <= 0 ) {
		// it should not get here but if it does return zero
		stress = 0.0;
		tangent = 0.0;
	} else if ( strain - eo <= x_p_cr*ec_p ) {
		// nonlinear pre-critical branch
		double y, z;
		tsaiEquation(x_p, r_p, n_p, y, z);
		stress = Fc_p * y;
		tangent = Ec * z;
	} else if ( strain - eo <= ecrk_p ) {
		// linear softening branch between critical and cracked
		double y_cr, z_cr;
		tsaiEquation(x_p_cr, r_p, n_p, y_cr, z_cr);
		stress = Fc_p * (y_cr + n_p * z_cr * (x_p - x_p_cr));
		tangent = Ec * z_cr;
	} else {
		// cracked, return positive small stress
		stress = SMALL_STRESS;
		tangent = SMALL_TANGENT;
	}
	return;
}


void
changManderConcrete01::transitionCurve(double Tstrain, double &Tstress, double &Ttangent,
		double ei, double fi, double Ei, double ef, double ff, double Ef, int rule)
{
	double R, A, Esec;

	if ( fabs (ei/ef - 1) < 1.0e-4) {
		// Inital and final strain are too close together, linear function
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		return;
	}

	Esec = ( ff - fi ) / ( ef - ei );

	if( fabs(Esec-Ei) <= SMALL_NUMBER ) {
		R = 0.0;
	} else {
		R = ( Ef - Esec ) / ( Esec - Ei );
	}

	if( R <= 0.0 ) {
		// Smooth transition curve without change of curvature not possible, linear function
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		return;
	}

	A = (Esec-Ei) * pow( fabs((Tstrain-ei)/(ef - ei)) , R );
	Tstress = fi + (Tstrain-ei) * ( Ei + A );
	Ttangent = Ei + A*(R+1);
	return;
}

void
changManderConcrete01::tsaiEquation(double x, double r, double n, double &y, double &z)
{
	double D;
	if( r == 1 ) {
	    D = 1 + ( n - 1 + log( x ) ) * x;
	} else {
	    D = 1 + ( n - r / ( r - 1 ) ) * x + pow ( x , r ) / ( r - 1 );
	}
	y = ( n * x / D );
	z = ( 1 - pow( x, r ) ) / ( D * D );
	return;
}

void
changManderConcrete01::backToCommitStateVar(void){
	isCracked = CisCracked;
	isSpalled = CisSpalled;
	Trule = Crule;
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
	eo = Ceo;
	er1 = Cer1;
	er3 = Cer3;
	er6 = Cer6;
	er13 = Cer13;
	er14 = Cer14;
	eb = Ceb;
	Esec_n = CEsec_n;
	Epl_n = CEpl_n;
	delta_f_n = Cdelta_f_n;
	delta_e_n = Cdelta_e_n;
	epl_n = Cepl_n;
	fnew_n = Cfnew_n;
	Enew_n = CEnew_n;
	ere_n = Cere_n;
	fnew_str_n = Cfnew_str_n;
	Enew_str_n = CEnew_str_n;
	ere_str_n = Cere_str_n;
	f71target = Cf71target;
	E71target = CE71target;
	e72target = Ce72target;
	f3target = Cf3target;
	f81target = Cf81target;
	E81target = CE81target;
	e82target = Ce82target;
	f4target = Cf4target;
	delta_e_p = Cdelta_e_p;
	delta_f_p = Cdelta_f_p;
	fnew_p = Cfnew_p;
	Enew_p = CEnew_p;
	ere_p = Cere_p;
	fnew_str_p = Cfnew_str_p;
	Enew_str_p = CEnew_str_p;
	ere_str_p = Cere_str_p;
	epl_p = Cepl_p;
	Epl_p = CEpl_p;
	er2 = Cer2;
	er4 = Cer4;
	er9 = Cer9;
	Esec_p = CEsec_p;
	er10 = Cer10;
	ea = Cea;
}

void
changManderConcrete01::commitStateVar(void){
	CisCracked = isCracked;
	CisSpalled = isSpalled;
	Crule = Trule;
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;
	Ceo = eo;
	Cer1 = er1;
	Cer3 = er3;
	Cer6 = er6;
	Cer13 = er13;
	Cer14 = er14;
	Ceb = eb;
	CEsec_n = Esec_n;
	CEpl_n = Epl_n;
	Cdelta_f_n = delta_f_n;
	Cdelta_e_n = delta_e_n;
	Cepl_n = epl_n;
	Cfnew_n = fnew_n;
	CEnew_n = Enew_n;
	Cere_n = ere_n;
	Cfnew_str_n = fnew_str_n;
	CEnew_str_n = Enew_str_n;
	Cere_str_n = ere_str_n;
	Cf71target = f71target;
	CE71target = E71target;
	Ce72target = e72target;
	Cf3target = f3target;
	Cf81target = f81target;
	CE81target = E81target;
	Ce82target = e82target;
	Cf4target = f4target;
	Cdelta_e_p = delta_e_p;
	Cdelta_f_p = delta_f_p;
	Cfnew_p = fnew_p;
	CEnew_p = Enew_p;
	Cere_p = ere_p;
	Cfnew_str_p = fnew_str_p;
	CEnew_str_p = Enew_str_p;
	Cere_str_p = ere_str_p;
	Cepl_p = epl_p;
	CEpl_p = Epl_p;
	Cer2 = er2;
	Cer4 = er4;
	Cer9 = er9;
	CEsec_p = Esec_p;
	Cer10 = er10;
	Cea = ea;
}
