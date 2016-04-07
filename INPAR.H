// inpar.h      new porous

#ifndef _INPAR_LIBSAM

#define _INPAR_LIBSAM

// Definition of auxiliary structures
#include <stdio.h>
#include "matr.h"

class MaterialConstants{
 public:
  //flags showing whether to calculate these features
  int damage_flag;
  int active_plasticity_flag;
  int microplasticity_flag;

  int kzer;
  int KVar; // numb. of transformation variants
  int KSlipTypes; //numb. of slip planes types
  int KSlipPlanesOfAllTypes;
  int KTauY;//numb. of points to specify TauYieldEq(iSlipType,Temper)
  int* KSlip;//[KSlipTypes]; //numb. of slip planes of given type
  Real TauYieldEq(int iSlipType, Real Temper);// function to calculate
			     //equilibrium yield limits for slip plane
			     // of iSlipSystem type and temperature Temper
  Real* HSlip;//[KSlipTypes];// hardening modules for slip
			     //in a plane of given type
  Real* R0Slip;//[KSlipTypes];// recovery coef. for TauDef --
			    // deformation part of flow stress
  Real* USlip;//[KSlipTypes]; // activation energy
  Real* PSlip;//[KSlipTypes]; // parameter for TauY rate dependence
  Real* AlfaSlip;//[KSlipTypes]; // power for TauY rate dependence

  Real* TemperTable;//[KTauY]=   (K)
  Real** TauYieldTable;//[KTauY][KSlipTypes]=  (MPa)

// constants derivative from others:
  Real Gm;
  Real Ftr;
// constants to be read from datafile:
  Real Ms;
  Real Mf;
  Real As;
  Real Af;
  Real T0;
  Real Q0;
  Real ScalingFactor;
  Matr33 Dtrans;
  Matr33 Dev_Dtrans;
  Real Gam_D;//intensity of Dev_Dtrans
  Real Spur_D;
  Real YoungA;
  Real YoungM;
  Real PuasA;
  Real PuasM;
  Real TECoefA;
  Real TECoefM;
  Real HeatConductivityA;
  Real HeatConductivityM;
  Real SpecificHeatA;
  Real SpecificHeatM;
  Real DensityA;
  Real DensityM;
  Real FYield;//Yield limit of resistance force for microplastic flow
  Real FYieldTwinning;//Yield limit of force for twinning
  Real Kmp;//fraction of plasic microdeformation
  Real Kap2002;//scaling factor for microplastic deformation
  Real Dolya; // коэффициент для FYield (0<=Dolya<=1)
  Real Rm; // Rm = -log(1-Dolya) priblizhenno= Dolya
  Real InitialOrder;
  Real CorrectRoGr;
  Real Tau_c; // normalizer for TauY_def
  Real Decrement_A; //internal friction in austenite (at Fi==0)
  Real Decrement_M; //internal friction in martensite (at Fi==1)
  Real Decrement_T; //internal friction in two-phase state (at Fi==0.5)
  
// Formulae for damage
  // db/dt = dbe/dt - (1/be*)(dbe/dt)*abs(b)*H(b*dbe/dt)
  // df/dt = (q2+q3*f)*abs(b)^alfa*abs(dbe/dt) - r(T)*f
  // where r(T) = r0*exp(-U/kT)
  //Condition of micro failure:
  // abs(b)*f = b_F, b_F being the critical value
// Constants for damage
  Real Al_d; //alfa
  Real be_zv; //be*
  //Real q2_F;
  Real q3_F;
  Real r0_F;//r0
  Real U_F; //U
  //Real b_F; //critical value for fracture connected with variant
  Real Tau_F; //critical stress for fracture connected with variant
  Real A_F; //parameter for damage criterion relating to f-defects
  Real K1_F; //parameter for damage criterion relating to SpSig
  Real K2_F; //parameter for damage criterion relating to Tau
  Real Pi_cr; // critical value for grain damage
  int kzer_cr; //critical numb of destroyed grains for macro-fracture

//Constructor reading constants from file
  MaterialConstants(void){}
  void Init(void);
//Destructor
  ~MaterialConstants(void);
}; //end of class MaterialConstants


class VariantsOrientations{
   Matr33 *Orvvar_mp;
   Matr33** OrSlip_mp;

  public:
   int kvar;

// Constructor allocating memory and reading 24 cubic
// matrixes from file "kubmatr.d" on disk
   VariantsOrientations(void){}
   void Init(MaterialConstants& MC);

   Matr33 *Dev_var;// deporoustors

// Destructor
  ~VariantsOrientations(void)
  {
    delete[]Orvvar_mp;
    delete[]OrSlip_mp;
    delete[]Dev_var;
  }
// Getting matrix with specified ivar
   Matr33 VarOr(int ivar){ return Orvvar_mp[ivar]; }
   Matr33 SlipOr(int iSlipType, int iPlane)
   {
      return(OrSlip_mp[iSlipType][iPlane]);
   }
   Matr33 operator[](int ivar){ return(Orvvar_mp[ivar]); }
}; // end of class Variants Orientations



// Definition of class GrainsParameters
class GrainsParameters
{
public:
    char GrainsOrientationsDistribution[16];
    int k_zer;
    Matr33 *Orvzer_mp;

//Constructor
    GrainsParameters(void){}
    int Init();
    //void InitMono(Real Alfa, Real Dist, Real Dist110);
//Destructor
   ~GrainsParameters(void) { delete[]Orvzer_mp; }
//Getting izer-th grain orientation
    Matr33 ZerOr(int izer) { return Orvzer_mp[izer]; }
}; //end of class GrainsParameters

// Definition of class InternalParameters
class InternalParameters
{
  public:
    int kloc;
    int kzer;
    int kzer_dmg; // number of destroyed grains
    int kslip;
    Real *Floc_dp;
    Real *mpDef_dp;
    Real *Fyloc_dp;
    Real *TYDef_dp; // deformational part of TauYield
    Real *Gb_dp;    //total slip (Odquist parameter) on a plane
    Real *Db_dp;    //rate of the total slip on a plane

    //variables for damage created by microplastic deformation
    Real *B_dmg;//[kloc] defects reversible by   --"--       --"--
    Real *F_dmg;//[kloc] defects irreversible by
    Real *pi_dmg;//[kloc] damage associated with each variant

	//variables for damage created by athermal plastic deformation
	Real *B1_pl_dmg; //[iPlaneGlob]  component of vector B_pl - defects reversible by plastic strain
	Real *B2_pl_dmg; //[iPlaneGlob]  component of vector B_pl - defects reversible by plastic strain
	Real Bmean_pl_dmg; // Mean on all grains absolute value of B1 and B2
	Real Bmax_pl_dmg; // Maximum on all grains absolute value of B1 and B2
    Real *F_pl_dmg; //[iPlaneGlob] increment of defects irreversible by plastic strain
	Real Fmean_pl_dmg; // Mean on all grains value of irreversible defects
	Real Fmax_pl_dmg; // Maximum value of irreversible defects
    Real *Pi1_pl_dmg; //[iPlaneGlob] component of vector Pi - damage due to plastic deformation
    Real *Pi2_pl_dmg; //[iPlaneGlob] component of vector Pi - damage due to plastic deformation

    //intergranular stress
    Matr33 *RoGrain; //[kzer]
    Matr33 *hRoGrain;
    Real RoGrainNorma; // Maximum norm of all RoGrain[izer]

    //damage of grains
    Real* PI_dmg; //[kzer] PI=0.0 - grain is undamaged, PI>0 - fracture of grain
    Real* bf_dmg; //[kzer] max(on variants) product of b_dmg*f_dmg
    //macro-fracture
    int PI_macro;// 0/1 - no/yes macro-fracture
    int izer_bf_max;


    Real Phase;
    Real dPhase;
    Matr33 *EpsUnGr;
    Matr33 PhaseStrain;
    Matr33 mpStrain;
    Matr33 TEStrain;
    Matr33 ElStrain;
    Matr33 PlasticStrain;
    Matr33 UnelasticStrain;
    Matr33 TotalStrain;
    Matr33 dTotalStrain;
    Real Order;

// Default constructor
    InternalParameters(void){}
// Constructor - unwritten
//  InternalParameters(MaterialConstants&,GrainsParameters&);
// Initializer:
    void Init(MaterialConstants&,GrainsParameters&);
//Destructor
   ~InternalParameters(void);

//Assignment
   void operator=(InternalParameters&);
}; // class InternalParameters

class Circular_beam //circular beam is a part of a slotted spring
{
	public:
	int i; // number of the beam in a slotted spring
	real thickness;
	real width;
	real length; 
	real column;
	real little_length;
	real band;
	int noV;
	real J_r;//inertia moment for bending
	real area; 
	InternalParameters XXold;
	InternalParameters XXnew;
	real K_sig; //constant to link force and stress
	real K_sig2;
	real hdisplac;
	real displac_old;
	real displac_new;
	real heps_elastic;
	real eps_elastic_old;
	real eps_elastic_new;
	real Check_destr_old;
	real Check_destr_new;
	//Default constructor
	//Circular_beam(void){}
	real heps_beam_old; //increment of strain in the beam
	real heps_beam_new;
	real e33_fix;
	real s33_fix;
	real e22_el_fix;
	real e22_mp_fix;
	real s22_fix;
	real force_fix;//force caused beams' connection
	real L_macro_fix;
	real sig_to_control;//to see s22_1 in .dat file
	real sigma33;//to see s33 in every beam in .dat file 
	real phoenix;//flag; =0 if the beams weren't connected, =1 if the beams were connected
	real loading_type;//flag to determine loading type (isothermal or non isothermal regime) to choose right value for sig(1,1) after beam's connection
	real L_macro_old;//length of working element in vertical direction 
	real L_macro_new;
	//Initializer
	void Init(int iuz, MaterialConstants& MC,GrainsParameters& GrPar);
	//Update XXold=XXnew
	void Update(void);
	//~Circular_beam(void){} //Default destructor
	//void destroy(void);//Procedure to release memory of XXold and XXnew
	real displ_old(void);
	real displ_new(void);
};


class Porous_exVIA
{//Vibroisolation element with two round slotted springs (RSS) and arreter
	public:
	int kc_beams; // number of circular beams
	Circular_beam* c_beam;
	//Circular_beam* c_beam2;
	Circular_beam* c_beam_max;//pointer to the circular beam with the maximum number (kc_beams-1) for the 1st RSS
	//Circular_beam* c_beam2_max;//pointer to the circular beam with the maximum number (kc_beams-1) for the 2nd RSS
	//int macro_failure1;
	//int macro_failure2;
	int macro_failure;

	Real mass;//mass to be isolated from external vibrations
	Real external_damping;

	real tem_old;//temperature of the 1-st RSS
	real tem_new;
	real htem;
	
	//real tem2_old;//temperature of the 2-nd RSS
	//real tem2_new;
	//real htem2;
	
	real force_old;// force on 1-st RSS, not on the shaft. The force on the shaft = (-1)*force1 
	real force_new;
	real hforce;
	real force_new_eff;

	//real force2_old;// force on 2-nd RSS,
	//real force2_new;
	//real hforce2;

	real force_reaction_old;//reaction force in the displacement limiter
	real force_reaction_new;
	real hforce_reaction;
	real displ0_old;// initial displacement of 1-st RSS
	real displ0_new;
	//real displ02_old;// initial displacement of 2-nd RSS
	//real displ02_new;
	
	real displ_old; //displacement of mass (attached to the housing)
	real displ_new;
	real hdispl;
	real displ_old_1;
	real displ_new_1;
	real displ_old_2;
	real displ_new_2;
	real displ_old_3;
	real displ_new_3;
	real displ_old_4;
	real displ_new_4; 
	real displ_old_5;
	real displ_new_5;
	real displ_old_6;
	real displ_new_6;

	real displ_new_total;
	real hdispl_1;
	real hdispl_2;
	real hdispl_3;
	real hdispl_4;
	real hdispl_5;
	real hdispl_6;

	real e_macro;
	real e_macro_1;
	real e_macro_2;
	real e_macro_elastic;
	real de_macro; //deformation increment on macroscale
	
	real check1;
	real check2;
	real check3;
	real check4;
	real check5;
	real k_parallel;
	real porosity;
    real stress_macro;
	real area_macro;

	real K_u_1;
	real K_u_2;
	real konst;

	real Check_destr_old;
	real Check_destr_new;
	
	real L_macro_new_overall;
	real L_macro_old_overall;

	real velocity_old;
	real velocity_new;
	real hvelocity;


//Magnetic force 
// F_mag = mf_factor*( Fa1*exp((u-u1)/xm1) - Fa2*exp((u-u2)/xm2) )
//Fa1=23.0g ,u1=-2.0mm ,xm1=2.0mm ,Fa2=167.0g ,u2=0.0mm ,xm2=17.0mm
	Real mf_factor_old;//for gradual adding of the magnetic force when mf_factor grows from 0.0 to 1.0
	Real mf_factor_new;
	Real hmf_factor;
	
	Real Fa1;
	Real u1;
	Real xm1;
	Real Fa2;
	Real u2;
	Real xm2;
//Magnetic_force_scale_1= 0.23 //Fa1 in Newtons
//Magnetic_force_shift_1= -2.0 //u1 in mm
//Magnetic_force_dist_1= 2.0 //xm1 in mm
//
//Magnetic_force_scale_2= 1.67 //Fa2 in Newtons
//Magnetic_force_shift_2= 0.0 //u2 in mm
//Magnetic_force_dist_2= 17.0 //xm2 in mm

	real displ_limit; //maximum displacement (limited by a detent)
	real F_mag(real displ, real factor);//The magnetic force for a given displacement
	//Default constructor
	Porous_exVIA(void){}
	//Initializer
	void Init(Real tem, Real force, MaterialConstants& MC ,GrainsParameters& GrPar);
	//Update XXold=XXnew
	void Update(void);
	//Destructor
	~Porous_exVIA(void);
};

class Isolation
{//meaning a helical SMA wire isolator
  public:
	Real	d1; //wire diameter
	Real	D1; // coil diameter
	Real	n1; //number of convolutions (coils)
	Real	d2; //wire diameter
	Real	D2; // coil diameter
	Real	n2; //number of convolutions (coils)
	Real	m; // isolated mass
	Real	K_F_tau1;//coef. for getting axial force applied to a helix from tangential stress
	Real	K_gam_u1;//coef. for getting shear from elongation of a helix
	Real	K_F_tau2;//coef. for getting axial force applied to a helix from tangential stress
	Real	K_gam_u2;//coef. for getting shear from elongation of a helix
	Real	external_damping; // damping coef.(Nu2) for external friction 
	//Default constructor
	Isolation(void){}
	//Initializer
	void Init(void);
	//Destructor
	~Isolation(void){}
};

class Majatnik{
 public:
  Real       Length,
           Diameter,
             Plecho,
         MassOfLoad,
       LengthOfLoad,
              DecrA,
              DecrM,
              DecrT,
                 Jp;
// Elastic unit
  char Pendilum[16];
  Real      Length2,
          Diameter2,
             Young2,
              Puas2,
                 K2;
  int  N;
  Real   *C,
	 *omega;
//Default constructor
  Majatnik(void){}
//Destructor
  ~Majatnik(void){}
//Initializer
  void Init(void);
  Real Psi(Real);
  Real dPsi(Real);
};


#endif