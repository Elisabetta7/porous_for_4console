#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defglob.h"
#include "diaproc.h"
#include "matr.h"
#include "inpar.h"
#include <conio.h>                                                                 
#include <iostream>                                                                 
#include <windows.h>
#define MALBUF 16
#define MAXBUF 128


void Porous_exVIA::Init(Real tem , Real force ,MaterialConstants& MC ,GrainsParameters& GrPar)   //добавлен real force
{
	FILE *filrd;
	//printf("\n\nVibro_isol_arret::Init() started.");
	if( (filrd=fopen("isolation.d","r"))==0){
		fprintf(stderr,"\n\nCannot open isolation.d");
		fprintf(stderr,"\n\nProgram will be terminated\n");
		pause(); exit(1);
	}
	fkeyread_int(filrd, "number_of_beams", &kc_beams);
	fclose(filrd);

	c_beam = new Circular_beam [kc_beams];
	//c_beam2 = new Circular_beam [kc_beams];
	for(int i=0;i<kc_beams;i++){
		c_beam[i].Init(i,MC,GrPar);
		//c_beam2[i].Init(i,MC,GrPar);
	 }//for(i)	
	c_beam_max = c_beam + (kc_beams-1);//pointer to the circular beam with the max number for the 1st RSS
	//c_beam2_max = c_beam2 + (kc_beams-1);//pointer to the circular beam with the max number for the 2nd RSS
	
	filrd=fopen("isolation.d","r");

	//fkeyread_double(filrd, "Magnetic_force_scale_1",&Fa1); //already in Newtons
	//fkeyread_double(filrd, "Magnetic_force_shift_1",&u1);		u1*=1e-3; // mm --> m
	//fkeyread_double(filrd, "Magnetic_force_dist_1",&xm1);		xm1*=1e-3; // mm --> m
	//
	//fkeyread_double(filrd, "Magnetic_force_scale_2",&Fa2); //already in Newtons
	//fkeyread_double(filrd, "Magnetic_force_shift_2",&u2);		u2*=1e-3; // mm --> m
	//fkeyread_double(filrd, "Magnetic_force_dist_2",&xm2);		xm2*=1e-3; // mm --> m

	fkeyread_double(filrd, "Displacement_limit",&displ_limit);	displ_limit*=1e-3; // mm --> m
	fkeyread_double(filrd, "isolated_mass", &mass);// in kilograms
	fkeyread_double(filrd, "damping_coef_of_external_friction", &external_damping);//in N/(m*s) == kg/s
	fkeyread_double(filrd, "porosity", &porosity);
	fkeyread_double(filrd, "parallel_elements", &k_parallel);//number of parallel chains of bulks

	fclose(filrd);

	//initial values for temperatures - they must be assigned after reading of [ini_section] of wzd-file
	tem_old=tem;// temperature of the 1-st RSS
	tem_new=tem;
	

	//zero initial values for forces and displacements:
	force_old=0.0;// on 1-st RSS. The force on the shaft = (-1)*force1 
	force_new=0.0;
	hforce=0.0;
	force_new_eff=0.0;
	
	
	force_reaction_old=0.0;//reaction force in the displacement limiter
	force_reaction_new=0.0;
	hforce_reaction=0.0;
	
	displ0_old=0.0;// initial displacement of 1-st RSS
	displ0_new=0.0;

	displ_old=0.0; //displacement - addition to displacement after the assembly (common for both RSS)
	displ_new=0.0; 
	hdispl=0.0;
	hdispl_1=0.0;
	hdispl_2=0.0;
	displ_new_1=0.0;
	displ_new_2=0.0;
	displ_old_1=0.0;
	displ_old_2=0.0;
	hdispl_3 = 0.0;
	hdispl_4 = 0.0;
	displ_new_3 = 0.0;
	displ_new_4 = 0.0;
	displ_old_3 = 0.0;
	displ_old_4 = 0.0;
	hdispl_5 = 0.0;
	displ_new_5 = 0.0;
	displ_old_5 = 0.0;
	hdispl_6 = 0.0;
	displ_new_6 = 0.0;
	displ_old_6 = 0.0;
    displ_new_total=0.0;
	e_macro=0.0;
	e_macro_elastic = 0.0;
	stress_macro=0.0;

	
	velocity_old=0.0;
	velocity_new=0.0;
	hvelocity=0.0;


// F_mag = mf_factor*( Fa1*exp((u-u1)/xm1) - Fa2*exp((u-u2)/xm2) )
	//mf_factor_old=0.0;//for gradual adding of the magnetic force when mf_factor grows from 0.0 to 1.0
	//mf_factor_new=0.0;
	//hmf_factor=0.0;
	//
	/*macro_failure1=0;
	macro_failure2=0;*/
	macro_failure=0;
	//length_overall =  3.5* c_beam[0].Length + 1.3* c_beam[1].Length +  0.35* c_beam[2].Length + 3.5 * c_beam[3].Length + 2* c_beam[4].Length;

	printf("\nVibro_isol_arret.Init:");
	//printf("iuz=%i thickness=%lg width=%lg radius=%lg",iuz,thickness ,width ,radius);
}

void Porous_exVIA::Update(){
	macro_failure=0;
	for(int i=0;i<kc_beams;i++){
		if(c_beam[i].XXnew.PI_macro == 1 ) macro_failure=1;
		//if(c_beam2[i].XXnew.PI_macro == 1 ) macro_failure2=1;
		c_beam[i].Update();
		c_beam[i].Check_destr_old = c_beam[i].Check_destr_new;
		c_beam[i].heps_beam_old = c_beam[i].heps_beam_new;
		c_beam[i].displac_old = c_beam[i].displac_new;
		c_beam[i].eps_elastic_old = c_beam[i].eps_elastic_new;
		c_beam[i].L_macro_old = c_beam[i].L_macro_new;
		//c_beam2[i].Update();
	}//for(i)	
	tem_old = tem_new;
	//tem2_old = tem2_new;
	force_old = force_new;
	//force2_old = force2_new;
	displ0_old = displ0_new;
	//displ02_old = displ02_new;
	displ_old_1 = displ_new_1;
	displ_old_2 = displ_new_2;
	displ_old = displ_new;
	//velocity_old=velocity_new;
	//mf_factor_old = mf_factor_new;
	macro_failure = macro_failure;
	force_reaction_old = force_reaction_new;
	//fi_angle_old = fi_angle_new;
	//radius_old =  radius_new;
}

Porous_exVIA::~Porous_exVIA()
{
	for(int i=0;i<kc_beams;i++){
		(c_beam+i)->Circular_beam::~Circular_beam();
		//(c_beam2+i)->Circular_beam::~Circular_beam();
	}//for(i)	
	delete[] c_beam;
	//delete[] c_beam2;
}

//class Circular_beam //circular beam is a part of a slotted spring
//{
//	public:
//	int i; // number of the beam in a slotted spring
//	real thickness;
//	real width;
//	real radius;
//	real J_r;//inertia moment for bending
//	real area;
//	//real mass; 
//	//real external_damping; 
//	InternalParameters XXold;
//	InternalParameters XXnew;
//	real K_sig;
//	real K_u;
//	//Default constructor
//	Circular_beam(void){}
//	//Initializer
//	void Init(MaterialConstants& ,GrainsParameters&);
//	void Update(void);
//	//Destructor
//	~Circular_beam(void){}
//	real displ(void){}
//};

void Circular_beam::Init(int iuz, MaterialConstants& MC ,GrainsParameters& GP) //circular beam is a part of a slotted spring
{
	FILE *filrd;
	//printf("\n\nCircular_beam::Init :");
	if( (filrd=fopen("isolation.d","r"))==0){
		fprintf(stderr,"\n\nCannot open isolation.d");
		fprintf(stderr,"\n\nProgram will be terminated\n");
		pause(); exit(1);
	}

	ffindkey(filrd, "beam_widths");
	for(int j=0;j<=iuz;j++){fscanf(filrd,"%lg",&width);} width*=1e-3; //mm --> m
	//we ignore all values until its number == iuz
	
	ffindkey(filrd, "beam_lengths");
	for(int j=0;j<=iuz;j++){fscanf(filrd,"%lg",&length);} length*=1e-3; //mm --> m

	ffindkey(filrd, "beam_thicknesses");
	for(int j=0;j<=iuz;j++){fscanf(filrd,"%lg",&thickness);} thickness*=1e-3; //mm --> m

	ffindkey(filrd, "beam_columns");
	for(int j=0;j<=iuz;j++){fscanf(filrd,"%lg",&column);} column*= 1e-3;

	ffindkey(filrd, "beam_little_lengths");
	for(int j=0;j<=iuz;j++){fscanf(filrd,"%lg",&little_length);} little_length*=1e-3;

	ffindkey(filrd, "beam_n");
	for (int j = 0; j <= iuz; j++){ fscanf(filrd, "%lg", &noV); } noV *= 1;
	//std::cout<<"fi_new="<<fi_angle;  //если не комментировать эти строки, 
    //system("pause");                      //то появится окно с новой фи и всё 
    //std::cin.get();
	//printf("\n iuz=%u, angle_old=%lg , angle_new=%lg", iuz   , fi_angle_old ,fi_angle_new);
	//pause();

	XXold.TotalStrain(2,2)=0.0;
	XXnew.TotalStrain(2,2)=0.0;
	XXold.TotalStrain(1,1) = 0.0;
	XXnew.TotalStrain(1, 1) = 0.0;
	e22_el_fix = 0.0;
	e22_mp_fix = 0.0;
	e33_fix = 0.0;
	s22_fix = 0.0;
	s33_fix = 0.0;
	force_fix = 0.0;

	Check_destr_new = 0.0;
	Check_destr_old = 0.0;
	phoenix = 0.0;
	loading_type = 0.0;
	sig_to_control = 0.0;
	sigma33 = 0.0;

	displac_old = 0.0;
	displac_new = 0.0;
	hdisplac = 0.0;
	heps_elastic = 0.0;
	eps_elastic_old = 0.0;
	eps_elastic_new = 0.0;

	L_macro_old = thickness+column;
	L_macro_new = 0.0;

	XXold.Init(MC,GP);
	XXnew.Init(MC,GP);
	fclose(filrd);
	//std::cout<<"fi_new="<<fi_angle;  // для проверки, читает ли.
    //system("pause");                       
    //std::cin.get();

	J_r = width*cube(thickness)/12.0;
	area = thickness*width;
	//area_macro_0 = pi*(sqr(diagonal + dia));
	//area_macro_2 = (thickness+0.00003)*(width+0.00003);//0.005 and 0.005 - pore's diametres, 0.5661 = (1-p)^2 
	//K_sig = ((-2*sin(fi_angle_z/2))/(fi_angle_z*width*thickness)
	//          +6*sqr(radius_z)*(-2*(sin(fi_angle_z/2))/fi_angle_z+1)/(width*sqr(thickness)*(radius_z+thickness/2))); //это мой
	K_sig = 3 * (length - little_length / 2) / (2*width*thickness*thickness);
	K_sig2 = 1/(little_length*width); 
	band = thickness + column;
	printf("\nCircular_beam.Init:");
	//pause();
	//printf("\nisolated mass=%lg external damping=%lg",mass,external_damping);
};

void Circular_beam::Update(void) { XXold = XXnew; }


void Isolation::Init(void)
{
	FILE *filrd;
	//printf("\n\nIsolation_INIT:");
	if( (filrd=fopen("isolation.d","r"))==0)
	{
	fprintf(stderr,"\n\nCannot open isolation.d");
	fprintf(stderr,"\n\nProgram will be terminated\n");
	pause(); exit(1);
	}
	
	fkeyread_double(filrd, "wire_diameter_1", &d1); d1*=1e-3; //mm --> m
	fkeyread_double(filrd, "coil_diameter_1", &D1); D1*=1e-3; //mm --> m
	fkeyread_double(filrd, "number_of_coils_1", &n1); 
	fkeyread_double(filrd, "wire_diameter_2", &d2); d2*=1e-3; //mm --> m
	fkeyread_double(filrd, "coil_diameter_2", &D2); D2*=1e-3; //mm --> m
	fkeyread_double(filrd, "number_of_coils_2", &n2); 
	fkeyread_double(filrd, "isolated_mass", &m); //already in kg
	fkeyread_double(filrd, "damping_coef_of_external_friction", &external_damping); 

	K_F_tau1 = pi*d1*d1*d1/(8.0*D1);//coef. for getting axial force applied to a helix from tangential stress
	K_gam_u1 = d1/(pi*D1*D1*n1);//coef. for getting shear from elongation of a helix
	K_F_tau2 = pi*d2*d2*d2/(8.0*D2);//coef. for getting axial force applied to a helix from tangential stress
	K_gam_u2 = d2/(pi*D2*D2*n2);//coef. for getting shear from elongation of a helix

	printf("\n\nIsolation.Init(): 1-st helical element:");
	printf("\nnumber of coils n1=%5.1lf, wire_diameter d1=%6.2lf mm, coil_diam D1=%5.1lf mm",n1,d1*1e3,D1*1e3);
	printf("\nK_F_tau1=%8.6lg, K_gam_u1=%8.6lg kg." ,K_F_tau1 ,K_gam_u1);
	printf("\n\n                  2-nd helical element:");
	printf("\nnumber of coils n2=%5.1lf, wire_diameter d2=%6.2lf mm, coil_diam D2=%5.1lf mm",n2,d2*1e3,D2*1e3);
	printf("\nK_F_tau2=%8.6lg, K_gam_u2=%8.6lg kg." ,K_F_tau2 ,K_gam_u2);
	printf("\n, isolated mass=%8.3lg kg. external_damping=%lg",m,external_damping);
}// end of Isolation::Init()


void Majatnik::Init(void)
{
 FILE *filrd;

 printf("\n\nMajatnik_INIT:");
 if( (filrd=fopen("majatnik.d","r"))==0)
 {
  fprintf(stderr,"Cannot open majatnik.d \n");
  pause(); exit(1);
 }
 frd1d(filrd,&Length); Length*=1e-3;  //from mm into m
 frd1d(filrd,&Diameter); Diameter*=1e-3;
 frd1d(filrd,&Plecho); Plecho*=1e-3;
 frd1d(filrd,&MassOfLoad); MassOfLoad*=1e-3; //from g into kg
 frd1d(filrd,&LengthOfLoad); LengthOfLoad*=1e-3; //from mm into m
 frd1d(filrd,&DecrA);
 frd1d(filrd,&DecrM);
 frd1d(filrd,&DecrT);

 frd1s(filrd,Pendilum);
 switch(Pendilum[0])
 {
  case 'S': K2=0.0;
            break;

  case 'D': frd1d(filrd,&Length2); Length2*=1e-3;  //from mm into m
            frd1d(filrd,&Diameter2); Diameter2*=1e-3;
            frd1d(filrd,&Puas2);
            frd1d(filrd,&Young2);Young2*=1E9;
	    if(Length2<=0.0 || Diameter2<=0.0 || Puas2<=0.0 || Young2<=0.0)
            {
             fprintf(stderr," Error in the steel element parameters definition\n");
             pause(); exit(1);
            }
            K2=pi*Diameter2*Diameter2*Diameter2*Diameter2*Young2/64/(1+Puas2)/Length2;
            break;
  default:  fprintf(stderr,"Error Pendilum parameter '%s' is bad",Pendilum);
            pause(); exit(1);
            break;
 }

 fclose(filrd);

 N=0;
 C=0; omega=0;
 Jp=2*MassOfLoad*Plecho*Plecho
     *(1+LengthOfLoad/Plecho+(LengthOfLoad*LengthOfLoad)/(3*Plecho*Plecho));

 printf("\nLength=%6.1lg mm, Diameter=%6.1lg mm,",Length*1e3,Diameter*1e3);
 printf("\nPlecho=%6.1lg mm, MassOfLoad=%6.1lg g,",Plecho*1e3,MassOfLoad*1e3);
 printf("\nJp=%6.1lg g*mm^2,",Jp*1e9);
 printf("\nDecrA=%8.5lg, DecrM=%8.5lg, DecrT=%8.5lg,",DecrA,DecrM,DecrT);
 printf("\nPendilum='%s'",Pendilum);
 if(Pendilum[0]=='D')
 {
  printf("\nLength2=%6.1lg mm, Diameter2=%6.1lg mm,",Length2*1e3,Diameter2*1e3);
  printf("\nPuas2=%6.1lg, Young2=%6.1lg GPa.",Puas2,Young2*1e-9);
  printf("\nK2=%6.1lg.",K2);
 }
}

Real Majatnik::Psi(Real t)
{
 Real Psi=0;
 for( int i=0;i<N;i++) Psi+=C[i]*sin(omega[i]*t);
 return(Psi);
}

Real Majatnik::dPsi(Real t)
{
 Real dPsi=0;
 for( int i=0;i<N;i++) dPsi+=C[i]*omega[i]*cos(omega[i]*t);
 return(dPsi);
}

void MaterialConstants::Init(void)
{
  FILE *fc;
  Real Dx;
  int i,j,Kx;
  if( (fc=fopen("cnst_key.d","r"))==0 ){ //почему-то не понимает NULL при запуске
    fprintf(stderr,"Cannot open cnst_key.d \n");
    pause(); exit(1);
  }
  //pause("file cnst_key.d successfully opened");
  //fnameread_int(fc,damage_flag);
  //printf("\inpar: damage_flag=%i",damage_flag);
  //pause();
  fkeyread_int(fc, "damage_flag", &damage_flag);
  //printf("\ndamage_flag=%1i",damage_flag);
  //pause();
  
  fkeyread_int(fc,"active_plasticity_flag", &active_plasticity_flag);
  fkeyread_int(fc,"microplasticity_flag", &microplasticity_flag);
  printf("\n InPar: damage_flag=%i active_plasticity_flag=%i microplasticity_flag=%i"
                     ,damage_flag, active_plasticity_flag, microplasticity_flag);
  //pause();

  ffindkey(fc, "Mf, Ms, As, Af, q0");
  fscanf(fc,"%lg%lg%lg%lg%lg",&Mf, &Ms, &As, &Af, &Q0);
  Q0*=1e6; // from MJ/m^3 into J/m^3

  Spur_D=0.0;
 
  ffindkey(fc, "Matr Dtrans");
  for(i=0;i<3;i++) for(j=0;j<3;j++){
   fscanf(fc,"%lg",&Dx); Dtrans.C(i,j)=Dx;
  }

  Spur_D=0.0;
  for(i=0;i<3;i++){
   Spur_D += Dtrans.C(i,i);
  }

  Dev_Dtrans = Dtrans;
  for(i=0;i<3;i++){
   Dev_Dtrans.C(i,i) -= Spur_D/3.0;
  }

  Gam_D = Dev_Dtrans.EDevInt();

  fkeyread_int(fc, "k_zer", &kzer); //number of grains
  fkeyread_int(fc, "KVar", &KVar); //number of transformation variants
  fkeyread_int(fc, "KSlipTypes", &KSlipTypes); //numb. of slip planes types

  if( !(KSlip = new int [KSlipTypes] )){
    fprintf(stderr," Cannot allocate KSlip\n");
    pause(); exit(1);
  }
  if( !(HSlip = new Real [KSlipTypes] )){
    fprintf(stderr," Cannot allocate HSlip\n");
    pause(); exit(1);
  }
  if( !(R0Slip = new Real [KSlipTypes] )){
    fprintf(stderr," Cannot allocate R0Slip\n");
    pause(); exit(1);
  }
  if( !(USlip = new Real [KSlipTypes] )){
    fprintf(stderr," Cannot allocate USlip\n");
    pause(); exit(1);
  }
  if( !(PSlip = new Real [KSlipTypes] )){
    fprintf(stderr," Cannot allocate PSlip\n");
    pause(); exit(1);
  }
  if( !(AlfaSlip = new Real [KSlipTypes] )){
    fprintf(stderr," Cannot allocate AlfaSlip\n");
    pause(); exit(1);
  }

  ffindkey(fc, "KSlip[KSlipType]");
  KSlipPlanesOfAllTypes=0;
  
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%u", &Kx); KSlip[i]=Kx;
   KSlipPlanesOfAllTypes += Kx;
  }

  ffindkey(fc, "HSlip[KSlipType]");
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%lg", &Dx); HSlip[i]=Dx*1e6;
  }

  ffindkey(fc, "R0Slip[KSlipType]");
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%lg", &Dx); R0Slip[i]=Dx;
  }

  ffindkey(fc, "USlip[KSlipType]");
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%lg", &Dx); USlip[i]=Dx * 1e3; // from kJ into J
  }

  ffindkey(fc, "PSlip[KSlipType]");
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%lg", &Dx); PSlip[i]=Dx*1e6; // in MPa*s^AlfaSlip
  }

  ffindkey(fc, "AlfaSlip[KSlipType]");
  for(i=0; i<KSlipTypes; i++){
   fscanf(fc, "%lg", &Dx); AlfaSlip[i]=Dx;
  }

  fkeyread_int(fc, "KTauY", &KTauY); //numb. of points to specify TauYieldEq
   TemperTable = new Real[KTauY];
   TauYieldTable = new Real*[KSlipTypes];
   for(i=0;i<KSlipTypes;i++)
     TauYieldTable[i] = new Real [KTauY];

  ffindkey(fc, "TemperTable[KTauY]");
   for(i=0;i<KTauY;i++)
     fscanf(fc,"%lg", TemperTable+i);

  ffindkey(fc, "TauYieldTable[KSlipTypes][KTauY]");
   for(j=0;j<KSlipTypes;j++)
     for(i=0;i<KTauY;i++){
       fscanf(fc,"%lg", &Dx);
       TauYieldTable[j][i]=Dx*1e6;//from MPa to Pa
       if(!active_plasticity_flag) TauYieldTable[j][i]=1e100;
         //setting big values to forbid active plasticity
       printf("\nTauYield[%u][%u]=%g",j,i,TauYieldTable[j][i] );
     }

  fkeyread_double(fc, "InitialOrder", &InitialOrder);
  fkeyread_double(fc, "ScalingFactor", &ScalingFactor);
  fkeyread_double(fc, "Young_A", &Dx); YoungA = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "Poisson_A", &PuasA);
  fkeyread_double(fc, "Young_M", &Dx); YoungM = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "Poisson_M", &PuasM);
  fkeyread_double(fc, "Heat expansion of austenite", &Dx); TECoefA = Dx * 1e-6;
  fkeyread_double(fc, "Heat expansion of martensite", &Dx); TECoefM = Dx * 1e-6;
  fkeyread_double(fc, "Heat conductivity of austenite", &HeatConductivityA); 
  fkeyread_double(fc, "Heat conductivity of martensite", &HeatConductivityM);
  fkeyread_double(fc, "Specific heat of austenite", &SpecificHeatA);
  fkeyread_double(fc, "Specific heat of martensite", &SpecificHeatM);
  fkeyread_double(fc, "Density of austenite", &Dx); DensityA = Dx *1e3; //into kg/cm^3
  fkeyread_double(fc, "Density of martensite", &Dx); DensityM = Dx *1e3;
  fkeyread_double(fc, "Decrement_A", &Decrement_A); //internal friction in austenite (at Fi==0)
  fkeyread_double(fc, "Decrement_M", &Decrement_M); //--------""------- in martensite (at Fi==1)
  fkeyread_double(fc, "Decrement_T", &Decrement_T); //--------""------- in two-phase state (at Fi==0.5)

  fkeyread_double(fc, "FYield", &Dx); FYield = Dx*1e6; // into Pa
  fkeyread_double(fc, "FYieldTwinning", &Dx); FYieldTwinning = Dx*1e6; // into Pa
  fkeyread_double(fc, "Kmp", &Kmp);
  fkeyread_double(fc, "Kap_MP", &Kap2002);
  fkeyread_double(fc, "Dolya", &Dolya);
  fkeyread_double(fc, "CorrectRoGr", &CorrectRoGr);
  fkeyread_double(fc, "Tau_c", &Dx); Tau_c = Dx*1e6; //into Pa - normalizer for TauYield


  if(!microplasticity_flag){//forbidding micro-plasticity
    FYield=1e100;
    Kap2002=0.0;
  }
  
  //read constants for damage
  // fkeyread_double(fc, "Al_d", &Al_d);//alfa (not used now)
  fkeyread_double(fc, "be_zv", &be_zv);

  //fkeyread_double(fc, "q2_F", &q2_F);
  fkeyread_double(fc, "q3_F", &q3_F);
  fkeyread_double(fc, "r0_F", &r0_F);
  fkeyread_double(fc, "U_F", &U_F);
  fkeyread_double(fc, "Tau_F", &Tau_F);  Tau_F*=1e6; // from Pa to MPa
  fkeyread_double(fc, "A_F", &A_F);
  fkeyread_double(fc, "K1_F", &K1_F);
  fkeyread_double(fc, "K2_F", &K2_F);
  fkeyread_double(fc, "Pi_cr", &Pi_cr);
  fkeyread_int(fc, "kzer_cr", &kzer_cr);
  

  if(!damage_flag){//forbidding damage accumulation
    be_zv = 1e100;
    //b_F*=1e6;
    //q2_F=0.0;
    q3_F=0.0;
  }

  fclose(fc);
          
  T0 = (Ms+Af)/2;
  Gm = -Q0/T0*(Ms-Mf);    // in J/m^3
  Ftr = -Q0*(Af-Ms) / (2.0*T0);    //in J/m^3

  if(kzer_cr==0) kzer_cr = (int)floor(My_pow(kzer,2.0/3.0));

  if (Dolya<1.0) Rm= -log(1.0-Dolya);
  else Rm=8.0;//this is very much (Dolya>0.999)

  printf("\n\nFazCnstINIT: Mf=%lg Ms=%lg t0=%lg q0=%lg ScalFac=%e"
	  ,Mf,Ms,T0,Q0,ScalingFactor);
  printf("\n Gm=%12.5lg Ftr=%12.5lg FYield=%12.5lg Kmp=%12.5lg"
	  ,Gm,Ftr,FYield,Kmp);
  printf("\n Dolya=%8.5lg Rm=%10.5lg"
	   ,Dolya,Rm);

  printf("\n YoungA=%9.5lg GPa YoungM=%9.5lg GPa"
	  ,YoungA*1e-9,YoungM*1e-9);
  printf("\n CorrectRoGr=%12.5lg Tau_c=%12.5lg"
	  ,CorrectRoGr,Tau_c);
  printf("\n\nConstants for damage:");
  printf("\n r0_F=%12.5lg U_F=%12.5lg be_zv=%12.5lg"
	    ,r0_F        ,U_F         ,be_zv);
  printf("\n  q3_F=%12.5lg Tau_F=%12.5lg A_F=%12.5lg K1_F=%12.5lg K2_F=%12.5lg"
	     ,q3_F        ,Tau_F,        A_F,        K1_F,        K2_F);
  printf("\n kzer=%3i kzer_cr=%3i" ,kzer ,kzer_cr);

  if(damage_flag && active_plasticity_flag && microplasticity_flag)
    printf("\nMicro-plasticity, active plasticity and damage: all are calculated");
  else{
    printf("\nNot calculated are:");
    if(!damage_flag) printf("\tdamage");
    if(!active_plasticity_flag) printf("\tactive_plasticity");
    if(!microplasticity_flag) printf("\tmicro-plasticity");
  }
  //pause();
  //pause("Init of Inpar successfully ended");
} //end of Init for MaterialConstants

  MaterialConstants::~MaterialConstants(void)
{
    delete[]KSlip;
    delete[]HSlip;
    delete[]USlip;
    delete[]R0Slip;
    delete[]TemperTable;//[KTauY]=   (K)
    //delete[]TauYieldTable;//[KTauY][KSlipTypes]=  (MPa)
    //this deletion was corrected
    for(int i=0;i<KSlipTypes;i++) delete[]TauYieldTable[i];
    delete[]TauYieldTable;
}


// Initializer reading cubic Matrixes from file "kubmatr.d" on disk
void VariantsOrientations::Init(MaterialConstants& MC)
{
   int l,k,i,j;
   FILE *fcube;
   FILE *fslip;
   Real Ox;

   kvar = MC.KVar;
   if( !(Orvvar_mp = new Matr33 [kvar] )){
     fprintf(stderr," Cannot allocate Orvvar_mp\n");
     pause(); exit(1);
   }

   if( !(Dev_var = new Matr33 [kvar] )){
     fprintf(stderr," Cannot allocate Dev_var\n");
     pause(); exit(1);
   }
   // Read cubic matrixes from file on disk
   if( (fcube=fopen("kubmatr.d","r"))==0){
	fprintf(stderr,"Cannot open kubmatr.d\n");
        pause();
   }
   for(k=0;k<kvar;k++){
     for(i=0;i<3;i++) for(j=0;j<3;j++){
       fscanf(fcube,"%lg",&Ox);
       Orvvar_mp[k].C(i,j)=Ox;
     }
     Dev_var[k] = MC.Dev_Dtrans.Povernut(1,Orvvar_mp[k]);
     Orvvar_mp[k] = MC.Dtrans.Povernut(1,Orvvar_mp[k]);
   }
   fclose(fcube);

   //Allocate memory for slip orientation matrixec
   OrSlip_mp = new Matr33*[MC.KSlipTypes];
   for(i=0;i<MC.KSlipTypes;i++){
     int k = MC.KSlip[i];
     OrSlip_mp[i] = new Matr33[k];
   }

   // Read matrixes of slip orientations from file on disk
   if( (fslip=fopen("slipmatr.d","r"))==0){
	fprintf(stderr,"Cannot open slipmatr.d\n");
        pause();
   }
   for(l=0;l<MC.KSlipTypes;l++)
   for(k=0;k<MC.KSlip[l];k++){
     for(i=0;i<3;i++) for(j=0;j<3;j++){
       fscanf(fslip,"%lg",&Ox);
       OrSlip_mp[l][k].C(i,j)=Ox;
     }
   }
   fclose(fslip);
}// End of initializer of VariantsOrientations

Real MaterialConstants::TauYieldEq(int iSlipType, Real Temper)// function to calculate
{                             //equilibrium yield limits for slip plane
                             // of iSlipSystem type and temperature Temper
  Real TY;
  int i;
  int Naiden=0;
  for(i=1; i<KTauY && Naiden==0 ;i++){
    if(Temper <= TauYieldTable[iSlipType][i]){
      TY=TauYieldTable[iSlipType][i-1]
         +(Temper-TemperTable[i-1])/(TemperTable[i]-TemperTable[i-1])
          *(TauYieldTable[iSlipType][i]-TauYieldTable[iSlipType][i-1]);
      Naiden=1;
    }
  }
  if(Naiden==0){
   int last=KTauY-1;
   TY=TauYieldTable[iSlipType][last-1]
      +(Temper-TemperTable[last-1])/(TemperTable[last]-TemperTable[last-1])
       *(TauYieldTable[iSlipType][last]-TauYieldTable[iSlipType][last-1]);
  }
  if(TY<0.0){
    fprintf(stderr,"\n\n TauYield<0 for Temperature=%lg .",Temper);
    fprintf(stderr," Check the table for TauYield in the file Mtrcnst.d");
    pause(); exit(2);
  }
  return(TY);
}// end of function MaterialConstants::TauYieldEq

// Initializer of class GrainsParameters
// by known kfi,kte,kps - uniform distribution
/*
void GrainsParameters::Init(int i1,int i2,int i3)
{
   Real fi,ps,hfi,hcte,hps,cf,sf,cte,ste,cp,sp;
   int i,j,k,izer;
   Matr33 R;

   kzer= (kfi=i1)*(kcte=i2)*(kps=i3);
*/

// Init for Grains Parameters
int GrainsParameters::Init()
{
   Real fi,ps,cf,sf,cte,ste,cp,sp;
   Real Alfa,Dist,Dist110;
// one must take Korobov's parameters for 4-th fold integrals
   Real a1=95.0,a2=122.0,a3=231.0;
   int p=307;
   int i,j,k,izer;
   int k_fi,k_cte,k_ps; //k_zer is a member of class GrainsParameters
   Matr33 R;
   FILE *f;

   if((f = fopen("cnst_key.d","r")) == 0){
     printf("Cannot open cnst_key.d");
     pause(); exit(1);
   }

  fkeyread_string(f, "GrainsOrientationsDistribution", GrainsOrientationsDistribution, 16);
   ffindkey(f, "k_fi,k_cte,k_ps");
   fscanf(f,"%i %i %i", &k_fi, &k_cte, &k_ps);
   fkeyread_int(f, "k_zer", &k_zer);     //for Korobov distribution this k_zer is taken
                                        //for Ravnomer. distribution it is calculated
                                       //for Monocrystal distributions it is =1
   fkeyread_double(f, "Alfa_mono", &Alfa); Alfa*=M_PI/180.0; //(from degrees to radians)
   fkeyread_double(f, "Dist_mono", &Dist);
   fkeyread_double(f, "Dist110_mono", &Dist110);
   fclose(f);
   
   printf("\nGrains orientations distribution= %s",GrainsOrientationsDistribution);
   //fprintf(freport,"\nGrains orientations distribution is %s\t",
   //           GrainsOrientationsDistribution);
  printf("\nGrainsINIT: k_fi=%5d k_cte=%5d k_ps=%5d k_zer=%5d"
          ,k_fi,k_cte,k_ps,k_zer);
  printf("\n Alfa=%12.5lg Dist=%12.5lg Dist110=%12.5lg"
	  ,Alfa,Dist,Dist110);
  //pause();

   if(GrainsOrientationsDistribution[0]=='K'){//Korobov - grains with orientations == Korobov's nodes
     Real intpart;
     if( !(Orvzer_mp = new Matr33 [k_zer] )){
       fprintf(stderr," Cannot allocate Orvzer_mp\n");
       pause(); exit(1);
     }
     izer=0;
     for(izer=0;izer<k_zer;izer++){
       fi= 2.0*M_PI*modf(a1*izer/p,&intpart);
       cte= 2.0*modf(a2*izer/p,&intpart)-1.0;
       ps= 2.0*M_PI*modf(a3*izer/p,&intpart);
       sf=sin(fi);cf=cos(fi);
       ste=sqrt(1.0-cte*cte);
       cp=cos(ps); sp=sin(ps);
 R.C(0,0)= cf*cp - sf*cte*sp; R.C(0,1)= -cf*sp - sf*cte*cp; R.C(0,2)= sf*ste;
 R.C(1,0)= sf*cp + cf*cte*sp; R.C(1,1)= -sf*sp + cf*cte*cp; R.C(1,2)= -cf*ste;
 R.C(2,0)= ste*sp;            R.C(2,1)= ste*cp;             R.C(2,2)= cte;
       Orvzer_mp[izer]=R;
     }
   }
   else if(GrainsOrientationsDistribution[0]=='M'){//Monocrystal
     k_zer = 1;
     fi=M_PI+Alfa;
     cte = cos(  2.0*atan2( Dist, Dist110/tan(M_PI/8.0))  );
     ps=0.0;
     if( !(Orvzer_mp = new Matr33 [k_zer] )){
       fprintf(stderr," Cannot allocate Orvzer_mp\n");
       pause(); exit(1);
     }
     sf=sin(fi);cf=cos(fi);
     ste=sqrt(1.0-cte*cte);
     cp=cos(ps); sp=sin(ps);
 //Matrix below is written in a transposed order because fi,te,ps
 //found above describe rotation transferring crystallagraphic basis
 //into laboratory one while the reverse rotation is needed
 R.C(0,0)= cf*cp - sf*cte*sp; R.C(1,0)= -cf*sp - sf*cte*cp; R.C(2,0)= sf*ste;
 R.C(0,1)= sf*cp + cf*cte*sp; R.C(1,1)= -sf*sp + cf*cte*cp; R.C(2,1)= -cf*ste;
 R.C(0,2)= ste*sp;            R.C(1,2)= ste*cp;             R.C(2,2)= cte;
     Orvzer_mp[0]=R;
   }
   else if(GrainsOrientationsDistribution[0]=='R'){//Ravnomernoe
     k_zer = k_fi*k_cte*k_ps;
     printf("k_fi=%d, k_cte=%d, k_ps=%d, k_zer=%d\n",k_fi,k_cte,k_ps,k_zer);
     // fprintf(freport,"k_fi=%d, k_cte=%d, k_ps=%d, k_zer=%d\n"
     //         ,k_fi,k_cte,k_ps,k_zer);
     if( !(Orvzer_mp = new Matr33 [k_zer] )){
       fprintf(stderr," Cannot allocate Orvzer_mp\n");
       pause(); exit(1);
     }
     Real hfi=2.0*M_PI/k_fi;
     Real hcte=2.0/k_cte;
     Real hps=2.0*M_PI/k_ps;
     izer=0;
     for(i=0;i<k_fi;i++){
       fi=hfi*(i+0.5); cf=cos(fi); sf=sin(fi);
       for(j=0;j<k_cte;j++){
         cte=-1.0 + hcte*(j+0.5);
	 ste=sqrt(1.0-cte*cte);
         for(k=0;k<k_ps;k++){
           ps=hps*(k+0.5); cp=cos(ps); sp=sin(ps);
 R.C(0,0)= cf*cp - sf*cte*sp; R.C(0,1)= -cf*sp - sf*cte*cp; R.C(0,2)= sf*ste;
 R.C(1,0)= sf*cp + cf*cte*sp; R.C(1,1)= -sf*sp + cf*cte*cp; R.C(1,2)= -cf*ste;
 R.C(2,0)= ste*sp;            R.C(2,1)= ste*cp;             R.C(2,2)= cte;
	   Orvzer_mp[izer++]=R;
	 }
       }
     }
   }
   else{
     fprintf(stderr
	   ,"\n Grains Orientations Distribution ""%s"" is not envisaged"
	   , GrainsOrientationsDistribution );
     return(1);
   }
   //fclose(f); already closed
   return(0);
}// End of constructor of class GrainsParameters


// Initializer of class InternalParameters
void InternalParameters::Init(MaterialConstants& MC
			      ,GrainsParameters& GP
			      //,int kl
			      //,int kz
			     )
{
  int k;
  kzer=GP.k_zer;
  kloc = MC.KVar * GP.k_zer;
  //kloc - total numb. of martensite variants in all grains
  kslip = 0;
  for(int i=0;i<MC.KSlipTypes;i++){
    kslip+=MC.KSlip[i];
  }
  kslip *= GP.k_zer;
  //kslip - total numb. of slip planes in all grains

  Bmean_pl_dmg=0.0; // Mean on all grains absolute value of B1 and B2
  Bmax_pl_dmg=0.0; // Maximum on all grains absolute value of B1 and B2
  Fmean_pl_dmg=0.0; // Mean on all grains value of irreversible defects
  Fmax_pl_dmg=0.0; // Maximum value of irreversible defects
  
  // allocation of parameters for plastic deformation and damage due to plastic deformation:
  if( !(TYDef_dp = new Real [kslip])){
    fprintf(stderr," Cannot allocate TYDef\n");
	pause(); exit(1);
  }
  if( !(Gb_dp = new Real [kslip])){
    fprintf(stderr," Cannot allocate Gb\n");
	pause(); exit(1);
  }
  if( !(Db_dp = new Real [kslip])){
    fprintf(stderr," Cannot allocate Db\n");
	pause(); exit(1);
  }
  if( !(B1_pl_dmg = new Real [kslip])){
    fprintf(stderr," Cannot allocate B1_pl_dmg\n");
	pause(); exit(1);
  }
  if( !(B2_pl_dmg = new Real [kslip])){
    fprintf(stderr," Cannot allocate B2_pl_dmg\n");
	pause(); exit(1);
  }
  if( !(F_pl_dmg = new Real [kslip])){
    fprintf(stderr," Cannot allocate F_pl_dmg\n");
	pause(); exit(1);
  }
  if( !(Pi1_pl_dmg = new Real [kslip])){
    fprintf(stderr," Cannot allocate Pi1_pl_dmg\n");
	pause(); exit(1);
  }
  if( !(Pi2_pl_dmg = new Real [kslip])){
    fprintf(stderr," Cannot allocate Pi2_pl_dmg\n");
	pause(); exit(1);
  }
  for(k=0;k<kslip;k++){
    TYDef_dp[k]=0.0;
    Gb_dp[k]=0.0;
    Db_dp[k]=0.0;
	B1_pl_dmg[k]=0.0;
	B2_pl_dmg[k]=0.0;
	F_pl_dmg[k]=0.0;
	Pi1_pl_dmg[k]=0.0;
	Pi2_pl_dmg[k]=0.0;
  }

//allocation of parameters for the transformation
  if( !(Floc_dp = new Real [kloc])){
    fprintf(stderr," Cannot allocate Floc\n");
    pause(); exit(1);
  }
  if( !(mpDef_dp = new Real [kloc])){
    fprintf(stderr," Cannot allocate mpDef\n");
    pause(); exit(1);
  }
  if( !(Fyloc_dp = new Real [kloc])){
    fprintf(stderr," Cannot allocate Fyloc\n");
    pause(); exit(1);
  }


  kzer_dmg=0; //initial number of destroyed grains
  PI_macro=0; //initial value for macro-fracture criterion - no fracture
  izer_bf_max=0; //initial number of grain with maximum fabs(b)*f (let it be 0)
//allocation of parameters for defects and damage caused by micro-plastic def-n
  if( !(B_dmg = new Real [kloc])){
    fprintf(stderr," Cannot allocate B_dmg\n");
    pause(); exit(1);
  }
  if( !(F_dmg = new Real [kloc])){
    fprintf(stderr," Cannot allocate F_dmg\n");
    pause(); exit(1);
  }
  if( !(pi_dmg = new Real [kloc])){
    fprintf(stderr," Cannot allocate pi_dmg\n");
    pause(); exit(1);
  }
  for(k=0;k<kloc;k++){
    Floc_dp[k]=0.0;
    mpDef_dp[k]=0.0;
    Fyloc_dp[k]=MC.FYield;
    B_dmg[k]=0.0;
    F_dmg[k]=0.0; //here initial value other F0_dmg than 0 can be assigned
    pi_dmg[k]=0.0;
  }

  if( !(RoGrain = new Matr33 [kzer])){
    fprintf(stderr," Cannot allocate RoGrain\n");
    pause(); exit(1);
  }
  if( !(hRoGrain = new Matr33 [kzer])){
    fprintf(stderr," Cannot allocate hRoGrain\n");
    pause(); exit(1);
  }
  if( !(PI_dmg = new Real [kzer])){
    fprintf(stderr," Cannot allocate PI_dmg\n");
    pause(); exit(1);
  }
  //Real* bf_dmg; //[kzer] max(on variants) product of b_dmg*f_dmg
  if( !(bf_dmg = new Real[kzer])){
    fprintf(stderr," Cannot allocate bf_dmg\n");
    pause(); exit(1);
  }
  for(k=0;k<kzer;k++){
    RoGrain[k]=0.0; //assign 0 to all elements of k-th matrix
    hRoGrain[k]=0.0;
    PI_dmg[k]=0.0;
    bf_dmg[k]=0.0;
  }

  if( !(EpsUnGr = new Matr33 [kzer])){
    fprintf(stderr," Cannot allocate EpsUnGr\n");
    pause(); exit(1);
  }
  for(k=0;k<kzer;k++){
    EpsUnGr[k]=0.0; //assign 0 to all elements of k-th matrix
  }

  Order = MC.InitialOrder;
  Phase=0.0;
  PhaseStrain=0.0;
  mpStrain=0.0;
  TEStrain=0.0;
  ElStrain=0.0;
  PlasticStrain=0.0;
  UnelasticStrain=0.0;
  TotalStrain=0.0;
  dTotalStrain=0.0; //time derivative of total strain
  return;
}
// End of initializer of class InternalParameters

 InternalParameters::~InternalParameters(void)
{
	delete[]Floc_dp;
	delete[]mpDef_dp;
	delete[]Fyloc_dp;
	delete[]TYDef_dp;
	delete[]Gb_dp;
	delete[]Db_dp;
	delete[]B1_pl_dmg;
	delete[]B2_pl_dmg;
	delete[]F_pl_dmg;
	delete[]Pi1_pl_dmg;
	delete[]Pi2_pl_dmg;
	delete[]RoGrain;
	delete[]hRoGrain;
	delete[]B_dmg;
	delete[]F_dmg;
	delete[]pi_dmg;
	delete[]PI_dmg;
	delete[]bf_dmg;
 }

void InternalParameters::operator=(InternalParameters& XX2)
{
  int k;
  if(kloc != XX2.kloc){
   fprintf(stderr
   ,"\nMembers kloc in two InternalParameters variables are not equal");
   pause(); exit(3);
  }

  for(k=0;k<kslip;k++){
    TYDef_dp[k]=XX2.TYDef_dp[k];
    Gb_dp[k]=XX2.Gb_dp[k];
    Db_dp[k]=XX2.Db_dp[k];
	B1_pl_dmg[k]=XX2.B1_pl_dmg[k];
	B2_pl_dmg[k]=XX2.B2_pl_dmg[k];
	F_pl_dmg[k]=XX2.F_pl_dmg[k];
	Pi1_pl_dmg[k]=XX2.Pi1_pl_dmg[k];
	Pi2_pl_dmg[k]=XX2.Pi1_pl_dmg[k];
  }

  for(k=0;k<kloc;k++){
    Floc_dp[k]=XX2.Floc_dp[k];
    mpDef_dp[k]=XX2.mpDef_dp[k];
    Fyloc_dp[k]=XX2.Fyloc_dp[k];
    B_dmg[k]=XX2.B_dmg[k];
    F_dmg[k]=XX2.F_dmg[k];
    pi_dmg[k]=XX2.pi_dmg[k];
  }

  for(k=0;k<kzer;k++){
    RoGrain[k]=XX2.RoGrain[k];
    hRoGrain[k]=XX2.hRoGrain[k];
    EpsUnGr[k]=XX2.EpsUnGr[k];
    PI_dmg[k]=XX2.PI_dmg[k];
    bf_dmg[k]=XX2.bf_dmg[k];
  }

  Bmean_pl_dmg=XX2.Bmean_pl_dmg; // Mean on all grains absolute value of B1 and B2
  Bmax_pl_dmg=XX2.Bmax_pl_dmg; // Maximum on all grains absolute value of B1 and B2
  Fmean_pl_dmg=XX2.Fmean_pl_dmg; // Mean on all grains value of irreversible defects
  Fmax_pl_dmg=XX2.Fmax_pl_dmg; // Maximum value of irreversible defects

  Phase=XX2.Phase;
  dPhase=XX2.dPhase;
  kzer_dmg=XX2.kzer_dmg;
  PI_macro=XX2.PI_macro;
  izer_bf_max=XX2.izer_bf_max;
  Order=XX2.Order;
  PhaseStrain=XX2.PhaseStrain;
  mpStrain=XX2.mpStrain;
  TEStrain=XX2.TEStrain;
  ElStrain=XX2.ElStrain;
  PlasticStrain=XX2.PlasticStrain;
  UnelasticStrain=XX2.UnelasticStrain;
  TotalStrain=XX2.TotalStrain;
  dTotalStrain=XX2.dTotalStrain;
}
