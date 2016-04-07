// taskporous.h
#ifndef _Task_LIBSAM

#define _Task_LIBSAM

// Definition of auxiliary structures
#include "defglob.h"
#include "matr.h"
#include "inpar.h"
#include "diaproc.h"
#include "reoeq1.h"


int Task_force_tem
  (
     int istep
	 ,Real htim
	 ,Real tem
	 ,Matr33& sig 
	 ,Matr33& hsig 
	 ,Matr33& sig_elastic
	 ,Matr33& hsig_elastic
	 ,Matr33& eps
	 ,Matr33& heps
	 ,MaterialConstants& MC
	 ,Porous_exVIA&  porous
	 ,VariantsOrientations&  Vrnts 
	 ,GrainsParameters& GrPar 
 );
  
int Task_force
  (
     int istep
    ,Real htim 
	//,Real hground_displ
	//,Real hdispl
	,Porous_exVIA& porous
    ,MaterialConstants& MC
    ,VariantsOrientations& Vrnts
    ,GrainsParameters& GrPar
	//InternalParameters& XXold,
 //   InternalParameters& XXnew
  );

int Task_force_2_slotted_elements
				(
					 Real htim
					,Real& hdispl
					,Porous_exVIA& porous
					,MaterialConstants& MC
					,VariantsOrientations& Vrnts
					,GrainsParameters& GrPar
				);

//int Task_force_2_slotted_elements
//	(
//		real htim		,real tem1		,real tem2  	,real htem1		,real htem2
//		,real force1	,real force2	,real hforce
//		,real& hdispl
//		,MaterialConstants& MC
//		,VariantsOrientations& Vrnts
//		,GrainsParameters& GrPar
//		,int Kuz
//		,Circular_beam* Cb1
//		,Circular_beam* Cb2
//	);

	int Task_VIA_operation
	( 
		 int istep	
		,int& VibrContinue	
		,Real htim	
		,Real hground_displ	
		,Porous_exVIA& porous	
		,MaterialConstants& MC
		,VariantsOrientations& Vrnts
		,GrainsParameters& GrPar
	);


int Task_Static_2_slotted_elements
	( 
		int istep	
		,Real htim
		,Porous_exVIA& porous	
		,class MaterialConstants& MC
		,class VariantsOrientations& Vrnts
		,class GrainsParameters& GrPar
	);

Real h_Force_on_hdispl
	//Return hforce1+hforce2 needed to produce specified hdispl 
	//						at specified porous.htem1 and porous.htem2
	(
		 Real htim
		,Porous_exVIA& porous
		,Real hdispl
		,class MaterialConstants& MC
		,class VariantsOrientations& Vrnts
		,class GrainsParameters& GrPar
	);

real Magnetic_force(real displacement, real scaling_factor);

int Task_eps
  ( int i ,int j//i,j denote the given strain component
   ,int istep       ,double htim
   ,double tem      ,double htem
   ,double fluence  ,double hfluence
   ,Matr33& sig   ,Matr33& hsig
   ,Matr33& heps
   ,MaterialConstants& MC
   ,VariantsOrientations& Vrnts
   ,GrainsParameters& GrPar
   ,InternalParameters& XXold
   ,InternalParameters& XXnew
  );

int Task_zhest
  ( int i,int j, int istep  ,double htim
   ,double tem ,double htem
   ,double fluence ,double hfluence ,double zhest
   ,Matr33& sig, Matr33& hsig
   ,Matr33& heps
   ,MaterialConstants& MC
   ,VariantsOrientations& Vrnts
   ,GrainsParameters& GrPar
   ,InternalParameters& XXold
   ,InternalParameters& XXnew
  );

int ReverseReology
  ( int i, int j ,int Continue,
    double htim,     double tem,  double htem,
    double fluence,  double hfluence,
    Matr33& sig,     Matr33& hsig,
	double force, double hforce,
    Matr33& heps,
    MaterialConstants& MC,
    VariantsOrientations& Vrnts,
    GrainsParameters& GrPar,
    InternalParameters& XXold,
    InternalParameters& XXnew
  );

int TaskVibr32(
    int     istep,    double  htim,
    double  tem,      double  htem,
    double  fluence,  double  hfluence,
    Matr33& sig,      Matr33& hsig,
    Matr33& eps,      Matr33& heps,
    class MaterialConstants& MC,
    class VariantsOrientations& Vrnts,
    class GrainsParameters& GrPar,
    class Majatnik& Maj1,
    class InternalParameters& XXold,
    class InternalParameters& XXnew
 );

int TaskVibr32_2elems( 
//                                                                                  //-->//
// excitation                                                                       //   //
// a*sin(wt)  /\    /\    /\          |-----------|       /\    /\    /\            //   //
//------->---/  \  /  \  /  \  /------|     m     |------/  \  /  \  /  \  /--------//   //
//               \/    \/    \/       |-----------|          \/    \/    \/         //   //
//           SMA helical coil 1           mass            SMA helical coil 2        //   //
//                                                                                  //-->//
//                                                                                  |    |
//                                                           preliminary_elongation-|<-->|
//                                                                                          
//Inrements of time, temperature are the same for the 1-st and the 2-nd elements
//                                                 
						int& VibrContinue   ,Real htim
						,Real tem           ,Real htem
						,Real fluence       ,Real hfluence
						,Matr33& sig1		,Matr33& hsig1
						,Matr33& sig2		,Matr33& hsig2
						,Real hground_displ ,Real hground_velocity
						,Real mass_displ    ,Real& hmass_displ
						,Real mass_velocity ,Real& hmass_velocity
						,class MaterialConstants& MC
						,class VariantsOrientations& Vrnts
						,class GrainsParameters& GrPar
						,class Isolation& Helix
						,class InternalParameters& XXold1
						,class InternalParameters& XXnew1
						,class InternalParameters& XXold2
						,class InternalParameters& XXnew2
);

int TaskVibr32_Isolation(
    int& VibrContinue   ,Real htim
   ,Real tem           ,Real htem
   ,Real fluence       ,Real hfluence
   ,Matr33& sig        ,Matr33& hsig
   ,Matr33& heps
   ,Real hground_angle ,Real hground_velocity
   ,Real element_displ ,Real& helement_displ
   ,Real element_velocity ,Real& helement_velocity
   ,class MaterialConstants& MC
   ,class VariantsOrientations& Vrnts
   ,class GrainsParameters& GrPar
   ,class Isolation& Helix
   ,class InternalParameters& XXold
   ,class InternalParameters& XXnew
);

#endif