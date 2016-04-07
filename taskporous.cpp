// Расчет декремента, обусловленного дефектами структуры,
// как функцмя фазы с максимумом при faza=1/2.
// Decr(0)=DecrA;Decr(1)=DecrM; Decr(1/2)=DecrT (DecrT>DecrM).
#include "DIAPROC.H"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "taskporous.h"
#define ITERVIBR 20
#define TOCHNOST 1E-4
#define ITER_MAX 20


int Task_force_tem

(
int istep,
Real htim
, Real tem
, Matr33& sig
, Matr33& hsig
, Matr33& sig_elastic
, Matr33& hsig_elastic
, Matr33& eps
, Matr33& heps
, class MaterialConstants& MC
, class Porous_exVIA& porous
, class VariantsOrientations& Vrnts
, class GrainsParameters& GrPar
)
{
	Real fluence = 0.0, hfluence = 0.0;
	const Real oshibka_max = 5e-8;//error of force in Newtons

	Real mass = porous.mass;

	Real Poly_A = 2 * MC.Decrement_M + 2 * MC.Decrement_A - 4 * MC.Decrement_T;
	Real Poly_B = 4 * MC.Decrement_T - MC.Decrement_M - 3 * MC.Decrement_A;
	Real Poly_C = MC.Decrement_A;
	Real elastic_compliance_total = 0.0;
	Real decr_total = 0.0;
	Real porosity = porous.porosity;
	int Kuz = porous.kc_beams;
	for (int i_beam = 0; i_beam < Kuz; i_beam++){
		Circular_beam *Cb = porous.c_beam + i_beam;//pointer to circular beam with number i_beam (for the 1-st RSS)
		Real faza = Cb->XXold.Phase;
		//Young moduli of two-phase material:
		Real E = 1.0 / ((1.0 - faza) / MC.YoungA + faza / MC.YoungM);
		//Decrements of two-phase material due to internal friction:
		Real decr_iuz = sqr(faza)*Poly_A + faza*Poly_B + Poly_C;
		real elastic_compliance_iuz = 1/*pi4pi * cube(Cb->radius) / Cb->J_r */ / E;
		decr_total += decr_iuz * elastic_compliance_iuz;
		elastic_compliance_total += elastic_compliance_iuz;
	}
	real K1 = 1.0 / elastic_compliance_total;
	real K = K1; //combined rigidity of both elements
	decr_total *= K;//decrement

	porous.hdispl_1 = 0.0;
	porous.hdispl_2 = 0.0;
	porous.hdispl_3 = 0.0;
	porous.hdispl_4 = 0.0;
	porous.hdispl_5 = 0.0;
	porous.hdispl_6 = 0.0;

	for (int iuz = 0; iuz < Kuz; iuz++){
		Real ks = porous.c_beam[iuz].K_sig;
		Real ka = porous.c_beam[iuz].K_sig2;
		Real df = porous.hforce;
		//sig(2, 2) = ks * porous.force_old;
		//hsig(2, 2) = ks * df;
	//	sig(1, 1) = -ka * porous.force_old;
	//	hsig(1, 1) = -ka * df;
		double Sigma_cr = 760e6; // critical value for stress of the beam, in Pa
		double Sigma_cr2 = 1000e6;// critical value for stress on tension, Pa
	//	real sig33 = sig(2, 2);
		Real faza = porous.c_beam[iuz].XXold.Phase;
		Real E = 1.0 / ((1.0 - faza) / MC.YoungA + faza / MC.YoungM);
if (porous.c_beam[iuz].Check_destr_old == 1.0){ //if the beam was destroyed
	fprintf(stderr, "\n  destroyed_beam");
	porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].XXold.TotalStrain(2, 2);
	porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = porous.c_beam[iuz].XXold.TotalStrain(1, 1);
	porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
	sig(1, 1) = 0.0;
											  }


else {
	
		if (porous.c_beam[iuz].phoenix == 1.0) {//if we already have had the contact
	
				if (((porous.tem_new-porous.tem_old) == 0.0)||(porous.c_beam[iuz].loading_type==1.0)){
					hsig(1, 1) = -ka*df;
					real force_im = porous.force_old - porous.c_beam[iuz].force_fix;
					sig(1, 1) = -ka*force_im;
					//sig(2, 2) = 0.0;
					hsig(2, 2) = 0.0;
					fprintf(stderr, "\n  connect by case 1");

					real raps = porous.c_beam[iuz].XXold.TotalStrain(1, 1);
					real meow = porous.c_beam[iuz].XXnew.TotalStrain(1, 1);
					porous.c_beam[iuz].XXold.TotalStrain(2, 2) = porous.c_beam[iuz].e33_fix;
					Reology(htim, porous.tem_old, porous.htem, fluence, hfluence, sig, hsig, MC, Vrnts, GrPar, porous.c_beam[iuz].XXold, porous.c_beam[iuz].XXnew, porous.c_beam[iuz].phoenix );
			//		porous.c_beam[iuz].XXnew.TotalStrain(1, 1) -= porous.c_beam[iuz].e22_el_fix;
					porous.c_beam[iuz].heps_elastic = porous.c_beam[iuz].XXnew.TotalStrain(1, 1) - porous.c_beam[iuz].XXold.TotalStrain(1, 1);

					porous.c_beam[iuz].hdisplac = 0.0;
					printf("\n sig(1,1)=%lg", sig(1, 1));
					printf("\n e22=%lg", porous.c_beam[iuz].XXnew.TotalStrain(1,1));
					printf("\n e22=%lg", porous.c_beam[iuz].XXold.TotalStrain(1, 1));
					//pause();
					porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
					porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].e33_fix;
					porous.c_beam[iuz].L_macro_new = porous.c_beam[iuz].L_macro_old *(1 + porous.c_beam[iuz].heps_elastic);
					porous.c_beam[iuz].loading_type = 1.0;

																									}
			else{
					sig(1, 1) = -ka*porous.force_old;
					sig(2, 2) = 0.0;
					hsig(2,2) = 0.0;
					hsig(1, 1) = -ka*df;
					fprintf(stderr, "\n  connect by case 2");
					pause();
					Reology(htim, porous.tem_old, porous.htem, fluence, hfluence, sig, hsig, MC, Vrnts, GrPar, porous.c_beam[iuz].XXold, porous.c_beam[iuz].XXnew, porous.c_beam[iuz].phoenix);
					porous.c_beam[iuz].heps_elastic = porous.c_beam[iuz].XXnew.TotalStrain(1, 1) - porous.c_beam[iuz].XXold.TotalStrain(1, 1);
					porous.c_beam[iuz].hdisplac = 0.0;
					porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
					porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].e33_fix;
					//porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = porous.c_beam[iuz].heps_elastic + porous.c_beam[iuz].eps_elastic_old;
					porous.c_beam[iuz].L_macro_new = porous.c_beam[iuz].L_macro_old *(1 + porous.c_beam[iuz].heps_elastic);
				}
		

		    if (abs(sig_elastic(1, 1)) >= Sigma_cr2){
				porous.c_beam[iuz].Check_destr_new = 1.0;
				fprintf(stderr, "\n  destroyed_beam");
				porous.c_beam[iuz].heps_elastic = 0.0; 
				porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = porous.c_beam[iuz].XXold.TotalStrain(1, 1);
				 porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
													 }
			//porous.c_beam[iuz].XXold.TotalStrain(1, 1) = porous.c_beam[iuz].XXold.TotalStrain(1, 1);
		//pause();
	                                       }// end of case when phoenix[iuz] =1.0



	  else { //if we haven't had the contact yet
			sig(1, 1) = 0.0;
			hsig(1, 1) = 0.0;
			hsig(2, 2) = ks * df;
			sig(2, 2) = ks * porous.force_old;
			
			//porous.c_beam[iuz].XXold.TotalStrain(1, 1) = 0.0; //added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//porous.c_beam[iuz].XXold.dTotalStrain(1, 1) = 0.0;
			//porous.c_beam[iuz].XXold.UnelasticStrain(1, 1) = 0.0;
			//porous.c_beam[iuz].XXold.ElStrain(1, 1) = 0.0;
			//porous.c_beam[iuz].XXold.mpStrain(1, 1) = 0.0;
			//porous.c_beam[iuz].XXold.PhaseStrain(1, 1) = 0.0;

		    Reology(htim, porous.tem_old, porous.htem, fluence, hfluence, sig, hsig, MC, Vrnts, GrPar, porous.c_beam[iuz].XXold, porous.c_beam[iuz].XXnew, porous.c_beam[iuz].phoenix);
			printf("\n sig(2,2)<500e6 and Check destr <1: sig(2,2)=%lg", sig(2, 2));
			porous.c_beam[iuz].heps_beam_new = porous.c_beam[iuz].XXnew.TotalStrain(2, 2) - porous.c_beam[iuz].XXold.TotalStrain(2, 2);
			porous.c_beam[iuz].Check_destr_new = 0.0;
			porous.c_beam[iuz].hdisplac = porous.c_beam[iuz].heps_beam_new*(-8 * cube(porous.c_beam[iuz].length) + 4 * porous.c_beam[iuz].length*sqr(porous.c_beam[iuz].little_length) - cube(porous.c_beam[iuz].little_length)) /
				(24*porous.c_beam[iuz].thickness*(2 * porous.c_beam[iuz].length - porous.c_beam[iuz].little_length));
		    porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old + porous.c_beam[iuz].hdisplac;
			porous.c_beam[iuz].L_macro_new = porous.c_beam[iuz].thickness + porous.c_beam[iuz].column - abs (porous.c_beam[iuz].displac_new);
			porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = 0.0;

				if (abs(porous.c_beam[iuz].displac_new) >= (porous.c_beam[iuz].column / 2)){
					porous.c_beam[iuz].phoenix = 1.0;
					printf("\n sig(2,2)=%lg", sig(2, 2));
					printf("\n eps(2,2)=%lg", porous.c_beam[iuz].XXnew.TotalStrain(2, 2));
					printf("\n eps(1,1)=%lg", porous.c_beam[iuz].XXnew.TotalStrain(1, 1));
					fprintf(stderr, "\n  connected");
					pause();
					porous.c_beam[iuz].s33_fix = sig(2, 2);
					//porous.c_beam[iuz].e22_el_fix = porous.c_beam[iuz].XXnew.ElStrain(1, 1);
					porous.c_beam[iuz].e33_fix = porous.c_beam[iuz].XXnew.TotalStrain(2, 2);
					porous.force_new = porous.force_old + porous.hforce;
					porous.c_beam[iuz].force_fix = porous.force_new;
					porous.c_beam[iuz].hdisplac = 0.0;
					porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
					porous.c_beam[iuz].L_macro_fix = porous.c_beam[iuz].L_macro_new;
					porous.c_beam[iuz].XXnew.dTotalStrain(1, 1) = 0.0;
					porous.c_beam[iuz].XXnew.UnelasticStrain(1, 1) = 0.0;
					porous.c_beam[iuz].XXnew.ElStrain(1, 1) = 0.0;
					porous.c_beam[iuz].XXnew.mpStrain(1, 1) = 0.0;
					porous.c_beam[iuz].XXnew.PhaseStrain(1, 1) = 0.0;
				//    porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = 0.0;///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
																							}

				if (abs(sig(2, 2)) >= Sigma_cr){
					porous.c_beam[iuz].Check_destr_new = 1.0;
					fprintf(stderr, "\n  destroyed_beam");
					porous.c_beam[iuz].heps_elastic = 0.0;
					porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].XXold.TotalStrain(2, 2);
					porous.c_beam[iuz].displac_new = porous.c_beam[iuz].displac_old;
					pause();
												}
			}			
	}
fprintf(stderr, "\n  first_circle_if_else");
//if ((porous.c_beam[iuz].phoenix == 1.0) && (porous.c_beam[iuz].XXnew.TotalStrain(1, 1)>porous.c_beam[iuz].e22_fix)){
if ((porous.c_beam[iuz].phoenix == 1.0) && (porous.c_beam[iuz].L_macro_new > porous.c_beam[iuz].L_macro_fix)){
    porous.c_beam[iuz].phoenix = 0.0;
	fprintf(stderr, "\n  disconnect");
	pause();
																											 }
porous.c_beam[iuz].sig_to_control = sig(1, 1)+hsig(1,1);
porous.c_beam[iuz].sigma33 = sig(2,2)+hsig(2,2);
real sig33 = sig(2, 2);

real fdshdf = porous.c_beam[iuz].XXnew.TotalStrain(1, 1);// delete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 }//end of 'iuz'
		
		porous.area_macro = (porous.c_beam[0].width*porous.c_beam[0].little_length)*(1 - porosity)*(-1 + porosity + porous.c_beam[0].noV) / (-1 + porous.c_beam[0].noV);
		porous.force_new = porous.force_old + porous.hforce;
		porous.force_new_eff = porous.force_new*0.01003*0.01086*0.151 / porous.area_macro;
		porous.L_macro_old_overall = porous.c_beam[0].L_macro_old/* + porous.c_beam[1].L_macro_old + 4*porous.c_beam[2].L_macro_old + 4*porous.c_beam[3].L_macro_old + 2*porous.c_beam[4].L_macro_old*/;
		porous.L_macro_new_overall = porous.c_beam[0].L_macro_new/* + porous.c_beam[1].L_macro_new + 4*porous.c_beam[2].L_macro_new + 4*porous.c_beam[3].L_macro_new + 2*porous.c_beam[4].L_macro_new*/;
		porous.de_macro = (porous.L_macro_new_overall - porous.L_macro_old_overall) / porous.L_macro_old_overall;
		porous.e_macro += porous.de_macro;
		porous.stress_macro = porous.force_new / porous.area_macro;
real rehtdkgl = porous.c_beam[0].XXnew.TotalStrain(1, 1);// delete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		return 0;
	}