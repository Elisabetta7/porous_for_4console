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
, class Porous_forVert& vert
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

	int Vuz = vert.kv_beams;
	for (int i_beam = 0; i_beam < Vuz; i_beam++){
		Vertical_beam *Vb = vert.v_beam + i_beam;//pointer to circular beam with number i_beam (for the 1-st RSS)
		Real faza = Vb->XXold.Phase;
		//Young moduli of two-phase material:
		Real E = 1.0 / ((1.0 - faza) / MC.YoungA + faza / MC.YoungM);
		//Decrements of two-phase material due to internal friction:
		Real decr_iuv = sqr(faza)*Poly_A + faza*Poly_B + Poly_C;
		real elastic_compliance_iuv = 1/*pi4pi * cube(Cb->radius) / Cb->J_r */ / E;
		decr_total += decr_iuv * elastic_compliance_iuv;
		elastic_compliance_total += elastic_compliance_iuv;
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

		double Sigma_cr = 600e6; // critical value for stress of the beam, in Pa
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
					sig(2, 2) = porous.c_beam[iuz].s33_fix;
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
					printf("\n hsig(1,1)=%lg", hsig(1, 1));
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
				//	pause();
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
			porous.c_beam[iuz].L_macro_new = porous.c_beam[iuz].thickness + porous.c_beam[iuz].column + porous.c_beam[iuz].displac_new;
			porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = 0.0;
			//porous.c_beam[iuz].XXold.PlasticStrain(1, 1) = 0.0;

				if (abs(porous.c_beam[iuz].displac_new) >= (porous.c_beam[iuz].column / 2)){
					porous.c_beam[iuz].phoenix = 1.0;
					printf("\n sig(2,2)=%lg", sig(2, 2));
					printf("\n eps(2,2)=%lg", porous.c_beam[iuz].XXnew.TotalStrain(2, 2));
					printf("\n eps(1,1)=%lg", porous.c_beam[iuz].XXnew.TotalStrain(1, 1));
					fprintf(stderr, "\n  connected");
	//				pause();
					porous.c_beam[iuz].s33_fix = (sig(2, 2)+hsig(2,2)/2);
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
					porous.c_beam[iuz].XXnew.PlasticStrain(1, 1) = 0.0;

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
	porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].XXold.TotalStrain(2, 2);
	hsig(1, 1) = 0.0;
	porous.c_beam[iuz].XXnew.TotalStrain(1, 1) = porous.c_beam[iuz].XXold.TotalStrain(1, 1);
	fprintf(stderr, "\n  disconnect");
//	pause();
																											 }
porous.c_beam[iuz].sig_to_control = sig(1, 1)+hsig(1,1);
porous.c_beam[iuz].sigma33 = sig(2,2)+hsig(2,2);
real sig33 = sig(2, 2);
real fdshdf = porous.c_beam[iuz].XXnew.TotalStrain(1, 1);// delete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 }//end of 'iuz'
		








	 for (int iuv = 0; iuv<Vuz; iuv++){
		 //int Check_destr = 0; // constant for checking of destroy
		 Real ksv = vert.v_beam[iuv].K_sig_old;
		 Real df = porous.hforce;
		 sig(2, 2) = ksv * porous.force_old+df;
		 hsig(2, 2) = ksv * df;
		 /*printf ("\n iuz=%u, sig(2,2)=%lg, hsig(2,2)=%lg", iuz , sig(2,2), hsig(2,2));
		 pause();*/
		 double Sigma_cr = 700e6; // critical value for stress of the beam, in Pa
		 real sig33 = sig(2, 2);

		 if ((sig(2, 2) >= Sigma_cr) && (vert.v_beam[iuv].Check_destr_old == 0.0)){
			 fprintf(stderr, "\n  destroyed1");
			 printf("\n sig(2,2) >= Sigma_cr: sig(2,2)=%lg", sig(2, 2));
			 pause();
			 vert.v_beam[iuv].Check_destr_new = 1.0;
			 vert.v_beam[iuv].heps_beam_new = 0.0;
			 vert.v_beam[iuv].XXnew.TotalStrain(2, 2) = vert.v_beam[iuv].XXold.TotalStrain(2, 2) + vert.v_beam[iuv].heps_beam_new;// +porous.c_beam[iuz].heps_beam_new;
			 vert.v_beam[iuv].phoenix = 1.0;
			 //porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].XXold.TotalStrain(2, 2) + porous.c_beam[iuz].heps_beam_old;
		 }
		 else {
			 if (vert.v_beam[iuv].Check_destr_old == 2.0){//╙┴╨└╥▄ 2!!!!!!!!!!!!!!!!!!!! ╥└╠ 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 vert.v_beam[iuv].heps_beam_new = 0.0;
				 vert.v_beam[iuv].phoenix = 0.0;
				 vert.v_beam[iuv].XXnew.TotalStrain(2, 2) = vert.v_beam[iuv].XXold.TotalStrain(2, 2) + vert.v_beam[iuv].heps_beam_new;
				 fprintf(stderr, "\n  destroyed2");
				 printf("\n Check destr>1: sig(2,2)=%lg", sig(2, 2));
				 //pause();
				 vert.v_beam[iuv].Check_destr_new = 1.0;
			 }
			 else {
				 Reology(htim, vert.tem_old, vert.htem, fluence, hfluence, sig, hsig
					 , MC, Vrnts, GrPar, vert.v_beam[iuv].XXold, vert.v_beam[iuv].XXnew, vert.v_beam[iuv].phoenix);
				 printf("\n sig(2,2)<500e6 and Check destr <1: sig(2,2)=%lg", sig(2, 2));
				 //porous.c_beam[iuz].heps_beam_new = porous.c_beam[iuz].XXnew.TotalStrain(2, 2) - porous.c_beam[iuz].XXold.TotalStrain(2, 2);
				 //porous.c_beam[iuz].XXnew.TotalStrain(2, 2) = porous.c_beam[iuz].heps_beam_old + porous.c_beam[iuz].XXold.TotalStrain(2, 2);
				 vert.v_beam[iuv].heps_beam_new = vert.v_beam[iuv].XXnew.TotalStrain(2, 2) - vert.v_beam[iuv].XXold.TotalStrain(2, 2);
				// vert.v_beam[iuv].XXnew.TotalStrain(2, 2) = vert.v_beam[iuv].heps_beam_new + vert.v_beam[iuv].XXold.TotalStrain(2, 2);
				 vert.v_beam[iuv].Check_destr_new = 0.0;
				 vert.v_beam[iuv].phoenix = 0.0;
			 }

		 }
		 if (vert.v_beam[iuv].XXnew.TotalStrain(2, 2) >= 1.0){
			 fprintf(stderr, "\n  destroyed3");
			 printf("\n eps's too large ");
			 //pause();
			 vert.v_beam[iuv].Check_destr_new = 1.0;
			 vert.v_beam[iuv].heps_beam_new = 0.0;
			 vert.v_beam[iuv].XXnew.TotalStrain(2, 2) = vert.v_beam[iuv].XXold.TotalStrain(2, 2) + vert.v_beam[iuv].heps_beam_new;
		 }

		 Real faza = vert.v_beam[iuv].XXold.Phase;
		 Real E = 1.0 / ((1.0 - faza) / MC.YoungA + faza / MC.YoungM);




		 vert.v_beam[iuv].radius_new = vert.v_beam[iuv].thickness / (2 * (vert.v_beam[iuv].heps_beam_new + vert.v_beam[iuv].thickness / (2 * vert.v_beam[iuv].radius_old)));
		 vert.v_beam[iuv].fi_angle_new = ((vert.v_beam[iuv].radius_old*vert.v_beam[iuv].radius_old_old*vert.v_beam[iuv].fi_angle_old*vert.v_beam[iuv].fi_angle_old_old)*(vert.v_beam[iuv].heps_beam_new) + sqr(vert.v_beam[iuv].radius_old*vert.v_beam[iuv].fi_angle_old)) / (vert.v_beam[iuv].fi_angle_old_old*vert.v_beam[iuv].radius_old_old*vert.v_beam[iuv].radius_new);

		 //	porous.c_beam[iuz].fix = (1 + porous.c_beam[iuz].heps_beam)*porous.c_beam[iuz].radius_old*log(porous.c_beam[iuz].radius_old / (porous.c_beam[iuz].radius_old - porous.c_beam[iuz].thickness));
		 //	porous.c_beam[iuz].radius_new = Newton(porous.c_beam[iuz].fix, porous.c_beam[iuz].thickness);
		 //	porous.c_beam[iuz].fi_angle_new = (1 + porous.c_beam[iuz].heps_beam)*porous.c_beam[iuz].radius_old*porous.c_beam[iuz].fi_angle_old / porous.c_beam[iuz].radius_new;
		 if (vert.v_beam[iuv].fi_angle_new >= 3.14){
			 fprintf(stderr, "\n fi_new > pi");
			 printf("\n force=%lg, temperature=%lg", vert.force_new, vert.tem_new);
			 vert.v_beam[iuv].Check_destr_new = 1.0;
			 vert.v_beam[iuv].fi_angle_new = vert.v_beam[iuv].fi_angle_old;
			 vert.v_beam[iuv].radius_new = vert.v_beam[iuv].radius_old;
			// pause(); 
		 }

	     vert.hdispl_1 = 2 * (2 * vert.v_beam[0].radius_old*sin(vert.v_beam[0].fi_angle_old / 2) - vert.v_beam[0].radius_new*sin(vert.v_beam[0].fi_angle_new / 2) - vert.v_beam[0].radius_old_old*sin(vert.v_beam[0].fi_angle_old_old / 2));
		 vert.hdispl_2 = 2 * (2 * vert.v_beam[1].radius_old*sin(vert.v_beam[1].fi_angle_old / 2) - vert.v_beam[1].radius_new*sin(vert.v_beam[1].fi_angle_new / 2) - vert.v_beam[1].radius_old_old*sin(vert.v_beam[1].fi_angle_old_old / 2));
		 //vert.hdispl_1 = 2 * (vert.v_beam[0].radius_old*sin(vert.v_beam[0].fi_angle_old / 2) - vert.v_beam[0].radius_new*sin(vert.v_beam[0].fi_angle_new / 2));
		 //vert.hdispl_2 = 2 * (vert.v_beam[1].radius_old*sin(vert.v_beam[1].fi_angle_old / 2) - vert.v_beam[1].radius_new*sin(vert.v_beam[1].fi_angle_new / 2));
		 //porous.hdispl_2 = 2*(2*porous.c_beam[1].radius_old*sin(porous.c_beam[1].fi_angle_old/2)-porous.c_beam[1].radius_new*sin(porous.c_beam[1].fi_angle_new/2)-porous.c_beam[1].radius_old_old*sin(porous.c_beam[1].fi_angle_old_old/2));
		 vert.hdispl_3 = 2 * (2 * vert.v_beam[2].radius_old*sin(vert.v_beam[2].fi_angle_old / 2) - vert.v_beam[2].radius_new*sin(vert.v_beam[2].fi_angle_new / 2) - vert.v_beam[2].radius_old_old*sin(vert.v_beam[2].fi_angle_old_old / 2));
		 vert.hdispl_4 = 2 * (2 * vert.v_beam[3].radius_old*sin(vert.v_beam[3].fi_angle_old / 2) - vert.v_beam[3].radius_new*sin(vert.v_beam[3].fi_angle_new / 2) - vert.v_beam[3].radius_old_old*sin(vert.v_beam[3].fi_angle_old_old / 2));
		 vert.hdispl_5 = 2 * (2 * vert.v_beam[4].radius_old*sin(vert.v_beam[4].fi_angle_old / 2) - vert.v_beam[4].radius_new*sin(vert.v_beam[4].fi_angle_new / 2) - vert.v_beam[4].radius_old_old*sin(vert.v_beam[4].fi_angle_old_old / 2));
		 vert.hdispl_6 = 2 * (2 * vert.v_beam[5].radius_old*sin(vert.v_beam[5].fi_angle_old / 2) - vert.v_beam[5].radius_new*sin(vert.v_beam[5].fi_angle_new / 2) - vert.v_beam[5].radius_old_old*sin(vert.v_beam[5].fi_angle_old_old / 2));

		 if ((vert.v_beam[0].K_sig_new <= 0.0) && (vert.v_beam[0].thickness / vert.v_beam[0].radius >= 0.1)){
			 fprintf(stderr, "\n bad bulk");
			 printf("\n beam=%u", 1);
			// pause(); 
		 }
		 if ((vert.v_beam[0].K_sig_new <= 0.0) && (vert.v_beam[0].thickness / vert.v_beam[0].radius <= 0.1)){
			 vert.hdispl_1 = abs(vert.hdispl_1);
		 }


		 if ((vert.v_beam[1].K_sig_new <= 0.0) && (vert.v_beam[1].thickness / vert.v_beam[1].radius >= 0.1)){
			 fprintf(stderr, "\n bad bulk");
			 printf("\n beam=%u", 2);
			 //pause(); 
		 }
		 if ((vert.v_beam[1].K_sig_new <= 0.0) && (vert.v_beam[1].thickness / vert.v_beam[1].radius <= 0.1)){
			 vert.hdispl_2 = abs(vert.hdispl_2);
		 }


		 if ((vert.v_beam[2].K_sig_new <= 0.0) && (vert.v_beam[2].thickness / vert.v_beam[2].radius >= 0.1)){
			 fprintf(stderr, "\n bad bulk");
			 printf("\n beam=%u", 3);
			 //pause(); 
		 }
		 if ((vert.v_beam[2].K_sig_new <= 0.0) && (vert.v_beam[2].thickness / vert.v_beam[2].radius <= 0.1)){
			 vert.hdispl_3 = abs(vert.hdispl_3);
		 }



		 if ((vert.v_beam[3].K_sig_new <= 0.0) && (vert.v_beam[3].thickness / vert.v_beam[3].radius >= 0.1)){
			 fprintf(stderr, "\n bad bulk");
			 printf("\n beam=%u", 4);
			 //pause(); 
		 }
		 if ((vert.v_beam[3].K_sig_new <= 0.0) && (vert.v_beam[3].thickness / vert.v_beam[3].radius <= 0.1)){
			 vert.hdispl_4 = abs(vert.hdispl_4);
		 }


		 if ((vert.v_beam[4].K_sig_new <= 0.0) && (vert.v_beam[4].thickness / vert.v_beam[4].radius >= 0.1)){
			 fprintf(stderr, "\n bad bulk");
			 printf("\n beam=%u", 5);
			 //pause(); 
		 }
		 if ((vert.v_beam[4].K_sig_new <= 0.0) && (vert.v_beam[4].thickness / vert.v_beam[4].radius <= 0.1)){
			 vert.hdispl_5 = abs(vert.hdispl_5);
		 }
		 

		 vert.hdispl = 8*vert.hdispl_1 + /* 5*vert.hdispl_2 +*/ vert.hdispl_3  /*+ porous.hdispl_4 + porous.hdispl_5 */+
			8*vert.v_beam[0].phoenix*(vert.v_beam[0].radius_old*sin(vert.v_beam[0].fi_angle_old / 2)) / (vert.k_parallel - 1 * 1)+
			/*5*vert.v_beam[1].phoenix*(vert.v_beam[1].radius_old*sin(vert.v_beam[1].fi_angle_old / 2)) / (vert.k_parallel - 1 * 1)+*/
			vert.v_beam[2].phoenix*(vert.v_beam[2].radius_old*sin(vert.v_beam[2].fi_angle_old / 2)) / (vert.k_parallel - 1 * 1)/*+
																																			3.5*porous.c_beam[3].phoenix*(2 * vert.v_beam[3].radius_old*sin(vert.v_beam[3].fi_angle_old / 2)) / (porous.k_parallel - 1 * 1)+
																																			2*porous.c_beam[4].phoenix*(2 * vert.v_beam[4].radius_old*sin(vert.v_beam[4].fi_angle_old / 2)) / (porous.k_parallel - 1 * 1)*/;

		 vert.radius_new_1 = vert.v_beam[0].radius_new;
		 vert.radius_old_1 = vert.v_beam[0].radius_old;
		 vert.fi_angle_new_1 = vert.v_beam[0].fi_angle_new;
		 vert.fi_angle_new_2 = vert.v_beam[1].fi_angle_new;
		 vert.fi_angle_old_1 = vert.v_beam[0].fi_angle_old;
		 vert.fi_angle_old_2 = vert.v_beam[1].fi_angle_old;
		 vert.fi_angle_new_3 = vert.v_beam[2].fi_angle_new;
		 vert.fi_angle_new_4 = vert.v_beam[3].fi_angle_new;
		 vert.fi_angle_old_3 = vert.v_beam[2].fi_angle_old;
		 vert.fi_angle_old_4 = vert.v_beam[3].fi_angle_old;


		 vert.v_beam[iuv].K_sig_new = (-1) / (vert.v_beam[iuv].width*vert.v_beam[iuv].thickness) +
			 2 * (vert.v_beam[iuv].fi_angle_new / 2 - sin(vert.v_beam[iuv].fi_angle_new / 2))*(vert.v_beam[iuv].radius_new - vert.v_beam[iuv].thickness / ((log(vert.v_beam[iuv].radius_new / (vert.v_beam[iuv].radius_new - vert.v_beam[iuv].thickness)))))
			 / (vert.v_beam[iuv].fi_angle_new*vert.v_beam[iuv].width*vert.v_beam[iuv].thickness*(vert.v_beam[iuv].radius_new - vert.v_beam[iuv].thickness / 2 - vert.v_beam[iuv].thickness / ((log(vert.v_beam[iuv].radius_new / (vert.v_beam[iuv].radius_new - vert.v_beam[iuv].thickness))))));

		 //real first_part_K_sig = porous.c_beam[iuz].fi_angle_new / (porous.c_beam[iuz].fi_angle_new - 2 * sin(porous.c_beam[iuz].fi_angle_new/2));
		 //real second_part_K_sig = (porous.c_beam[iuz].radius_new - porous.c_beam[iuz].thickness / ((log(porous.c_beam[iuz].radius_new / (porous.c_beam[iuz].radius_new - porous.c_beam[iuz].thickness)))))
		 /// (porous.c_beam[iuz].radius_new - porous.c_beam[iuz].thickness / 2 - porous.c_beam[iuz].thickness / ((log(porous.c_beam[iuz].radius_new / (porous.c_beam[iuz].radius_new - porous.c_beam[iuz].thickness)))));

		 //if (porous.c_beam[iuz].K_sig_new<=0.0){
		 //	//porous.c_beam[iuz].K_sig_new <= 0.0){
		 //	fprintf(stderr, "\n K_sig < 0");
		 //	printf("\n beam=%u", iuz);
		 //	pause(); 
		 //}


		 vert.dlina1 = vert.v_beam[0].Length;
		 vert.dlina2 = vert.v_beam[1].Length;
		 vert.dlina3 = vert.v_beam[2].Length;
		 vert.dlina4 = vert.v_beam[3].Length;
		 vert.dlina5 = vert.v_beam[4].Length;
		 vert.dlina6 = vert.v_beam[5].Length;

		 vert.check1 = vert.v_beam[0].Check_destr_new;
		 vert.check2 = vert.v_beam[1].Check_destr_new;
		 vert.check3 = vert.v_beam[2].Check_destr_new;
		 vert.check4 = vert.v_beam[3].Check_destr_new;
		 vert.check5 = vert.v_beam[4].Check_destr_new;

		 vert.dlina_nov1 = 2 * vert.v_beam[0].radius_new*sin(vert.v_beam[0].fi_angle_new / 2);
		 vert.dlina_nov2 = 2* vert.v_beam[1].radius_new*sin(vert.v_beam[1].fi_angle_new / 2);
		 vert.dlina_nov3 = 2 * vert.v_beam[2].radius_new*sin(vert.v_beam[2].fi_angle_new / 2);
		 vert.dlina_nov4 = 2 * vert.v_beam[3].radius_new*sin(vert.v_beam[3].fi_angle_new / 2);
		 vert.dlina_nov5 = 2 * vert.v_beam[4].radius_new*sin(vert.v_beam[4].fi_angle_new / 2);

		 vert.dlina_star1 = 2 * vert.v_beam[0].radius_old*sin(vert.v_beam[0].fi_angle_old / 2);
		 vert.dlina_star2 = 2 * vert.v_beam[1].radius_old*sin(vert.v_beam[1].fi_angle_old / 2);
		 vert.dlina_star3 = 2 * vert.v_beam[2].radius_old*sin(vert.v_beam[2].fi_angle_old / 2);
		 vert.dlina_star4 = 2 * vert.v_beam[3].radius_old*sin(vert.v_beam[3].fi_angle_old / 2);
		 vert.dlina_star5 = 2 * vert.v_beam[4].radius_old*sin(vert.v_beam[4].fi_angle_old / 2);



		 real pore_diam1 = vert.v_beam[0].width*porosity / (1 - 1 * porosity);
		 real pore_diam2 = vert.v_beam[0].thickness*porosity / (1 - 1 * porosity);
		 vert.area_macro = (vert.v_beam[0].width*vert.v_beam[0].thickness)/**(1-porosity)*/;
		 //porous.area_macro = (porous.c_beam[0].width + pore_diam1)*(porous.c_beam[0].thickness + pore_diam2);
		 //}	
	 }//end for(int iuz=0; iuz<Kuz ;  iuz++)


	 vert.displ_new = vert.displ_old + vert.hdispl;






		porous.area_macro = (porous.c_beam[0].width*porous.c_beam[0].little_length)*sqr(porous.c_beam[0].noV - 1 + porosity)/sqr(-1 + porous.c_beam[0].noV);
		porous.force_new = porous.force_old + porous.hforce;
		porous.force_new_eff = porous.force_new*0.01003*0.01086*0.151 / porous.area_macro;
		porous.L_macro_old_overall = porous.c_beam[0].L_macro_old + 8*vert.dlina_star1/* + porous.c_beam[1].L_macro_old */+ 2*porous.c_beam[2].L_macro_old + 0.3*porous.c_beam[3].L_macro_old/* + 5*vert.dlina_star2*/ + vert.dlina_star3;
		porous.L_macro_new_overall = porous.c_beam[0].L_macro_new + 8*vert.dlina_nov1 /*+ porous.c_beam[1].L_macro_new */+ 2*porous.c_beam[2].L_macro_new + 0.3*porous.c_beam[3].L_macro_new/* + 5*vert.dlina_nov2 */+ vert.dlina_nov3;
		porous.de_macro = (porous.L_macro_new_overall - porous.L_macro_old_overall) / porous.L_macro_old_overall;
		porous.e_macro += porous.de_macro;
		porous.stress_macro = porous.force_new / porous.area_macro;
real rehtdkgl = porous.c_beam[0].XXnew.TotalStrain(1, 1);// delete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		return 0;
	}