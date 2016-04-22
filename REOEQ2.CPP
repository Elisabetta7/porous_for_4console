// Reology with RoGrain equilibration
// 12/09/06 new (stress) condition of micro failure
// Reverse transformation causes relaxation of microplastic flow stress 
// Accurate calculation of Twinning direction.
// Strain rate dependence of the flow stress is accounted for.

// Formulae for damage
  // db/dt = dbe/dt - (1/be*)(dbe/dt)*abs(b)*H(b*dbe/dt)
  // df/dt = (q2+q3*f)*abs(b)^alfa*abs(dbe/dt) - r(T)*f
  // where r(T) = r0*exp(-U/kT)
  //Condition of micro failure:
  // tau_concentrated being the critical value
  // tau_concentrated=fabs(b)*(1.0+A_F*f)*tau_loc,
  // tau_loc=SpSig+tauzer:Dn/Gam(Dn)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <alloc.h>
#include <time.h>
#include "defglob.h"
#include "matr.h"
#include "inpar.h"
#include "reoeq1.h"
#include "diaproc.h"
#define VNIZ 1
#define VVERH 0
#define IATTEMPTMAX 10
#define ROOT 2
#define MAX_FORBIDDEN_TWDIR 4

//#define ROOT 1.414213562
//#define ROOT4 1.183207115

Matr33 Hook(Real E, Real Nu, Matr33 &S)
{
 Matr33 Eps;
 int i;
 Real da=Nu/E*( S.C(0,0)+S.C(1,1)+S.C(2,2) );
 Eps = S; Eps *= (1+Nu)/E;
 for(i=0;i<3;i++) Eps.C(i,i) -= da;
 return Eps;
}//hook

Matr33 HookS(Real E, Real Nu, Matr33 &Eps)
{
 Matr33 Sig(Eps);
 int i;
 Real da=Nu/(1.0-2.0*Nu)*( Eps.C(0,0)+Eps.C(1,1)+Eps.C(2,2) );
 for(i=0;i<3;i++) Sig.C(i,i) += da;
 Sig *= E/(1+Nu);
 return Sig;
}//HookS

Real E_Young(Real phase, MaterialConstants& MC)
{
   return 1.0/( (1.0-phase)/MC.YoungA + phase/MC.YoungM );
}

void CalculationOfHysterisisParameters
              (
	       MaterialConstants& MC
              ,Real fluence
              ,Real& T0
              ,Real& Ftr
              ,Real& Gm
              ,Real& Q0
              )
{
  Real kT0_fluence=0.0, kFtr_fluence=0.0;
  T0 = MC.T0 - kT0_fluence * fluence;
  Ftr= MC.Ftr - kFtr_fluence*fluence;
  Gm = MC.Gm;
  Q0 = MC.Q0;
}

Real ThermodynamicForce(int n
                       ,VariantsOrientations& VO
                       ,MaterialConstants& MC
                       ,Real T
                       ,Matr33& Sig
                       )
{
  if(n<0 || n>VO.kvar){
    fprintf(stderr,"\n\tThermodynamicForce: N of variant = %d out of range"
            ,n);
    //pause();
  }
  Real SpSig=Sig.C(0,0)+Sig.C(1,1)+Sig.C(2,2);
  Real Teta = T - MC.T0;
  Real S0_AM = MC.Q0/MC.T0;
  Real Cs_AM = MC.SpecificHeatM - MC.SpecificHeatA;
  Real Alfa_AM = MC.TECoefM - MC.TECoefA;
  Matr33 Eps_elA = Hook(MC.YoungA, MC.PuasA, Sig);
  Matr33 Eps_elM = Hook(MC.YoungM, MC.PuasM, Sig);
  Eps_elM -= Eps_elA;
  Real Ft =  S0_AM*Teta + Cs_AM*(T*log(T/MC.T0)-Teta)
           + Alfa_AM*Teta*SpSig + Sig.Svertka(VO[n])
           + 0.5 * Eps_elM.Svertka(Sig);
  return(Ft);
}

Real ReducedGibbsForceForTwinning(int n
            ,VariantsOrientations& VO
            ,Matr33& sig)
{
  return(sig.Svertka(VO[n]));
}


void Reology(
              Real htim
              ,Real Tem
	      ,Real hTem
              ,Real Fluence
              ,Real DeltaFluence
              ,Matr33 &Sigma
              ,Matr33 &DeltaSigma
			  /*,Real force
			  ,Real hforce*/
              ,MaterialConstants& MC
              ,VariantsOrientations& Vrnts
	      ,GrainsParameters& GrPar
              ,InternalParameters &XXold
	      ,InternalParameters &XXnew
	,Real phoenix
          )
{

 if(htim<=0.0) Exit_error("\nReo: invalid htim=%lg", htim);
     //printf("\nReo_100:htim=%lg",htim);
     //pause();

 int i, isis, izer;
 int k_sis=Vrnts.kvar;
 int k_zer=GrPar.k_zer;
 const Real RBolts = 8.31e-3; // KJ/(K*mol) Universal gas constant
 Real A_rad_creep=0.0;
 Matr33 hrad_strain;
 //unsigned k_loc=k_sis*k_zer;
 int isisglob=0;
 Real Fs = 0.0;
 Real hFazTw = 0.0;
 Real FazaNew=0.0, Fazer,FazerNew,hfazp,hfazo,fzl,hfzl;
// Real Tem = Tem+hTem;
 //Temperature Arrhenius factor for defects F_dmg recovery:
 Real r_dmg = MC.r0_F*exp(-MC.U_F/(RBolts*Tem));
 Real F_dmg_glob=0.0;
 Real T0,Ftr,Gm,Q0;
 //Real mn = MC.Ms, ak = MC.Af;
 //Real mn_mk = mn - MC.Mf, ak_an = ak - MC.As;

     //printf("\nReo_200");
     //pause();

 CalculationOfHysterisisParameters(MC,Fluence,T0,Ftr,Gm,Q0);

    // printf("\nReo_300:T0=%lg Ftr=%lg Gm=%lg Q0=%lg",T0,Ftr,Gm,Q0);
    // pause();


 Real Hm=Gm*(1.0-MC.Kmp)/MC.Kmp;
 Real One_minus_k = 1.0-MC.Kmp;
 //Real Gm=MC.Gm, Ftr=MC.Ftr;
 Real FYield;
 Real cx ,corr;
 Matr33 sig(Sigma);
 sig+=DeltaSigma;
 Real SpSig=sig.C(0,0)+sig.C(1,1)+sig.C(2,2);


 Matr33 rzer ,tauzer ,tauzerold ,temp;
 Matr33 hbetp ,hbzph ,hbzphp ,hbzmp ,hbzap ,hbzmpp ,hbzapp;
 Matr33 dEPh ,dEmp ,dETE ,dEAP;
 Real *hFloc;     hFloc     = new Real [k_sis];
 Real *hFlocTw;   hFlocTw   = new Real [k_sis];
 Real *hFi_MP;    hFi_MP    = new Real [k_sis];
 Real *hFi_MP_Tw; hFi_MP_Tw = new Real [k_sis];
 Real *KMP_Tw;    KMP_Tw    = new Real [k_sis];
 Real *TwDir;     TwDir     = new Real [k_sis];
 Real *F_mix;     F_mix     = new Real [k_sis];
 Real *FTw_n;     FTw_n     = new Real [k_sis];
 int *tw_dir_is_allowed; tw_dir_is_allowed = new int [k_sis];

 Real FTnew, FTold, delFT, delFspr,delFsobr, Kmp;
  // FazaNew = 0; in the declaration

 //Prepare for RoGrain finding self consistent procedure:
// randomize();
 Real AdmissibleErrRo=5.0e4; //Pa = 0.05 MPa
 Real ErrRo, ErrRop;
 Real LamRo=0.9;
 Real LamRop=0.9;
 Real MaxRoNorma;
 Matr33 *RoGrProb, *RoGrProbp, *RoGrzv;//suffix p means previous
 RoGrProb = new Matr33[k_zer];
 RoGrProbp = new Matr33[k_zer];
 RoGrzv = new Matr33[k_zer];
 for(izer=0;izer<k_zer;izer++){
   //RoGrProb[izer] = XXold.RoGrain[izer];
   RoGrProb[izer] = Komb(1.0, XXold.RoGrain[izer], LamRo, XXold.hRoGrain[izer]);
 }
 int RoGrainIsFound = 0;
 //int FixPlasticVariants = 0;
    // FixPlasticVariants may be put =1 to avoid uncertainty in
    // plastic variants determination and to reach convergence
 int IAttempt=0;
 int IterRo=1;
 int NumOfGrWithTwin=0;

   //  printf("\nReo_400");
   //  pause();

 for(;RoGrainIsFound==0 && IAttempt<IATTEMPTMAX;IterRo++){
 //while(RoGrainIsFound==0 && IAttempt<IATTEMPTMAX){
   //IterRo++;
   FazaNew=0.0;
   dEPh = 0.0;
   dEAP = 0.0;
   dEmp = 0.0;
   NumOfGrWithTwin=0;

   XXnew.Bmax_pl_dmg=0.0;
   XXnew.Bmean_pl_dmg=0.0;
   XXnew.Fmax_pl_dmg=0.0;
   XXnew.Fmean_pl_dmg=0.0;
// ================================= open cycle on izer ===========
   for(izer=0;izer<k_zer;izer++){
    Fazer=0.0;
    int smesh=izer*k_sis;//grain numb. times numb. of syst. in a grain
    int smesAP=izer*MC.KSlipPlanesOfAllTypes;
    //shift for Slip Plane number
    for(isis=0;isis<k_sis;isis++){
     Fazer += XXold.Floc_dp[smesh+isis];
    }
                                //volume fract of M in a zer
    Fazer /= k_sis;

    if(Fazer<0.0 || Fazer>1.0+1.0e-9){
      fprintf(stderr,"\nBad value of Fazer=%lg",Fazer);
      exit(2);
    }

    rzer= GrPar.ZerOr(izer);
              temp = Sigma+XXold.RoGrain[izer];
    tauzerold=temp.Povernut(VNIZ,rzer);
              temp= sig+RoGrProb[izer]; //+XXold.hRoGrain[izer];
    tauzer=temp.Povernut(VNIZ,rzer);
    hfazp = 0.0; hfazo = 0.0;
  // -------------open main cycle on isis (loc.s) ----------------

     //printf("\nReo_500");
     //pause();

    for(isis = 0;isis<k_sis;isis++){
      int isisglob=smesh+isis;
      fzl = XXold.Floc_dp[isisglob];
      Real fzlmp = *(XXold.mpDef_dp + isisglob);
      temp = Vrnts[isis];
      FTnew = ThermodynamicForce(isis,Vrnts,MC,Tem,tauzer);
      FTold = ThermodynamicForce(isis,Vrnts,MC,Tem,tauzerold);
      delFT = FTnew - FTold;
      Fs = Gm * (fzl - fzlmp);  //resistance force (self-variant term)
      delFspr = FTnew - Ftr - Fs;
      delFsobr = FTnew + Ftr - Fs;
      FYield = *(XXold.Fyloc_dp + smesh + isis);

      hFi_MP[isis] = Kmp =
        Fs+delFT>FYield && delFspr>0  ||  Fs+delFT < -FYield && delFsobr<0 ?
        // if this condition of microplastic flow holds for variant # isis
        // then we assign to the incr. of Fi_MP the value of the multiplier:
        MC.Kmp :
        // otherwise this multiplier is 0.0
        0.0;

      hfzl = (delFspr > 0 ? delFspr / (Gm * (1-Kmp)) : 0)
           + (delFsobr < 0 ? delFsobr / (Gm * (1-Kmp)) : 0) ;

      if( fzl + hfzl <0 ) hfzl = -fzl;
      hFloc[isis] = hfzl;
                                     // memorized uncorrected hfazl
      if(hfzl<0) hfazo += hfzl;
      if(hfzl>0) hfazp += hfzl;
    }//for(isis) - close main cycle on isis -------------------

     //printf("\nReo_600");
     //pause();

    hfazo /= k_sis;
    hfazp /= k_sis;
    // --------------------begin correction----------------------------
    if( hfazp>0 && Fazer+hfazo+hfazp>1.0 ){
      cx = Fazer + hfazo;
      if(cx>0.9999)
        corr = 0.0;
      else
        corr = (1.0 - cx) / hfazp ;
      for (isis = 0; isis < k_sis; isis++) {
        if(hFloc[isis] > 0) hFloc[isis] *= corr;
      }//for isis
    }//if -----------------end correction--------------------------------

    //Calcul.incr.of Fi_MP.  Before Kmp was assigned to hFi_MP
    // and now it is multiplied onto hFloc:
    for (isis = 0; isis < k_sis; isis++) {
      hFi_MP[isis] *= hFloc[isis];
    }

    //==================================================================
    //================= begin calculation of twinning ==================
    //==================================================================
    int twinning_may_be=1;
    //Intermediate calculation of new phase amount in the grain under consideration
    //Futher on this operation is done again
    FazerNew=0.0;
    for(isis = 0 ; isis < k_sis ; isis++){
      //int isisglob=smesh+isis;
      FazerNew += XXold.Floc_dp[smesh+isis] + hFloc[isis];
      hFlocTw[isis] = 0.0;
      hFi_MP_Tw[isis] = 0.0;
      //here we shall store the increments due to consequent steps of twinning
      //calculated in the cycle while(twinning_may_be)
    }
    if(FazerNew<1e-5){ twinning_may_be=0; }//no twinning when little martensite
    int iter_twin=0;
    double TwdNorm=0.0;
    int n0,n1;
    //Open cycle on steps corresponding to several displacements
    //in the plane SUM{Fi_n}=const :
    while(twinning_may_be && iter_twin++<100){  //iter_twin++<20
      //if(iter_twin>8){ - someties it happens that iter_twin==9
        //printf("\n iter_twin=%i",iter_twin);
        //pause();
      //}
      Real TwdAv = 0.0;
	  Real Ftw_norm=0.0;
      for (isis = 0; isis < k_sis; isis++) {
        isisglob=smesh+isis;
        //9.04.2006: possible increments due to phase trans-ns are taken into account
        //12.04.2006 added are additions due to twinning (found on previous iterations
        //           of cycle while(twinning_may_be...
        Real fzl = XXold.Floc_dp[isisglob] + hFloc[isis] + hFlocTw[isis];
        Real fzlmp = XXold.mpDef_dp[isisglob] + hFi_MP[isis] + hFi_MP_Tw[isis];
        //  this version was before 9.04.2006
        //Real fzl = XXold.Floc_dp[isisglob];
        //Real fzlmp = XXold.mpDef_dp[isisglob];
	Real Fs =
		F_mix[isis] = Gm * (fzl - fzlmp);  //resistance force

	FTw_n[isis] = TwDir[isis]
	= ReducedGibbsForceForTwinning(isis,Vrnts,tauzer) - Fs;

        TwdAv += TwDir[isis];
		Ftw_norm += sqr(FTw_n[isis]);
      }//for isis
      TwdAv /= k_sis;
	  Ftw_norm = sqrt(Ftw_norm);
      
      //if the force Ftw_norm is too small we need not make any projections
	  //as twinning cannot occur
	  Ftw_norm -= MC.FYieldTwinning;
      if(Ftw_norm<=0){ twinning_may_be=0; break; }//break from twinning
	  
	  
	  //Projecting TwDir on the octahedral plane:
      for (isis = 0; isis < k_sis; isis++) { TwDir[isis] -= TwdAv; }//for isis

 	  //Begin 2nd projection to exclude twinning,which makes negative variants
      //Initially we assume that all twinning directions are allowed:
      for (isis = 0; isis < k_sis; isis++){ tw_dir_is_allowed[isis]=1; }

      int there_are_forbidden_tw_dir = 1;//flag to stop projecting
      int projection_count;
      for(projection_count=0; there_are_forbidden_tw_dir; projection_count++){
      //while(there_are_forbidden_tw_dir){
	TwdAv=0.0; //TwdAv is used in new sense for the second projection
        n0=n1=0;
        for (isis = 0; isis < k_sis; isis++){
          if(tw_dir_is_allowed[isis]==0){n0++; continue;}//pass to next isis
          //(this forbidden direction was excluded in the previous iteration) 
	  if(XXold.Floc_dp[smesh+isis]+hFloc[isis]<1e-3 && TwDir[isis]<0.0){
          //projection of tw.direction can not be <0 if variant_vol.fract.== 0
            tw_dir_is_allowed[isis]=0;
	    TwdAv += TwDir[isis];//summing up increments in forbidden directions
	    TwDir[isis]=0.0;
	    n0++;
	  }
          else{
            n1++;
          }
        }
	if(n0>MAX_FORBIDDEN_TWDIR){//too many directions are forbidden
          twinning_may_be=0; break;// breaking from cycle for projecting
        }

        TwdAv/=n1;
        TwdNorm = 0.0;
        //calcul.for projecting on the boundary of triangle on octahedral plane
        for (isis = 0; isis < k_sis; isis++){
	  if(tw_dir_is_allowed[isis]){
	    TwDir[isis]+=TwdAv;
	    TwdNorm += my_square(TwDir[isis]);//my_square(x)=x*x - my function
	  }
        }
        TwdNorm=sqrt(TwdNorm);
        //Check if any forbidden directions are still present:
        there_are_forbidden_tw_dir=0;
        for (isis = 0 ; isis<k_sis ; isis++){
	  if(XXold.Floc_dp[smesh+isis]+hFloc[isis]<1e-3 && TwDir[isis]<0.0){
            there_are_forbidden_tw_dir=1;
            break;//at least one unchecked forbidden twinning direction is left
          }
        }//for isis
      }//while there_are_forbidden_tw_dir
      if (twinning_may_be==0) break;//break from cycle for twinning

      for (isis = 0 ; isis<k_sis ; isis++){
	if(XXold.Floc_dp[smesh+isis]+hFloc[isis]<1e-3 && TwDir[isis]<0.0){
	  fprintf(stderr,"\nReo: ERROR: Twinning dir makes negative variants");
	  pause();
	  exit(2);
		}
      }//for isis

      Real Twd_sum = 0.0;
      for (isis = 0; isis < k_sis; isis++) Twd_sum += TwDir[isis];
      if(fabs(Twd_sum)>1e-6){
	fprintf(stderr,"\nReo: ERROR: TwDir does not correspond to twinning");
	fprintf(stderr,"\nprojection_count=%i Twd_sum=%15.12lg  Ftw_norm=%15.12lg",projection_count,Twd_sum,Ftw_norm*1e-6);
	fprintf(stderr,"\nsig_norma=%lg ",sig.Norma()*1e-6);
	pause("\nProgram aborted.");
	exit(2);
      }

      if(TwdNorm>1e-3){
	    for (isis = 0; isis < k_sis; isis++){ TwDir[isis] /= TwdNorm; }
      }
      else { twinning_may_be=0; break; }//breaking the twinning cycle
      //................................................................
      //::::::::::reorientation direction TwDir is found::::::::::::::::
      //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      //find twinning force FTw_lam in the direction TwDir:
      Real FTw_lam=0.0;
      for (isis = 0; isis < k_sis; isis++) {
	FTw_lam += FTw_n[isis]*TwDir[isis];
      }//for isis

      //now the condition of twinning:
      FTw_lam -= MC.FYieldTwinning;
      if(FTw_lam<=0){ twinning_may_be=0; break; }//break from twinning
      else { NumOfGrWithTwin++; }
      // from this place FTw_lam has the role of its increment

      //Question:
	//Is it good to assume that hFloc[isis] is proportional to TwDir
	//rather than that hFloc[isis]-hFlocP[isis] is proportional to TwDir?

	  //find the amount (norm) of the phase increment vector:
      Real hFi_norm=0.0;
      for (isis = 0; isis < k_sis; isis++) {
	int isisglob=smesh+isis;
	Real FYield = *(XXold.Fyloc_dp + isisglob);
	if( fabs(F_mix[isis] - FYield) >= 0.0){
	  hFi_norm += my_square(TwDir[isis]) * One_minus_k;//One_minus_k = 1.0-MC.Kmp;
	  KMP_Tw[isis] = MC.Kmp;
	}
	else{
	  hFi_norm += my_square(TwDir[isis]);
	  KMP_Tw[isis]=0.0;
	  //KMP_Tw = stores the multiplier for incr.of hFi_MP_Tw
	}
      }//for isis
      hFazTw = FTw_lam / (Gm*hFi_norm);
      // hFazTw is the thermodyn. equlibr. length of the Fi-vector increment

      //Find the minimum of all admissible values of hFazTw
      Real hFazTw_min = hFazTw;
      Real hFazTw_admissible = hFazTw;
      for (isis = 0; isis < k_sis; isis++) {
	if(TwDir[isis]<0.0){
	  hFazTw_admissible = - XXold.Floc_dp[smesh+isis] / TwDir[isis];
        }
        if(hFazTw_admissible<hFazTw_min) hFazTw_min = hFazTw_admissible;
      }//for isis

      if(hFazTw == hFazTw_min) twinning_may_be=0;
      //when this is true all possible twinning will be calculated in the
      //present iteration of the cycle "while(twinning_may_be)"

      hFazTw = hFazTw_min;
      for (isis = 0; isis < k_sis; isis++) {
        Real incr_hFloc = hFazTw*TwDir[isis];
        hFlocTw[isis] += incr_hFloc;
        hFi_MP_Tw[isis] += KMP_Tw[isis] * incr_hFloc;
        // Here we are summing up the increments to hFlocTw and hFi_MP_Tw
        // due to the current step in the cycle while(twinning_may_be...)
      }
    }//end of while-cycle (twinning_may_be && iter_twin++<?)
   //======================================================================
   //================= end of the calculation of twinning =================
   //======================================================================
     //printf("\nReo_700");
     //pause();


    XXnew.PI_dmg[izer]=XXold.PI_dmg[izer]; 
	//if XXold.PI_dmg[izer]was==0 and at least one variant is with damage in cycle on isis or iplane
    //XXnew.PI_dmg[izer] will be assigned  1.0
	
//-------------------------------------------------------------------------
//------------- Calculation of athermal plastic deformation ---------------
//-------------------------------------------------------------------------
    if(MC.active_plasticity_flag){ //if flag is not ==0 we calculate active plasticity
      hbzap = 0.0; //Matrix to sum up Active PLASTIC strains for all systems
      int iPlaneGlob=smesAP; //this variable will follow global numbers of planes
      Real Cond_dmg_1 = 0.0;
	  Real Cond_dmg_2 = 0.0;
	  // open cycle on slip planes
      for(int iSlipType=0;iSlipType<MC.KSlipTypes;iSlipType++){
        int KPlanes=MC.KSlip[iSlipType];
        for(int iPlane=0;iPlane<KPlanes;iPlane++,iPlaneGlob++){
			Matr33 rPlane = Vrnts.SlipOr(iSlipType,iPlane);
			Matr33 TauPlane = tauzer.Povernut(VNIZ,rPlane);//tauzer is the new stress
			Real TTau = TauPlane.SIntensivityOn3Plane();
			Real hGb=1.7e308, RelaxTY=1.7e308;
			//Matr33 TauPlaneOld=tauzerold.Povernut(VNIZ,rPlane);//tauzer is the new stress
			//Real TTauOld = TauPlaneOld.SIntensivityOn3Plane();

			//Calculate the relaxation term of the yield stress (RelaxTY)
			if(MC.R0Slip[iSlipType]>0.0){
				RelaxTY = MC.R0Slip[iSlipType]
					* exp( -MC.USlip[iSlipType]/(RBolts*Tem) )
					* My_pow(XXold.TYDef_dp[iPlaneGlob]/MC.Tau_c, 4.0)
					*htim;
				//NSlip must be included into class MaterialConstants
				//RBolts = 8.31e-3 KJ/(K*mol) is defined in "defglob.h"
			}
			else RelaxTY=0.0;

			if(MC.PSlip[iSlipType]>0.0){//there is strain rate
							//dependence of the flow stress
				XXnew.TYDef_dp[iPlaneGlob] = XXold.TYDef_dp[iPlaneGlob]
				+ (MC.HSlip[iSlipType]*XXold.Db_dp[iPlaneGlob] - RelaxTY)*htim;

				Real Db_dp_temp = ( TTau - MC.TauYieldEq(iSlipType, Tem)
						- XXnew.TYDef_dp[iPlaneGlob] ) /MC.PSlip[iSlipType];

				if(Db_dp_temp<=0.0) //the flow condition does not hold
							XXnew.Db_dp[iPlaneGlob]=0.0;

				else{ //the flow condition holds
				XXnew.Db_dp[iPlaneGlob] =
				My_pow(Db_dp_temp, 1.0/MC.AlfaSlip[iSlipType]);
				//printf("\n Db_dp_temp =%lg",Db_dp_temp);
				//printf(" MC.AlfaSlip=%lg",MC.AlfaSlip[iSlipType]);
				//printf("\n XXnew.Db_dp[iPlaneGlob] =%lg",XXnew.Db_dp[iPlaneGlob]);
				//pause();
				}
				hGb = XXnew.Db_dp[iPlaneGlob]*htim;
				XXnew.Gb_dp[iPlaneGlob] = XXold.Gb_dp[iPlaneGlob] + hGb;
			}
			else{ //there is no strain rate dependence of the flow stress
				hGb = (TTau - MC.TauYieldEq(iSlipType, Tem)
						- XXold.TYDef_dp[iPlaneGlob]
						//- XXold.hLatent[izer] latent hardening is not accounted
						+ RelaxTY ) /MC.HSlip[iSlipType];
				if(hGb<0.0) hGb=0.0; //this plays the role of flow condition
				XXnew.TYDef_dp[iPlaneGlob] = XXold.TYDef_dp[iPlaneGlob]
							+ MC.HSlip[iSlipType]*hGb - RelaxTY;
			}
			// end of hGb calculation
			Matr33 hBALoc;
			hBALoc=0.0;
			Real hB1_pl_dmg=0.0;
			Real hB2_pl_dmg=0.0;
			Real hF_pl_dmg=0.0;
				Real B1_0 = XXold.B1_pl_dmg[iPlaneGlob];
				Real B2_0 = XXold.B2_pl_dmg[iPlaneGlob];
				Real F_pl_0 = XXold.F_pl_dmg[iPlaneGlob];

				//if(fabs(B1_0) > 0.0 || fabs(B2_0) > 0.0){
				//	printf("\n Br710: F_pl_0=%lg, B1_0=%lg, B2_0=%lg"
				//		,F_pl_0, B1_0, B2_0);
				//	pause();
				//}
				//	

			if(TTau>0.0){
				hBALoc(2,0) = hBALoc(0,2)= 0.5*hGb*TauPlane(2,0)/TTau;
				hBALoc(2,1) = hBALoc(1,2)= 0.5*hGb*TauPlane(2,1)/TTau;
				hbzap += hBALoc.Povernut(VVERH,rPlane);
			//----------Calculation of defect density at athermal plastic deformation----------		
				Real b_na_hBA = B1_0*hBALoc(0,2) + B2_0*hBALoc(1,2); // scalar product B on hBALoc
				Real abs_B=sqrt(B1_0*B1_0 + B2_0*B2_0);
				if(abs_B>0){
					hB1_pl_dmg = hBALoc(2,0)- (1/MC.be_zv)*XHeavisideX(b_na_hBA)*B1_0/abs_B;
					hB2_pl_dmg = hBALoc(2,1)- (1/MC.be_zv)*XHeavisideX(b_na_hBA)*B2_0/abs_B;
					hF_pl_dmg = (1+MC.q3_F*F_pl_0)*abs_B*hGb - r_dmg*F_pl_0*htim;
				}
				else{
					hB1_pl_dmg = hBALoc(2,0);
					hB2_pl_dmg = hBALoc(2,1);
					hF_pl_dmg = 0.0;
				}
				XXnew.B1_pl_dmg[iPlaneGlob] = B1_0 + hB1_pl_dmg;
				XXnew.B2_pl_dmg[iPlaneGlob] = B2_0 + hB2_pl_dmg;
				Real B_norm = sqrt(sqr(XXnew.B1_pl_dmg[iPlaneGlob])+sqr(XXnew.B2_pl_dmg[iPlaneGlob]));
				XXnew.F_pl_dmg[iPlaneGlob] = F_pl_0 + hF_pl_dmg;
				if(XXnew.Bmax_pl_dmg<B_norm) XXnew.Bmax_pl_dmg=B_norm;
				if(XXnew.Fmax_pl_dmg<XXnew.F_pl_dmg[iPlaneGlob]) XXnew.Fmax_pl_dmg=XXnew.F_pl_dmg[iPlaneGlob];
				XXnew.Bmean_pl_dmg+=B_norm;//after end of cycle on grains /=k_zer
				XXnew.Fmean_pl_dmg+=XXnew.F_pl_dmg[iPlaneGlob];//after end of cycle on grains /=k_zer

					//if(sig1.C(2,2)>400e6 && 
					//	(XXnew.F_pl_dmg[iPlaneGlob] > 0.0 || fabs(XXnew.B1_pl_dmg[iPlaneGlob]) > 0.0 || fabs(XXnew.B2_pl_dmg[iPlaneGlob]) > 0.0))
					//{
					//	printf("\n Br510: iPlaneGlob =%i, \n hBALoc(2,0)=%lg, hBALoc(2,1)=%lg, \nF_pl_0=%lg, B1_0=%lg, B2_0=%lg, \nF_pl_dmg=%lg, B1_pl_dmg=%lg, B2_pl_dmg=%lg"
					//				,iPlaneGlob, hBALoc(2,0), hBALoc(2,1), 
					//				F_pl_0, B1_0, B2_0,XXnew.F_pl_dmg[iPlaneGlob], XXnew.B1_pl_dmg[iPlaneGlob],XXnew.B2_pl_dmg[iPlaneGlob]);
					//	pause();
					//}

			//---------Control of micro-failure condition due to athermal plastic deformation -------------------------
    			Cond_dmg_1 = (1+MC.A_F*XXnew.F_pl_dmg[iPlaneGlob])
								*(MC.K1_F*SpSig + (1.0 + MC.K2_F*fabs(XXnew.B1_pl_dmg[iPlaneGlob])/MC.be_zv)*fabs(TauPlane(2,0)))
								-MC.Tau_F;
				Cond_dmg_2 = (1+MC.A_F*XXnew.F_pl_dmg[iPlaneGlob])
								*(MC.K1_F*SpSig + (1.0 + MC.K2_F*fabs(XXnew.B2_pl_dmg[iPlaneGlob])/MC.be_zv)*fabs(TauPlane(2,1)))
								-MC.Tau_F; 

				//if(sig1.C(2,1)>338.0e6 && Cond_dmg_1>0.0){
				//	printf("\n Br110: izer = %i, Cond_dmg_1=%lg, F_pl_dmg=%lg, \nB1_pl_dmg=%lg, hBeta31=%lg"
				//		        ,izer	,Cond_dmg_1   ,XXnew.F_pl_dmg[iPlaneGlob], XXnew.B1_pl_dmg[iPlaneGlob],hBALoc(2,0));
				//	pause();
				//}
				//if(sig1.C(2,1)>338.0e6 && Cond_dmg_2>0.0){
				//	printf("\n Br120: izer = %i, Cond_dmg_2=%lg, F_pl_dmg=%lg, \nB2_pl_dmg=%lg, hBeta32=%lg"
				//		        ,izer	,Cond_dmg_2   ,XXnew.F_pl_dmg[iPlaneGlob], XXnew.B2_pl_dmg[iPlaneGlob],hBALoc(2,1));
				//	pause();
				//}

				//Real Pi1_0 = XXold.Pi1_pl_dmg[iPlaneGlob];
				//Real Pi2_0 = XXold.Pi2_pl_dmg[iPlaneGlob];
				//Real hPi1 = 0.0;
				//Real hPi2 = 0.0;
				//if(Pi1_0 == 0.0 && Cond_dmg_1>0.0) {
				//	hPi1=fabs(XXnew.B1_pl_dmg[iPlaneGlob]);
				//}
			 //   if(Pi1_0>0.0)
				//{
				//	if(XXnew.B1_pl_dmg[iPlaneGlob]>0.0 && hB1_pl_dmg>0.0)		hPi1 = hB1_pl_dmg;
				//	else if(XXnew.B1_pl_dmg[iPlaneGlob]<0.0 && hB1_pl_dmg<0.0) hPi1 = -hB1_pl_dmg; 			  
				//}
				//
				//if(Pi2_0 == 0.0 && Cond_dmg_2>0.0){
				//	hPi2=fabs(XXnew.B2_pl_dmg[iPlaneGlob]);
				//}
            // 	if(Pi2_0>0.0){
				//	if(XXnew.B2_pl_dmg[iPlaneGlob]>0.0 && hB2_pl_dmg>0.0)		hPi2 = hB2_pl_dmg; 
				//	else if(XXnew.B2_pl_dmg[iPlaneGlob]<0.0 && hB2_pl_dmg<0.0)  hPi2 = -hB2_pl_dmg; 			  
				//}       		
				//XXnew.Pi1_pl_dmg[iPlaneGlob] = Pi1_0 + hPi1; 
				//XXnew.Pi2_pl_dmg[iPlaneGlob] = Pi2_0 + hPi2; 

				//if(Cond_dmg_1 > 0.0 || Cond_dmg_2 > 0.0){
				//	printf("\n Br200: XXnew.Pi1_pl_dmg=%lg, XXnew.Pi2_pl_dmg=%lg"
				//		        ,XXnew.Pi1_pl_dmg[iPlaneGlob],XXnew.Pi2_pl_dmg[iPlaneGlob]);
				//	pause();
				//}
				
			}// end if(TTau>0.0)
		   
			if(XXold.PI_dmg[izer] == 0.0)
			{
				//if(XXnew.Pi1_pl_dmg[iPlaneGlob] >= MC.Pi_cr || XXnew.Pi2_pl_dmg[iPlaneGlob] >= MC.Pi_cr)
				if(Cond_dmg_1 > 0.0 || Cond_dmg_2 > 0.0)
				{
					XXnew.PI_dmg[izer] = 1.0;
					//printf("\nFracture in grain %i occurred due to active plastic def-n", izer);
					//if(Cond_dmg_1 > 0.0){
					//	printf("\n Br110: Cond_dmg_1=%lg, F_pl_dmg=%lg, \nB1_pl_dmg=%lg, hBeta31=%lg"
					//				,Cond_dmg_1   ,XXnew.F_pl_dmg[iPlaneGlob], XXnew.B1_pl_dmg[iPlaneGlob],hBALoc(2,0));
					//	pause();
					//}
					//if(Cond_dmg_2 > 0.0){


					//	printf("\n Br120: Cond_dmg_2=%lg, F_pl_dmg=%lg, \nB2_pl_dmg=%lg, hBeta32=%lg"
					//				,Cond_dmg_2   ,XXnew.F_pl_dmg[iPlaneGlob], XXnew.B2_pl_dmg[iPlaneGlob],hBALoc(2,1));
					//}
					//pause();
				}
			}

        }
		//end for(iPlane)
      }
	  //end for(iSlipType)
      hbzapp = hbzap.Povernut(VVERH,rzer);
      //hbzap is the increment of plastic deformation in a grain
    }
    else{
      hbzapp = 0.0;
    }
//-------------------------------------------------------------------------
//-------- End of the calculation of athermal plastic deformatiom ---------
//-------------------------------------------------------------------------
     //printf("\nReo_800");
     //pause();

// Begin calculation of grain strain increment, new fraction of martensite
// and defects accumulation.
    FazerNew = 0.0;
    hbzph = 0.0; //Matrix to sum up PHASE strains for all systems
    hbzmp = 0.0; //Matrix to sum up MICROPLASTIC strains for all systems
    for(isis = 0; isis < k_sis; isis++) {
      int isisglob=smesh+isis;
      Real hFloc_total = hFloc[isis] + hFlocTw[isis];
      *(XXnew.Floc_dp + isisglob) = *(XXold.Floc_dp + isisglob)
				       + hFloc_total;
      FazerNew += *(XXnew.Floc_dp + isisglob);

      Real hFlocmp = hFi_MP[isis] + hFi_MP_Tw[isis];
      XXnew.mpDef_dp[isisglob] = XXold.mpDef_dp[isisglob] + hFlocmp;

    // Calculation of new defect parameters
      Real B0 = XXold.B_dmg[isisglob];
      Real F0 = XXold.F_dmg[isisglob];
      XXnew.B_dmg[isisglob] = B0 + hFlocmp*(1.0 - Hev2(B0,hFlocmp)*B0/MC.be_zv );
      XXnew.F_dmg[isisglob] = F0 + (1.0+MC.q3_F*F0)*fabs(B0)* fabs(hFlocmp) - r_dmg*F0* htim;
      F_dmg_glob += XXnew.F_dmg[isisglob];
      // Real r_dmg = MC.r0_F*exp(-MC.U_F/(RBolts*Tem)); - defined in the beginning of this program

     //Prepare for check of micro-failure condition
     Real Tau_loc_forvariant = tauzer.Svertka(Vrnts.Dev_var[isis])/MC.Gam_D;
     Real Concentration_factor = MC.K2_F*fabs(XXnew.B_dmg[isisglob])/MC.be_zv;

     //	if(XXnew.B_dmg[isisglob] > 0.0){
     //	  printf("\n B_dmg =%9.5lg" ,XXnew.B_dmg[isisglob]);
     //     pause();
     //   }
     //if(Concentration_factor > 0.0){
     //     printf("\n SpSig =%9.5lg Tau_loc=%9.5lg Concentration_factor=%9.5lg Tau_F_=%9.5lg"
     //		,SpSig     ,Tau_loc_forvariant, Concentration_factor   ,MC.Tau_F/(1.0+MC.A_F*XXnew.F_dmg[isisglob]));
     //     pause();
     //   }


     //Check the condition of micro-failure
     if(XXold.pi_dmg[isisglob]==0.0){
      //	if(MC.K1_F*SpSig + Tau_loc_forvariant * Concentration_factor >
      //			      MC.Tau_F/(1.0+MC.A_F*XXnew.F_dmg[isisglob])){//condition of micro-failure
			if(MC.K1_F*SpSig + Tau_loc_forvariant * Concentration_factor >
					MC.Tau_F/(1.0+MC.A_F*F_dmg_glob)) //condition of micro-failure
			{
					XXnew.pi_dmg[isisglob]=1.0;
			}
			else XXnew.pi_dmg[isisglob]=0.0;
      }
     //end if XXold.pi_dmg[isisglob]==0.0
	 else XXnew.pi_dmg[isisglob]=1.0;

     //Register the failure of the current grain [izer] (find XXnew.PI_dmg)
     //before evaluating XXnew.PI_dmg[izer] it was made equal  XXold.PI_dmg[izer]
      if(XXnew.PI_dmg[izer] == 0.0){
		if(XXnew.pi_dmg[isisglob]>0.0){
			XXnew.PI_dmg[izer] = 1.0;
			printf("\nFracture in grain %i occurred due to microplastics", izer);
			//pause();
		}
      }

    //Calculation of the new yield limit for microplastic flow.
    //Its recovery by the reverse transform-n is described by the second term
	Real F_Yield_n = XXold.Fyloc_dp[isisglob];
	XXnew.Fyloc_dp[isisglob] = F_Yield_n + Hm*fabs(hFlocmp)
	   -MC.Rm*XHeavisideX(F_Yield_n-MC.FYield)*XHeavisideX(-hFloc[isis]);
       //  printf("\n %6.1lg,",MC.Rm*XHeavisideX(F_Yield_n-MC.FYield)*XHeavisideX(-hFloc[isis]));
       //  pause();

       //}
//    }
      hbetp = Vrnts[isis];
      //hbetp = (hFloc[isis] + hFlocmp) ;
      hbetp.Multiplication(hFloc_total);
      hbzph += hbetp; //includes only PHASE strain (including twinning)

      hbetp = Vrnts[isis];
      hbetp.Multiplication(hFlocmp);
      hbzmp += hbetp;
    }//for isis
    
	//...........if(FazerNew > 0.0 || Fazer > 0.0) {
     FazaNew += (FazerNew /= k_sis);
     if(FazerNew<0.0 || FazerNew>1.0+1e-9){
       fprintf(stderr,"\nBad value of FazerNew=%30.28le",FazerNew);
       exit(2);
     }
     hbzph.Multiplication(1.0 / k_sis);
     hbzmp.Multiplication(MC.Kap2002 / k_sis);

     hbzphp= hbzph.Povernut(VVERH,rzer);
     hbzmpp= hbzmp.Povernut(VVERH,rzer);

     temp= hbzphp+hbzmpp;
     temp.Multiplication(MC.ScalingFactor);

     XXnew.EpsUnGr[izer] = XXold.EpsUnGr[izer] + temp;
     XXnew.EpsUnGr[izer] += hbzapp;
         // !!! if heat expansion is anisotropic we must add hear
         // heat expansion strain in a grain to XXnew.EpsUnGr[zer]
	 dEPh+= hbzphp;
     dEmp += hbzmpp;
     dEAP += hbzapp;


  //...........}

   }//for izer - close cycle on grains ----------------------------
   //==============================================================
     //printf("\nReo_900");
     //pause();


   XXnew.Phase = (FazaNew /= k_zer);
   dEPh.Multiplication(MC.ScalingFactor / k_zer);
   dEmp.Multiplication(MC.ScalingFactor / k_zer);
   dEAP.Multiplication(1.0 / k_zer);
   if (sig(2, 2)>400e6){
	   dEPh = 0.0;

   }
   XXnew.PhaseStrain = XXold.PhaseStrain + dEPh;
   XXnew.mpStrain = XXold.mpStrain + dEmp;
   XXnew.PlasticStrain = XXold.PlasticStrain + dEAP;


	XXnew.Bmean_pl_dmg/=k_zer;
	XXnew.Fmean_pl_dmg/=k_zer;
   //----------------------
   // Determine dETE
   dETE = 0.0;
   if(fabs(hTem) > 0.0) {
     cx =  MC.TECoefA * hTem * (1.0 - XXold.Phase)
         + MC.TECoefM * hTem * XXold.Phase
         +(MC.TECoefM-MC.TECoefA)*(Tem-MC.Ms)*(FazaNew-XXold.Phase);
     //cx = MC.TECoefA * hTem + (MC.TECoefM - MC.TECoefA)
        //* (FazaNew * hTem + (Tem-MC.Mf)
        //* (FazaNew-XXold.Phase));
        for(i=0; i < 3; i++) dETE(i,i) = cx;
   }//if
   XXnew.TEStrain = XXold.TEStrain + dETE;
   //-----------------------

   //Strain of heat expansion must be also added to XXnew.EpsUnGr[izer]
   for (izer=0;izer<k_zer;izer++){
     XXnew.EpsUnGr[izer] += dETE;
   }

   //Count the number of destroyed grains
   XXnew.kzer_dmg = 0;
   for(izer=0;izer<k_zer;izer++){
     if(XXnew.PI_dmg[izer]>0.0) XXnew.kzer_dmg++;  
   }
   if(XXnew.kzer_dmg < XXold.kzer_dmg){
     fprintf(stderr,"\n Reo: New number of destroyed grains %i < old %i"
                   ,XXnew.kzer_dmg,XXold.kzer_dmg);
     pause();
   }

   // Calculation of max b_F on grains
   //int izer_maxbf=0;
   //Real max_bfzer=0.0;//for maximum on all grains
   //for (izer=0;izer<k_zer;izer++){
   //  Real max_bf=0.0;//for maximum on variants in each grain
   //  for(isis = 0; isis < k_sis; isis++) {
   //    int smesh=izer*k_sis;//grain numb. times numb. of syst. in a grain
   //    int isisglob=smesh+isis;
   //    Real bf = fabs(XXnew.B_dmg[isisglob])*XXnew.F_dmg[isisglob];
   //    if(bf>max_bf)max_bf=bf;
   //  }
   //  XXnew.bf_dmg[izer]=max_bf;
   //  if(max_bf>max_bfzer){
   //    max_bfzer=max_bf;
   //    XXnew.izer_bf_max=izer;
   //  }
   //}

   //Check the condition of macro-fracture
   if(XXnew.kzer_dmg >= MC.kzer_cr) XXnew.PI_macro = 1;
   else XXnew.PI_macro = 0;

   real siggg = sig(2, 2);
   // Determine Elastic strain
   if (phoenix == 1.0){
	   sig(2, 2) = 0.0;
   }
   XXnew.ElStrain = Komb(1.0-FazaNew, Hook(MC.YoungA, MC.PuasA, sig),
	   FazaNew, Hook(MC.YoungM, MC.PuasM, sig));

   sig(2, 2) = siggg;


   //------------------------------------------------------------------
   //Now A-rad_creep = 0.0 in the beginning of this file
   if(A_rad_creep!=0.0) hrad_strain = A_rad_creep*DeltaFluence;
   else hrad_strain=0.0;
   XXnew.UnelasticStrain = hrad_strain;
   XXnew.UnelasticStrain+=XXnew.PhaseStrain;
   XXnew.UnelasticStrain+=XXnew.mpStrain;
   XXnew.UnelasticStrain+=XXnew.TEStrain;
   XXnew.UnelasticStrain+=XXnew.PlasticStrain;

  //----------- open cycle on grains and calculate XXnew.RoGrain------
   ErrRo = 0.0;
   MaxRoNorma = 0.0;
   for (izer=0;izer<k_zer;izer++){
     Matr33 Tmp = XXnew.UnelasticStrain - XXnew.EpsUnGr[izer];
     //XXnew.RoGrain[izer] = HookS(MC.YoungA, MC.PuasA, Tmp);
     //XXnew.RoGrain[izer] *= MC.CorrectRoGr;
     RoGrzv[izer] = HookS(MC.YoungA, MC.PuasA, Tmp);
     RoGrzv[izer] *= MC.CorrectRoGr;
     Tmp = RoGrProb[izer] - RoGrzv[izer];
     Real ErrTmp = Tmp.Norma();
     if(ErrTmp > ErrRo) ErrRo = ErrTmp;
     Real RoNorma=RoGrProb[izer].Norma();
     if(MaxRoNorma<RoNorma) MaxRoNorma = RoNorma;

    // Matr33 Tmp = XXnew.UnelasticStrain - XXnew.EpsUnGr[izer];
    // XXnew.RoGrain[izer] = HookS(MC.YoungA, MC.PuasA, Tmp);
    // Tmp = XXnew.RoGrain[izer] *= MC.CorrectRoGr;
    // XXnew.hRoGrain[izer] = Tmp - XXold.RoGrain[izer];
    // Real RoNorma=Tmp.Norma();
    } //--------------------end cycle on grains------------------------
   //XXnew.RoGrainNorma = MaxRoNorma;

   if(ErrRo<AdmissibleErrRo) {
        RoGrainIsFound=1;
        for (izer=0;izer<k_zer;izer++){
          XXnew.RoGrain[izer] = RoGrProb[izer];
          XXnew.hRoGrain[izer] = RoGrProb[izer] - XXold.RoGrain[izer];
        }
        XXnew.RoGrainNorma=MaxRoNorma;
    break;//from cycle to determine RoGrain
   }
   if(IterRo==1){
     //printf("\n First: LamRo=%lg IterRo=%3d ErrRo=%lg MaxRoNorma=%lg"
       //           ,LamRo,IterRo,ErrRo*1e-6,MaxRoNorma*1e-6);
     LamRop=LamRo;
     for(izer=0;izer<k_zer;izer++){
       Matr33 Tmp = Komb(1.0-LamRo, RoGrProb[izer], LamRo, RoGrzv[izer]);
       RoGrProbp[izer]=RoGrProb[izer];
       RoGrProb[izer]=Tmp;
     }
     ErrRop=ErrRo;
     continue; //cycle to determine RoGrain
   }//if(IterRo==1)
   if(ErrRo<0.60*ErrRop){
     //printf("\n  Good: LamRo=%lg IterRo=%3d ErrRo=%lg MaxRoNorma=%lg"
       //          ,LamRo,IterRo,ErrRo*1e-6,MaxRoNorma*1e-6);
     LamRop=LamRo;
     for(izer=0;izer<k_zer;izer++){
       Matr33 Tmp = Komb(1.0-LamRo, RoGrProb[izer], LamRo, RoGrzv[izer]);
       RoGrProbp[izer]=RoGrProb[izer];
       RoGrProb[izer]=Tmp;
     }
     ErrRop=ErrRo;
     continue; //cycle to determine RoGrain
   }//if(ErrRo<0.60*ErrRop)
   if(ErrRo<ErrRop){ //but ErrRo>0.80*ErrRop
     //printf("\n So-so: LamRo=%lg IterRo=%3d ErrRo=%lg MaxRoNorma=%lg"
       //           ,LamRo,IterRo,ErrRo*1e-6,MaxRoNorma*1e-6);
     LamRop=LamRo;
     if(LamRo>0.02) LamRo/=ROOT;
     else LamRo = 0.9;
     //else LamRo = 0.4+random(4000)/10000.0;
     //else LamRo = 1e-5 + random(800)/1000.0;
     //LamRo = 1e-6 + random(800)/1000.0;
     for(izer=0;izer<k_zer;izer++){
       Matr33 Tmp = Komb(1.0-LamRo, RoGrProb[izer], LamRo, RoGrzv[izer]);
       RoGrProbp[izer]=RoGrProb[izer];
       RoGrProb[izer]=Tmp;
     }
     ErrRop=ErrRo;
     continue; //cycle to determine RoGrain
   }//else if(ErrRo<ErrRop)
   if(ErrRo>ErrRop){
     //printf("\n   Bad: LamRo=%lg IterRo=%3d ErrRo=%lg MaxRoNorma=%lg"
       //           ,LamRo,IterRo,ErrRo*1e-6,MaxRoNorma*1e-6);
     if(LamRo>0.02) LamRo/=ROOT;
     //else LamRo = 1e-5 + random(800)/1000.0;
     else LamRo = 0.9;
     for(izer=0;izer<k_zer;izer++){
       Matr33 Tmp = Komb(1.0-LamRo/LamRop, RoGrProbp[izer], LamRo/LamRop, RoGrProb[izer]);
       RoGrProb[izer]=Tmp;
     }
     if(IterRo>20){
       for(izer=0;izer<k_zer;izer++){
         //Real Stir= 0.25+random(5000)/10000.0;
         //for(int ii=0;ii<3;ii++) for(int jj=0;jj<3;jj++){
           //Real Stir= 0.0001+random(10000)/10000.0;
           //RoGrProb[izer].C(ii,jj)*=Stir;
         //}
         RoGrProb[izer]=XXold.RoGrain[izer];
         LamRo=0.9/My_pow(ROOT,IAttempt);
         //if(IAttempt==1) RoGrProb[izer]=XXold.RoGrain[izer];
         //else{
           //Real Stir= 0.5+random(10000)/10000.0;
           //LamRop = -M_PI/10.0 + random(10000)/10000.0;
           //RoGrProb[izer] = Komb(1.0, XXold.RoGrain[izer]
           //                   , LamRop, XXold.hRoGrain[izer]);
           //RoGrProb[izer]*=Stir;
         //}
         //RoGrProbp[izer]*=Stir;
         //ErrRop=1e100;
         IterRo=0;//this will start a new attempt to find RoGrain
       }
       //ErrRo=ErrRop=1.0e100;
       //printf("\n New attempt started");
       printf("\n Difficult to equilibrate RoGrain. Please, wait.\n");
       IAttempt++;
     }
     ErrRo=ErrRop;
     LamRop=LamRo;
   }//if(ErrRo>ErrRop)
   //printf("\n LamRo=%lg IterRo=%d ErrRo=%lg MaxRoNorma=%lg ",LamRo,IterRo,ErrRo,MaxRoNorma);
   //if(IterRo==5) FixPlasticVariants=1;
 }// end of for- or while- cycle equlibrating of RoGrain

 if(RoGrainIsFound==1 && ErrRo>0.0){
   //fprintf(stderr,"\n RoGrain is found. LamRo=%lg ErrRo=%lg IterRo=%3d MaxRoNorma=%lg\n"
		   //,LamRo ,ErrRo*1e-6,IterRo,MaxRoNorma*1e-6);
   //pause();
 }

 if(IAttempt==IATTEMPTMAX){
   fprintf(stderr,"\nCan not equilibrate RoGrain. ErrRo=%lgkPa\n",ErrRo*1e-3);
   //pause();
 }

 //printf("\nNumOfGrWithTwin=%u",NumOfGrWithTwin);
 //printf(" MaxRo=%lg \n "  ,MaxRoNorma)
 XXnew.dPhase=(XXnew.Phase-XXold.Phase)/htim;
 XXnew.TotalStrain = XXnew.UnelasticStrain;
 XXnew.TotalStrain += XXnew.ElStrain;
 XXnew.dTotalStrain=XXnew.TotalStrain;
 XXnew.dTotalStrain-=XXold.TotalStrain;
 //XXnew.dTotalStrain/=htim;
 XXnew.dTotalStrain.Division(htim);

 delete [] hFloc;
 delete [] hFlocTw;
 delete [] hFi_MP;
 delete [] hFi_MP_Tw;
 delete [] KMP_Tw;
 delete [] TwDir;
 delete [] F_mix;
 delete [] FTw_n;
 delete [] RoGrProb;
 delete [] RoGrProbp;
 delete [] RoGrzv;
 delete [] tw_dir_is_allowed;
 return;
}//end of March