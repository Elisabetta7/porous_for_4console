// Reoproc.h
#ifndef _REO_SAM
#define _REO_SAM

#include <stdio.h>
#include "matr.h"
#include "inpar.h"

Real ThermodynamicForce(int //number of variant
                       ,VariantsOrientations& //VO
                       ,MaterialConstants& //MC
                       ,Real //Temperature
                       ,Matr33& //Sig
                       );

Matr33 Hook(
              double //E
             ,double //Nu
             ,class Matr33& //Sigma
           );
Matr33 HookS(Real E, Real Nu, Matr33 &Eps);

Real E_Young(Real phase, MaterialConstants& MC);

void CalculationOfHysterisisParameters
              (
               MaterialConstants& MC
              ,Real Fluence
              ,Real& T0
              ,Real& Ftr
              ,Real& Gm
              ,Real& Q0
              );
Real ThermodynamicForce(int n
                       ,VariantsOrientations& VO
                       ,MaterialConstants& MC
                       ,Real T
                       ,Matr33& Sig
                       );

Real ReducedGibbsForceForTwinning(int n
            ,VariantsOrientations& VO
            ,Matr33& sig);


void Reology(
            Real //DeltaTime
            ,double //Temperature
            ,double //DeltaTemperature
            ,Real Fluence
            ,Real hFluence
            ,class Matr33& // Sigma
            ,class Matr33& // DeltaSigma
			//,double //force
			//  ,double //hforce
            ,class MaterialConstants& //TMC
            ,class VariantsOrientations& //Vrnts
            ,class GrainsParameters& //GrPar
            ,class InternalParameters& //XXold
            ,class InternalParameters& //XXnew
            //,class Matr33& //dEPh
            //,class Matr33& //dETE
            //,class Matr33& //dEEl
			, Real // phoenix
          );
#endif