#ifndef _GLOB_CONSTANTS_AND_TYPES

#define _GLOB_CONSTANTS_AND_TYPES

#define ITERMAX 50
#define PODBORMAX 80
#define ERRDPABS_E 5.0e-5
#define ERRDPREL_W 1.0e-2
#define ERRDPABS_W 1.0e6
#define SUBITERMAX 2
#define FILENAME_LENGTH 64
#define SEVERAL_FILENAME_LENGTH 128


#define EXP_POW_MAX  709.7
#define EXP_POW_MIN  -708.6
#define BIGINT 10000
#define M_PI 3.14159265
#define MAXDOUBLE 1e300
#define SUCCESS 0
#define FAILURE 1


//typedef unsigned int;
//typedef int int;
typedef double Real;
typedef double real;


const Real pi=    3.14159265359;
const real pi_2= 1.57079632679;
const Real two_pi=6.28318530718;
const Real pi4pi = 0.25*(pi/4.0 -1.0/pi);
const Real RBolts = 8.31; // J/(K*mol) Universal gas constant 
#endif