//diaproc.h

#ifndef _DIALOG_LIBSAM

#define _DIALOG_LIBSAM

#include <stdio.h>
#include "defglob.h"

#define MALBUF 16
#define MAX_KEY_LEN 80
#define MAX_SECTION_NAME_LEN 32
#define LINE_LEN 80
//#define fnameread_int(fc,id_name) fkeyread_int((fc)," id_name ",&id_name)
//#define fnameread_double(f,r) fkeyread_double(f, "r", &r)

struct DATA_SECTION
{
  char name[MAX_SECTION_NAME_LEN];
  long pos0;
  long pos1;
  //Constructor
  DATA_SECTION(void);
};

enum TASK {load, strain33, strain32, zhestkost33, zhestkost32 ,vibr32 ,vibr32_quake
			,load_2_elements ,vibr_2_elements 
			,preload ,add_magnetic_force ,vibr_trig_polynom ,idle};

enum ETAP_TYPE {etap,cycle,endcycle};

//enum Operation_type{add_magnetic_force ,trig_polynom_excitation};

struct ETAP
{
  ETAP_TYPE type; // enum ETAP_TYPE {etap,cycle,endcycle};
  int num;
  //FILE* etap_res_file;
};

struct CYCLE
{
  int activated;
  long k; //number of iterations
  long i;  //current iteration
  long pos0;//starting position in wzd-file, i.e. position after
  //Constructor setting all members to zero
  CYCLE(void);
};

 void pause(void);
 void pause(char* msg);
 void pause(int pause_factor);
 void Exit_error(const char* msg);
 void Exit_error(const char* msg, int n);
 void Exit_error(const char* msg, double x);
 void Exit_error(const char* msg, char* s);
 void Exit_error(const char* msg_with_int, const char* msg, int n);

 void clear_buffer(int buffer_length, char * buffer);
 int fexist(char* name);
 FILE* fopen_next_dat(char* file_name, char* fmode);

 int freadnextkey(FILE *fc, char* found_key);
 void ffindkey(FILE *fc, const char* sample_key);
    //This utility finds and reads the next key including '=' and then stops.
    //Spaces and tabs inside a key are ignored.
    //all characters from // until the end of the line are skipped as comments.

 void fkeyread_string(FILE *fc, const char* sample_key, char* res, int res_length);
 double fkeyread_double(FILE *fc, const char* sample_key, double *d);
 int fkeyread_int(FILE *fc, const char* sample_key, int *d);
 long fkeyread_long(FILE *fc, const char* sample_key, long *d);

 int freadnextsection(FILE *fc, DATA_SECTION& sect); //auxiliary utility
 int ffindsection(FILE *fc, const char* sample_name, DATA_SECTION& current_section);
 //int ffindsection(FILE *fc, const char* sample_name);//looks for section with name "sample_name" (0 - success)
 int ffindkey_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key);
 int fkeyread_word_in_section(FILE *fc, DATA_SECTION& sect ,const char* sample_key,char* res, int res_length);
 //three functions reading upon key "sample_key" from DATA_SECTION "sect" (return 0 on success)
 int fkeyread_double_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, double *d);
 int fkeyread_int_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, int *d);
 int fkeyread_long_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, long *d);

 //this routine reads the string until the end of the line
 //for example this string may contain several names separated by spaces
   //int fkeyread_names_in_section(FILE *fc, DATA_SECTION& sect ,const char* sample_key
   //                              ,char* many_names, int res_length);
 //then these names may be taken one by one from string many_names by routine:
   //void take_name(char* first_name, char* many_names);
 //the first name is cut from string many_names and pasted into string first_name

 //routine reads a sequence of chars until ' ','\t' or '\n' is met
 // and removes comments like //abc...
 void fread_word(FILE* fc, char* res, int res_length);

 int recognize_etap(ETAP& et, DATA_SECTION& sect);

 void skiplines(FILE *fc ,int k_line);
 double frd1d(FILE*, double*);
 int frd1i(FILE*, int*);
 int frd1u(FILE *fc, unsigned *k);
 void frd3i(FILE*, int*,int*,int*);
 char* frd1s(FILE*, char*);

Real my_square(Real x);
Real sqr(Real x);
Real cube(Real x);
Real XHeavisideX(Real x);
//These functions may not be made inline: they are not class-members

Real min(Real a, Real b);
Real min(Real a, Real b, Real c);
Real min(Real a, Real b, Real c, Real d);
Real max(Real a, Real b);
Real max(Real a, Real b, Real c);
Real max(Real a, Real b, Real c, Real d);
int is_numeral(char c);
int numeral_follows(FILE* fc);//returns 1 if next non-blank char is a numeral (in FILE fc)
void get_number(char* s);//procedure truncates string s leaving char-s that suit for a number
void remove_comment(char* s); //in the result string comment is removed
                              //,i.e. slash '/' if followed by another '/' is replaced by '\0'

//Function aliquot(x,y,eps) for positive x,y,eps
//returns 1 if exists N : fabs(x - N*y) < eps*y
//oterwise returns 0
int aliquot(double x, double y, double eps);

//function ceil_2 returns minimum n=2**k, such that n>=x
int ceil_2(double x);

//Function Hev = Heaviside(x)
//returns 1.0 if x>0 , 0.0 if x<0, and second argument if x==0
//by default Heaviside(0) is 0.0
double Hev(double x, double Hev_at_0=0.0);

//Function Hev2 = Heaviside(x*y) 
//returns 1.0 if x and y have the same sign and 0.0 otherwise
double Hev2(double x, double y);

//Functions My_pow(x,y)
//returns x**y when x>0 and returns 0 when x==0 && y>=0 
double My_pow(double x, double y);

//Function My_pow(x,n)
//returns x**n for integer n and double x
double My_pow(double x, int n);

//Function My_power(x,n)
//returns k**n for integer k and n>=0
int My_pow(int k, int n);

double My_pow(int x, double n);

double Newton (double fix, double thickness_n);//, double fi_angle);


#endif