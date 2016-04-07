//diaproc.cpp

#include <stdio.h>
#include <stdlib.h>    // atof
#include <math.h>
#include "diaproc.h"

//defined in "diaproc.h"
//#define MALBUF 16
//#define MAX_KEY_LEN 80
//#define MAX_SECTION_NAME_LEN 32
//#define LINE_LEN 80
CYCLE::CYCLE(void)
{
  activated=0;
  k=0L; //number of iterations
  i=0L;  //current iteration
  pos0=0L;//starting position in wzd-file, i.e. position after
}

DATA_SECTION::DATA_SECTION(void)
{
  clear_buffer(MAX_SECTION_NAME_LEN,name);
}

Real my_square(Real x){return(x*x);}
Real sqr(Real x){return(x*x);}
Real cube(Real x){return(x*x*x);}
Real XHeavisideX(Real x){ return( x>=0 ? x : 0.0 ); }
Real min(Real a, Real b){return(a<b ? a : b); }
Real min(Real a, Real b, Real c)
{
 Real m12 = a<b ? a : b;
 return (m12<c ? m12 : c);
}

     double Func(double radius_new, double fix, double thickness_n){
		//return sin(radius_new/2)-fix*radius_new;}
		 return radius_new*log(radius_new/(radius_new-thickness_n))-fix;
	 }
     //double Func1(double fix, double fi_angle_new){
     //   return cos(fi_angle_new/2)/2-fix;}
	 double Fp(double radius_new, double Ek, double fix, double thickness_n){
		 return (Func(radius_new + Ek / 2, fix, thickness_n) - Func(radius_new - Ek / 2, fix, thickness_n)) / Ek;
	 }


double Newton (double fix, double thickness_n){//, double fi_angle){
const double v = -1.0;
   const double w = 1.0;
   const double epsilon = 0.0000001;
  // double fix;
    //int main ()
    {
    double radius_new;
    system("cls");
    radius_new=v;
    double radius_newt;
    int count=0;
    do
    {
    count++;
    radius_newt=radius_new;
    {
        radius_new=-Func(radius_new, fix, thickness_n)/Fp(radius_new,epsilon, fix, thickness_n)+radius_new;
} 
}
   while(fabs(radius_new-radius_newt)>=epsilon);
   radius_new=fabs(radius_new);
   return radius_new;
   // std::cout<<"fi_new="<<fi_angle_new;  //если не комментировать эти строки, 
   //system("pause");                      //то появится окно с новой фи и всё 
   // std::cin.get();

	}}

Real min(Real a, Real b, Real c, Real d)
{
 Real m12 = a<b ? a : b;
 Real m34 = c<d ? c : d;
 return (m12<m34 ? m12 : m34);
}

Real max(Real a, Real b){return(a>b ? a : b); }
Real max(Real a, Real b, Real c)
{
 Real m12 = a>b ? a : b;
 return (m12>c ? m12 : c);
}

Real max(Real a, Real b, Real c, Real d)
{
 Real m12 = a>b ? a : b;
 Real m34 = c>d ? c : d;
 return (m12>m34 ? m12 : m34);
}


int is_numeral(char c)
{
  if(  c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||c=='5'||c=='6'||c=='7'||c=='8'
     ||c=='9'||c=='+'||c=='-'||c=='.'||c=='e'||c=='E'||c=='L'||c=='U') return 1;
  else return 0;
}

int numeral_follows(FILE* fc)
{//returns 1 if next non-blank char is a numeral (in FILE fc)
  int res;
  char c=' ';
  long pos=ftell(fc);
  while(c==' ' || c=='\t'|| c=='\n') c=fgetc(fc);//skip blank char-s
  if( c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||c=='5'||c=='6'||c=='7'||c=='8'
     ||c=='9'||c=='+'||c=='-') res=1;
  else res=0;
  fseek(fc,-1L,SEEK_CUR);
  return res; 
}

void get_number(char* s)
{//procedure truncates string s leaving char-s that suit for a number
  char *p = s;
  while(is_numeral(*p)) p++;
  *p='\0';
  return;
}

void remove_comment(char* s)
{//procedure places '\0' just before symbols // are encountered
  char*p = s;
  char c=' ';
  while( !(*p=='\0' || *p=='/'&& c=='/') ) c = *p++;
  if(*p!='\0') *(p-1)='\0';
  return;
}

//Function Hev = Heaviside(x)
//returns 1.0 if x>0 , 0.0 if x<0, and Hev_at_0 if x==0 (by default 0.0)
double Hev(double x, double Hev_at_0){
  if(x>0.0) return 1.0;
  if(x<0.0) return 0.0;
  return Hev_at_0;
}



//Function aliquot(x,y,eps) for positive x,y,eps
//returns 1 if exists N : fabs(x - N*y) < eps * y
//oterwise returns 0
int aliquot(double x, double y, double eps){
  if(x<=0.0 || y<=0.0 || eps<=0.0 )
    Exit_error("Arguments of function aliquot must be positive");
  double r = fmod(x,y), eps_y=eps*y;
  if(r<eps_y || y - r <eps_y) return 1;
  else return 0;
}

//function ceil_2 returns minimum n=2**k, such that n>=x
int ceil_2(double x){
  if(x<=0.0)
    Exit_error("Argument of function ceil_2 must be positive");
  int n=1;
  while(n<x) n*=2;
  return n;
}

//Function Hev2 = Heaviside(x*y)
//returns 1.0 if x and y have the same sign and 0.0 otherwise
double Hev2(double x, double y){
  if(x>0.0 && y>0.0 || x<0.0 && y<0.0) return 1.0;
  return 0.0;
}

//Function My_pow(x,y)
//returns x**y when x>0 and returns 0 when x==0 && y>=0 
double My_pow(double x, double y){
  if(x<0){fprintf(stderr,"\nMy_pow: x<0 in x**y for double y"); pause(); exit(2);}
  if(x==0){
    if(y<0){fprintf(stderr,"\nMy_pow: x==0 && y<0 in x**y"); pause(); exit(2);}
    return 0.0;
  }
  double exp_pow = y*log(x);
  if(exp_pow>EXP_POW_MAX){ // define EXP_POW_MAX  709.7
    fprintf(stderr,"\nMy_pow: overflow in x**y:");
    fprintf(stderr,"\n x=%le, y=%le", x,y);
    pause(); exit(2);
  }
  if(exp_pow<EXP_POW_MIN){ // define EXP_POW_MIN  -708.6
    fprintf(stderr,"\nMy_pow: underflow in x**y:");
    fprintf(stderr,"\n x=%le, y=%le", x,y);
	pause(15);
    return 0.0;
  }
  return exp(exp_pow);
}

//Function My_pow(x,n)
//returns x**n for integer x and double n
double My_pow(int x, double n){
  return My_pow((double)x,n);
}

double My_pow(double x, int n){
  if(x==0){
    if(n<0){fprintf(stderr,"\nMy_pow: x==0 && n<0 in x**n"); pause(); exit(2);}
    else return 0.0;
  }
  double res=1.0;
  if(n>0) while(n-- > 0) res *= x;
  if(n<0) while(n++ < 0) res /= x;
  return res;
}

//Function My_pow(k,n)
//returns k**n for integer k and n>=0
int My_pow(int k, int n){
  if(n<0){
	fprintf(stderr,"\nInteger version of My_pow does not allow n<0 in k**n");
	pause("\nProgram will be terminated. "); 
	exit(2);
  }
  if(k==0){
    if(n==0){
		fprintf(stderr,"\n0**0 occured in function My_pow"); 
		pause("\nProgram will be terminated. "); 
		exit(2);
	}
    else return 0;
  }
  int res=1;
  while(n-- > 0) res *= k;
  //if(n<0) while(n++ < 0) res /= k; in this routine it must be n>=0
  return res;
}

 void pause(void)
{
 int n;
 fprintf(stdout, "\nPAUSE. Continue?(1/0)");
 fscanf_s(stdin,"%i",&n);
 //a=fgetc(stdin);
 //fgetc does not work correctly as well as fscanf(stdin,"%c",&a);
 //some other char-s may appear in stdin (when ENTER is pressed).

 //fprintf(stdout,"Answer is %i",n);
 if(n==0) exit(0);
 return;
}

 void pause(char* msg)
{
 int n;
 fprintf(stdout, "\nPAUSE. ");
 fprintf(stdout, "%s",msg);
 fprintf(stdout, "\nContinue?(1/0)");
 fscanf(stdin,"%i",&n);
 if(n==0) exit(0);
 return;
}

 void pause(int pause_factor)
{
  if(pause_factor==0) return;
  if(pause_factor<0){
    int n;
    fprintf(stdout, "\nPAUSE. Continue?(1/0)");
    fscanf(stdin,"%i",&n);
    if(n==0) exit(0);
    return;
  }
  if(pause_factor>0){
    int i,j,k;
    printf("\n PAUSE %i ",pause_factor);
    for(k=0;k<pause_factor;k++){
      printf("o");
      for(j=0;j<BIGINT;j++){
	for(i=0;i<10000;i++){}
      }
    }
    return;
  }
}

 void Exit_error(const char* msg)
 {
   fprintf(stderr,"\n\nERROR: ");
   fprintf(stderr,msg);
   fprintf(stderr,"\n Program is terminated.");
   pause(); exit(1);
 }

 void Exit_error(const char* msg, int n)
 {
   fprintf(stderr,"\n\nERROR: ");
   fprintf(stderr,msg,n);
   fprintf(stderr,"\n Program is terminated.");
   pause(); exit(1);
 }

 void Exit_error(const char* msg_with_int, const char* msg, int n)
 {
   fprintf(stderr,"\n\nERROR: ");
   fprintf(stderr,msg_with_int,n);
   fprintf(stderr,msg);
   fprintf(stderr,"\n Program is terminated.");
   pause(); exit(1);
 }

 void Exit_error(const char* msg, char* s)
 {
   fprintf(stderr,"\n\nERROR: ");
   fprintf(stderr,msg,s);
   fprintf(stderr,"\n Program is terminated.");
   pause(); exit(1);
 }

 void Exit_error(const char* msg, double x)
 {
   fprintf(stderr,"\n\nERROR: ");
   fprintf(stderr,msg,x);
   fprintf(stderr,"\n Program is terminated.");
   pause(); exit(1);
 }

 void clear_buffer(int buffer_length, char * buffer)
{
  int i;
  for(i=0;i<buffer_length-1;i++) buffer[i]=' '; buffer[i]='\0';
  return;
}

/*
void take_name(char* first_name, char* many_names)
{
  remove_comment(many_names);
  char *r=first_name, *p=many_names;
  while(*p==' ') p++;
  while(*p!=' ') *r++ = *p++;
  *r='\0';
  r=many_names;//now r points to beginning of many_names
  //and tail of many_names will be copied to its beginning
  while(*p!='\0')*r++ = *p++;
  *r='0';
  return;
}
*/

int fexist(char* name)
{
  FILE* ftemp=fopen(name,"r");
  if(ftemp==NULL) return 0;
  else{ fclose(ftemp); return 1;}
}

//this function opens file_name.dat for writing or appending
//file_name is only the name. suffix .dat is appended by this routine 
FILE* fopen_next_dat(char* file_name, char* fmode)
{
  int i;
  FILE *ftemp;
  ftemp=NULL;
  char buf[FILENAME_LENGTH];//FILENAME_LENGTH is defined in defglob.h
  //if file_name.dat does not exist we open this file and return
  sprintf(buf,"%s.dat",file_name);
  if(fopen(buf,"r")==NULL) return fopen(buf,fmode);
  //if file_name.dat already exists we open file file_name_???.dat
  //where ??? is a serial number
  clear_buffer(FILENAME_LENGTH,buf);
  int file_exist=1;
  for(i=0;i<999 && file_exist;i++){
    sprintf(buf,"%s_%03i.dat",file_name,i);
    //printf("\n buf=%s",buf); pause();
    if( file_exist=fexist(buf) ) clear_buffer(FILENAME_LENGTH,buf);
       //this assignment is correct
  }
  if(i==999) Exit_error("too many files %s...dat already exist",file_name);

    //printf("\n buf=%s",buf); pause();

  return fopen(buf,fmode);
}//end of routine fopen_next_dat

int freadnextkey(FILE *fc, char* found_key)
{//This utility finds and reads the next key including '=' and then stops.
 //Spaces and tabs inside a key are ignored.
 //found_key must contain '\0' only at its end
  char *p;
  char c=' '; //c - char from the stream.
  //int key_ready = 0;
  p=found_key;
  New_char: 
    if(*p=='\0'){
      fprintf(stderr,"\n\nERROR: too long word %s... in a data file",found_key);
	  pause(); exit(1);
    }
    c=fgetc(fc);
    if(c==EOF) return EOF;
    if(c=='\n'){p=found_key; goto New_char;}
    if(c==' '|| c=='\t') goto New_char;
    if(c=='#'){ while(c!='\n')c=fgetc(fc); p=found_key; goto New_char; } //skip comment

    if(c=='/'){ *p++=c; c=fgetc(fc);
      if(c=='/'){while(c!='\n'&&c!=EOF)c=fgetc(fc); p=found_key; goto New_char;} //skip comment
    }
    if(c != '='){*p++ = c; goto New_char;}
    *p='\0';
    //fprintf(stdout,"\nFound key %s",found_key);
    //pause();
    return 0;
}//freadnextkey

void ffindkey(FILE *fc, const char* sample_key)
{
  int key_is_found, key_is_wrong;
  char current_key[MAX_KEY_LEN];

  rewind(fc);

  
  key_is_found=0;
  while(clear_buffer(MAX_KEY_LEN,current_key),//buffer must be cleared before each iteration
        key_is_found==0 && freadnextkey(fc,current_key)!=EOF){
    //freadnextkey(fc,current_key);
    char *p_current_key;
    const char *p_sample_key;
    p_current_key=current_key;
    p_sample_key=sample_key;

    key_is_wrong=0;
    while(*p_sample_key !='\0' && key_is_wrong==0){
      while(*p_sample_key == ' ' || *p_sample_key == '\t') p_sample_key++;
      if(*p_current_key++ != *p_sample_key++) key_is_wrong=1;
    }
    //Now if key_is_wrong==0 then *p_sample_key must be =='\0'.
    //So for matching of the two keys also *p_current_key must be =='\0'.
    //Thus:
    if(*p_current_key !='\0') key_is_wrong=1;
    if(key_is_wrong==0) key_is_found=1;

    if(key_is_wrong==1){//skip a word belonging to this key
       char c=' ';
       while(c==' ' || c=='\t') c=fgetc(fc);//skip blank char-s
       while( !(c==' ' || c=='\t' || c=='\n' ) ) c=fgetc(fc);//skip not blank char-s
    }
  }//while

  if(key_is_found==0){
    fprintf(stderr,"\n\nERROR in a data file: key %s was not found",sample_key);
	pause(); exit(1);
  }
  return;
}//ffindkey

 void fkeyread_string(FILE *fc, const char* sample_key, char* res, int res_length)
{//in the result string comment is removed
 //,i.e. slash '/' if followed by another '/' is replaced by '\0'
 char c;
 char *p_res;
 clear_buffer(res_length,res);
 //find the key, which is needed:
 ffindkey(fc, sample_key);//now the last character read is '='
 //skip possible blank characters after '=' :
 c=' ';
 while(c==' ' || c=='\t' || c=='\n') c=fgetc(fc);
 //now c == the first character of the result string

 p_res = res;
 while(c!=' ' && c!='\t' && c!='\n' && c!=EOF){
   if(*p_res=='\0'){
     fprintf(stderr,"\nWARNING: string %s... in a data file is too long",res);
     break;
   }
   *p_res++ = c;
   c=fgetc(fc);
 }
 *p_res='\0';
 remove_comment(res);
 return;
}// fkeyread_string

 double fkeyread_double(FILE *fc, const char* sample_key, double *d)
{
  char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = strtod(buf, &endptr);
}//fkeyread_double

 int fkeyread_int(FILE *fc, const char* sample_key, int *d)
{
  //char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = atoi(buf);
}//fkeyread_int

 long fkeyread_long(FILE *fc, const char* sample_key, long *d)
{
  //char *endptr;
  char buf[MALBUF];
  fkeyread_string(fc,sample_key,buf,MALBUF);
  get_number(buf);
  return *d = atol(buf);
}//fkeyread_int

//---------------------------------------------------------------------------------
//functions for input from section
int freadnextsection(FILE *fc, DATA_SECTION& sect)
{//This utility finds the next section - expression like: [name] and assigns
 //sect.name = string inside [],
 //pos0 = position just after ']'
 //pos1 = position just before '[' of the next section
 // string section_name must contain '\0' only at its end
 //returns 0 on success and EOF if the end of file is reached

  //printf("\nfreadnextsection:");
  //printf("\nStart from position %li",ftell(fc));

  char *p;
  char c=' '; //c - char from the stream.
  long end_position_in_section = 0L;
  clear_buffer(MAX_SECTION_NAME_LEN,sect.name);
  p=sect.name;
  while(c!='['){
    c=fgetc(fc);
    if(c==EOF) return EOF;
    //---skip comment---------
   if(c=='/'){ c=fgetc(fc);
      if(c=='/'){while(c!='\n'&&c!=EOF)c=fgetc(fc);} //skip comment
   }
    //---skip comment---------
  }
  while(c!=']'){
    c=fgetc(fc);
    if(c==EOF) return EOF;
    if(*p=='\0'){
      fprintf(stderr,"\n\nERROR: too long section_name %s... ",sect.name);
	  pause(); exit(1);
    }
    if(c!=']') *p++ = c;
  }
  *p='\0';//specify the ending char. of sect.name
  sect.pos0 = ftell(fc);
  end_position_in_section = sect.pos0;
  while(c!='[' && c!=EOF){
	c=fgetc(fc);//looking for '[' opening next section
	if(c!=EOF) end_position_in_section++;
     //---skip comment---------
    if(c=='/'){ 
		c=fgetc(fc); if(c!=EOF) end_position_in_section++;
		if(c=='/'){
			while(c!='\n'&&c!=EOF){ 
				c=fgetc(fc); 
				if(c!=EOF) end_position_in_section++;
			}
		} //skip comment
    }
     //---skip comment---------
    //if(c==EOF) return EOF;
  }
  //sect.pos1 = ftell(fc)-1;
  sect.pos1 = end_position_in_section;
  //printf("\nfreadnextsection:");
  //printf("\nFound section [%s] with positions %li, %li",sect.name,sect.pos0,sect.pos1);
  //pause();
  return 0;
}//freadnextsection

int ffindsection(FILE *fc, const char* sample_name, DATA_SECTION& current_section)
{ //this utility looks for section, whose name is "sample_name"
  //blank char-s inside names are ignored
  //returns 1 if section is found and 0 if it is not found
  int section_is_wrong;
  rewind(fc);

  while(clear_buffer(MAX_SECTION_NAME_LEN,current_section.name),//name must be cleared before freadnextsection
        freadnextsection(fc,current_section)!=EOF){
    char *p_current_section;
    const char *p_sample_name;
    p_current_section=current_section.name;
    p_sample_name=sample_name;

    //compare sample_name and current_section.name :
    section_is_wrong=0;
    while(*p_sample_name !='\0' && section_is_wrong==0){
      while(*p_sample_name == ' ' || *p_sample_name == '\t') p_sample_name++;
      while(*p_current_section == ' ' || *p_current_section == '\t') p_current_section++;
      if(*p_current_section++ != *p_sample_name++) section_is_wrong=1;
    }
    //Now if section_is_wrong==0 then *p_sample_name must be =='\0'.
    //So for matching of the two names also *p_current_section.name must be =='\0'.
    //Thus:
    if(*p_current_section !='\0') section_is_wrong=1;
    if(section_is_wrong==0) break;//section is found and we stop reading sections
  }//while

    // printf("\nffindsection:");
    // printf("\ncurrent section name is %s",current_section.name);
    // if(section_is_wrong){ printf("\nSection with sample_name %s not found",sample_name);pause();}
    // else { printf("\nFound section with sample name %s",sample_name);pause();}

  if(section_is_wrong) return 1;
  else return 0; //return 0 on success
}//ffindsection

int ffindkey_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key)
{
  int key_not_found, key_is_wrong;
  char current_key[MAX_KEY_LEN];

  fseek(fc,sect.pos0,SEEK_SET);//go to the initial position of section "sect"

  key_not_found=1;
  while(clear_buffer(MAX_KEY_LEN,current_key)//here operator , is used
	,key_not_found && freadnextkey(fc,current_key)!=EOF && ftell(fc)<sect.pos1){
    //freadnextkey(fc,current_key);
    char *p_current_key;
    const char *p_sample_key;
    p_current_key=current_key;
    p_sample_key=sample_key;

    key_is_wrong=0;
    while(*p_sample_key !='\0' && key_is_wrong==0){
      while(*p_sample_key == ' ' || *p_sample_key == '\t') p_sample_key++;
      while(*p_current_key == ' ' || *p_current_key == '\t') p_current_key++;
      if(*p_current_key++ != *p_sample_key++) key_is_wrong=1;
    }
    //Now if key_is_wrong==0 then *p_sample_key must be =='\0'.
    //So for matching of the two keys also *p_current_key must be =='\0'.
    //Thus:
    if(*p_current_key !='\0') key_is_wrong=1;
    if(key_is_wrong==0) key_not_found=0;

    if(key_is_wrong==1){//skip a word belonging to this key
       char c=' ';
	   while((c==' ' || c=='\t') && c!=EOF) c=fgetc(fc); //skip blank char-s
       while( c!=' ' && c!='\t' && c!='\n' && c!=EOF) c=fgetc(fc);//skip non-blank char-s
	if(c==EOF) break;//from while-cycle
    }
  }//while

  if(key_not_found){
    //fseek(fc,sect.pos1,SEEK_SET);//go to the last position of section "sect"
    fseek(fc,sect.pos0,SEEK_SET);//go to the first position of section "sect"
    return 1;
  }
  //if key is found data after this key will be read by fkeyread_word_in_section
  else return 0;//return 0 on success
}// ffindkey_in_section

 int fkeyread_word_in_section(FILE *fc, DATA_SECTION& sect ,const char* sample_key
                                 ,char* res, int res_length)
{
 char c;
 char *p_res;
 clear_buffer(res_length,res);
 //find the key, which is needed:
 if( ffindkey_in_section(fc, sect, sample_key) ) return 1;//if there is no such key return 1
 //the last character read is '='
 //skip possible blank characters after '=' :
 c=' ';
 while(c==' ' || c=='\t' || c=='\n') c=fgetc(fc);
 //now c == the first character of the result string

 p_res = res;
 while(c!=' ' && c!='\t' && c!='\n'&& c!=EOF){
   if(*p_res=='\0'){
     fprintf(stderr,"\nWARNING: string %s... in a data file is too long",res);
	 pause();
     break;
   }
   *p_res++ = c;
   c=fgetc(fc);
 }
 *p_res='\0';
 remove_comment(res);
 //fseek(fc,sect.pos1,SEEK_SET);//go to the last position of section "sect"
 fseek(fc,sect.pos0,SEEK_SET);//go to the first position of section "sect"
 return 0;
}// fkeyread_word_in_section

/*
 this routine must be checked 
 int fkeyread_names_in_section(FILE *fc, DATA_SECTION& sect ,const char* sample_key
                                 ,char* res, int res_length)
 //routine reads the string until the end of the line
 //for example this string may contain several names separated by spaces
{
 char c;
 char *p_res;
 clear_buffer(res_length,res);//now the char at position #res_length is '\0'
 //find the key, which is needed:
 if( ffindkey_in_section(fc, sect, sample_key) ) return 1;//if there is no such key return 1
 //the last character read is '='
 //skip possible blank characters after '=' :
 c=' ';
 while(c==' ' || c=='\t' || c=='\n') c=fgetc(fc);
 //now c == the first character of the result string

 p_res = res;
 while(c!='\n'){
   if(*p_res=='\0'){
     fprintf(stderr,"\nWARNING: string %s... in a data file is too long",res);
     break;
   }
   *p_res++ = c;
   c=fgetc(fc);
 }
 *p_res='\0';
 remove_comment(res);
 //fseek(fc,sect.pos1,SEEK_SET);//go to the last position of section "sect"
 fseek(fc,sect.pos0,SEEK_SET);//go to the first position of section "sect"
 return 0;
}// fkeyread_names_in_section

*/

void fread_word(FILE* fc, char* res, int res_length)
{
 char c;
 char *p_res;
 clear_buffer(res_length,res);//now the char at position #res_length is '\0'
 p_res = res;
 c=fgetc(fc);
 while(c==' ' || c=='\t' || c=='\n') c=fgetc(fc);
 //now c == the first non-blank character of the result string

 while(c!=' ' && c!='\t' && c!='\n'&& c!=EOF){
   if(*p_res=='\0'){
     fprintf(stderr,"\nWARNING: string %s... in a data file is too long",res);
     pause();
     break;
   }
   *p_res++ = c;
   c=fgetc(fc);
 }
 *p_res='\0';
 remove_comment(res);
 //fprintf(stderr,"\nfread_word: word  %s  has been read", res);
 //pause();
 return;
}

 int fkeyread_double_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, double *d)
{
  char *endptr;//it is needed for calling strtod(...)
  char buf[MALBUF];
  if(fkeyread_word_in_section(fc,sect,sample_key,buf,MALBUF)) return 1;
  get_number(buf);
  *d = strtod(buf, &endptr);
  return 0; //returns 0 on success
}//fkeyread_double_in_section

 int fkeyread_int_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, int *d)
{
  char buf[MALBUF];
  if(fkeyread_word_in_section(fc,sect,sample_key,buf,MALBUF)) return 1;
  get_number(buf);
  *d = atoi(buf);
  return 0;
}//fkeyread_int_in_section

 int fkeyread_long_in_section(FILE *fc, DATA_SECTION& sect, const char* sample_key, long *d)
{
  char buf[MALBUF];
  if(fkeyread_word_in_section(fc,sect,sample_key,buf,MALBUF)) return 1;
  get_number(buf);
  *d = atol(buf);
  return 0;
}//fkeyread_long_in_section

int recognize_etap(ETAP& et, DATA_SECTION& sect)
{//returns 0 on success
  int i,j;
  char buf[MALBUF];
  char* p=sect.name;
  if( p[0]=='e'&& p[1]=='t' && p[2]=='a' && p[3]=='p' && p[4]==' ' ||
      p[0]=='E'&& p[1]=='T' && p[2]=='A' && p[3]=='P' && p[4]==' '
    ){
    et.type=etap;
    for(i=0,j=5 ; i<MALBUF-1 && p[j]!='\0' ; i++, j++) buf[i]=p[j];
    buf[i]='\0';
    et.num=atoi(buf);
    return 0;
  }

  if( p[0]=='c'&& p[1]=='y' && p[2]=='c' && p[3]=='l' && p[4]=='e' && p[5]==' ' ||
      p[0]=='C'&& p[1]=='Y' && p[2]=='C' && p[3]=='L' && p[4]=='E' && p[5]==' '
    ){
    et.type=cycle;
    for(i=0,j=6 ; i<MALBUF-1 && p[j]!='\0' ; i++, j++) buf[i]=p[j];
    buf[i]='\0';
    et.num=atoi(buf);
    return 0;
  }

  if( p[0]=='e'&& p[1]=='n' && p[2]=='d' && p[3]=='c' && p[4]=='y' && p[5]=='c' && p[6]=='l' && p[7]=='e' && p[8]==' ' ||
      p[0]=='E'&& p[1]=='N' && p[2]=='D' && p[3]=='C' && p[4]=='Y' && p[5]=='C' && p[6]=='L' && p[7]=='E' && p[8]==' '
    ){
    et.type=endcycle;
    for(i=0,j=9 ; i<MALBUF-1 && p[j]!='\0' ; i++, j++) buf[i]=p[j];
    buf[i]='\0';
    et.num=atoi(buf);
    return 0;
  }
  return 1;
}//recognize_etap

//old-style functions
 void skiplines(FILE *fc, int k_line)
{ char buf[LINE_LEN];
  while(k_line-- >0) fgets(buf,LINE_LEN,fc);
}

 double frd1d(FILE *fc, double *d)
{ char buf[LINE_LEN],buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,LINE_LEN,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t')) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  return *d = atof(buf1);
}
 int frd1i(FILE *fc, int *k)
{
  char buf[LINE_LEN], buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,LINE_LEN,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  return *k = atoi(buf1);
}

 int frd1u(FILE *fc, unsigned *k)
{
  char buf[LINE_LEN], buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,LINE_LEN,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  return *k = atoi(buf1);
}

void frd3i(FILE *fc, int *k1, int *k2, int* k3)
{
  char buf[LINE_LEN], buf1[MALBUF];//MALBUF = 16
  char *pb=buf, *pb1=buf1;
  int i=0;
  fgets(buf,LINE_LEN,fc);
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  *k1 = atoi(buf1);

  pb1=buf1;
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  *k2 = atoi(buf1);

  pb1=buf1;
  while(i<MALBUF && (*pb==' ' || *pb=='\t') ) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
  *k3 = atoi(buf1);
  return;
}

 char* frd1s(FILE *fc, char *s)
{ char buf[LINE_LEN];
  char *pb=buf, *pb1=s;
  int i=0;
  fgets(buf,LINE_LEN,fc);
//  printf("buf=%s",buf);
  while(i<MALBUF && (*pb==' ' || *pb=='\t')) { i++; pb++;}
  while(i<MALBUF && *pb!=' ' && *pb!='\t' ) { *pb1++ = *pb++; i++;}
  *pb1 = '\0';
//  printf("2:s=%s",s);
  return s;
}
