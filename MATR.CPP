// Matr33.cpp

#include <stdio.h>
#include <math.h>
#include "defglob.h"
#include "matr.h"
// Methods of class Matr33

void Matr33::operator=(Matr33& M2)
{
  short i,j;
  for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j]=M2.M[i][j];
}

Matr33 Matr33::operator=(double d)
{
  short i;
  double *pa=M[0];
  for(i=0;i<9;i++) *pa++=d;
  return *this;
}

Matr33 Matr33::operator+(Matr33& M2)
{
  short i;
  Matr33 Res;
  double *pa=M[0], *pb=(M2.M)[0], *pres=(Res.M)[0];
  for(i=0;i<9;i++) *pres++ = *pa++ + *pb++;
  return Res;
}

Matr33 Matr33::operator-(Matr33& M2)
{
  short i;
  Matr33 Res;
  double *pa=M[0], *pb=(M2.M)[0], *pres=(Res.M)[0];
  for(i=0;i<9;i++) *pres++ = *pa++ - *pb++;
  return Res;
}

Matr33 Matr33::operator*(Matr33& M2)
{
  Matr33 Res;
  double d, *sa, *cb;
  short i, j, k;

  for(i=0;i<3;i++)
  { for(j=0;j<3;j++)
    { sa=M[i];
      cb=&(M2.M)[0][j];
      d=0.0;
      for(k=0;k<3;k++,sa++,cb+=3)
          d+=(*sa)*(*cb);
      (Res.M)[i][j]=d;
    }
  }
  return Res;
}

Matr33 Matr33::operator+=(Matr33& M2)
{
  short i;
  double *pa=M[0], *pb=(M2.M)[0];
  for(i=0;i<9;i++) *pa++ += *pb++;
  return *this;
}

Matr33 Matr33::operator-=(Matr33& M2)
{
  short i;
  double *pa=M[0], *pb=(M2.M)[0];
  for(i=0;i<9;i++) *pa++ -= *pb++;
  return *this;
}

Matr33 Matr33::operator*=(double d)
{
  short i;
  double *pa=M[0];
  for(i=0;i<9;i++) *pa++ *= d;
  return *this;
}

Matr33 Matr33::operator/=(double d)
{
  short i;
  double *pa=M[0];
  for(i=0;i<9;i++) *pa++ /= d;
  return *this;
}

Matr33 Matr33::Povernut(short tomicro, Matr33& R)
{
 double d, *r1, *r2;
 short i,j,p,q;
 Matr33 Res;
 if(tomicro)
 {for(i=0;i<3;i++)
  {r1=(R.M)[i];
   for(j=0;j<3;j++)
   { r2=(R.M)[j];
    d=0;
//  printf("a[%1d][%1d]=%g\",i,j,a[i][j]);
    for(p=0;p<3;p++)
    {
      for(q=0;q<3;q++) d+= r1[p] * r2[q] * M[p][q];
    }
    (Res.M)[i][j] = d;
   }
  }
 }
 else
 {for(i=0;i<3;i++)
  {r1=&(R.M)[0][i];
   for(j=0;j<3;j++)
   { r2=&(R.M)[0][j];
    d=0;
    for(p=0;p<3;p++)
    {
     for(q=0;q<3;q++) d+= *(r1+3*p) * (*(r2+3*q)) * M[p][q];
    }
    (Res.M)[i][j] = d;
   }
  }
 }
 return Res;
}

void Matr33::Division(double d)
{
 short i;
 double *pa=&M[0][0];
 for(i=0;i<9;i++) *pa = (*pa++) / d;
}
void Matr33::Multiplication(double d)
{
 short i;
 double *pa=&M[0][0];
 for(i=0;i<9;i++) *pa = (*pa++) * d;
}

void Matr33::Add_to_diagonal(double d)
{                
 for(int i=0;i<3;i++) M[i][i]+=d;
}

double Matr33::Svertka(Matr33& M2)
{
  short i;
  double *pa=M[0], *pb=(M2.M)[0];
  double d=0.0;
  for(i=0;i<9;i++) d+= (*pa++) * (*pb++);
  return(d);
}

double Matr33::Norma(void)
{
  short i;
  double *pa=M[0];
  double d=0.0;
  for(i=0;i<9;i++){
    d+= (*pa) * (*pa);
    pa++;
  }
  return(sqrt(d));
}

double Matr33::SDevInt(void) //Matrix should be a deporoustor
{
  short i;
  double *pa=M[0];
  double d=0.0;
  for(i=0;i<9;i++){d+= (*pa) * (*pa); pa++;}
  return sqrt(d/2.0);
}

double Matr33::EDevInt(void) //Matrix should be a deporoustor
{
  short i;
  double *pa=M[0];
  double d=0.0;
  for(i=0;i<9;i++){d+= (*pa) * (*pa); pa++;}
  return( sqrt(2.0*d) );
}

Real Matr33::SIntensivityOn3Plane(void)
{
  return(  sqrt( 0.5*(M[2][0]*M[2][0]+M[2][1]*M[2][1]) )   );
}

void Matr33::frd1m(FILE* fc)
{
  short i;
  double *pa=&M[0][0];
  for(i=0;i<9;i++)fscanf_s(fc,"%lg",pa++);
}

void Matr33::ftpr1m(FILE *fil, const char* inf)
{ short i;
  fprintf(fil,"%s\n",inf);
   for(i=0;i<3;i++)
    fprintf(fil,"%-25.17lE %-25.17lE %-25.17lE\n"
              ,M[i][0],M[i][1],M[i][2]);
}
void Matr33::fpr1m(FILE *fil, const char* inf)
{ short i;
  fprintf(fil,"%s\n",inf);
   for(i=0;i<3;i++)
    fprintf(fil,"%-15.7lE %-15.7lE %-15.7E\n"
              ,M[i][0],M[i][1],M[i][2]);
}
void Matr33::pr1m(const char* inf)
{ short i;
  printf("%s\n",inf);
   for(i=0;i<3;i++)
    printf("%-15.7lE %-15.7lE %-15.7E\n"
              ,M[i][0],M[i][1],M[i][2]);
}
Matr33 Komb(double d1, Matr33& M1, double d2, Matr33& M2)
{
  short i;
  Matr33 Res;
  double *pa=(M1.M)[0], *pb=(M2.M)[0], *pres=(Res.M)[0];
  for(i=0;i<9;i++) *pres++ = d1*(*pa++) + d2*(*pb++);
  return Res;
}
