// Simple C example program that performs a Metropolis Monte Carlo Simulation with "Kawasaki" dynamics of 2d Conserved-Order-parameter Ising model with 2 Temp on a Rectangular lattice with periodic boundary conditions.
// Version 2: Includes a External field in the simulation.
// by Yi-Kai Mo, originally written 10/2022, last modified Oct. 20, 2022
//
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and uses simple graphics to show the evolution of the lattice.
//
// Parameters:
// Lx - linear extent of the lattice at x direction.
// Ly - linear extent of the lattice at y direction.
// N = Lx*Ly total holes/Lattice size
// pn - particle numbers
// Tx - temperature (real number)
// Ty - temperature (real number)
// h - External field (real number)
// seed - random number seed (should be a [large] negative integer)
//
//

#include <stdio.h>	  // standard io library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	  // standard math library
#define NRANSI		  // needed for NR
#include "my_nrutil.h"    // include NR header files

long L, N, pn;
double Tx, Ty;
long seed;
int **A;		// dimension pointer to pointer to int (needed for NR imatrix subroutine)
double *Boltzx, *Boltzy;

void initA();
void plotA();
void printA();
int HeatBath(int i, int j);
int MCmove(int i, int j);
float ran2(long *idum);	  // typecast ran2

int main(int argc, char *argv[]){   // argc and argv used for command line input
  int rint, i, j;
  int accepted;
  long t, mcs;

  FILE *outfile;	// pointer to filename

  if (argc == 7){		// require argc be equal to number of command line entries
    L = atol(argv[1]);		// read long variable from command line spot 1 (Lattice site)
    pn = atol(argv[2]);   // read long variable  (particle numbers)
    Tx = atof(argv[3]);   // read double variable (temperature for x axis)
    Ty = atof(argv[4]);		// read double variable (temperature for y axis)
    seed = atol(argv[5]);
    mcs = atol(argv[6]);  // maximum MCS
  }
  else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: 2Tmodel.x L pn Tx Ty seed mcs \n");  // correct input syntax
    
    return 1;
  }
  
  Boltzx = dvector(-3,3);
  for (i=-3; i<=3; ++i) {Boltzx[i] = exp(-4.0*i/Tx);}
  Boltzy = dvector(-3,3);
  for (i=-3; i<=3; ++i) {Boltzy[i] = exp(-4.0*i/Ty);}

  N = L*L;					// number of total nodes
  A = imatrix(0,L-1,0,L-1);			// use NR subroutine to allocate memory for the array A[i][j], where i and j in range [0, N-1]
  initA();					// initialize lattice A
  printf("Initial configuration: \n");
  printA();

  // Random shuffle the lattice A
   for (t=1;t<=N*100000; ++t){
    rint = N*ran2(&seed);
    i = rint/L;
    j = rint%L;
    accepted = HeatBath(i,j);
  }
  printf("Random shuffle final state. Numbers of MCS %ld: \n", 100000);
  printA();

  // Simulation starts
  for (t=1; t<=mcs*N; ++t){
    rint = N*ran2(&seed);
    i = rint/L;
    j = rint%L;
    accepted = MCmove(i,j);
  }
  printf("Final configuration. Numbers of MCS %ld: \n", mcs);
  printA();

  free_imatrix(A,0,L-1,0,L-1);			// NR subroutine to free allocated memory 
  free_dvector(Boltzx,-3,3);
  free_dvector(Boltzy,-3,3);

  return 0;
}

int HeatBath(int i, int j){
  int k,l,dir;
  int ip,im,jp,jm;
  int kp,km,lp,lm;
  int temp;
  dir = (int)4*ran2(&seed);
  if (dir==1){k=(i-1+L)%L;l=j;}      // dir == 0 up
  else if (dir==2){k=(i+1)%L;l=j;}   // dir == 1 down
  else if (dir==3){k=i;l=(j+1)%L;}   // dir == 2 right
  else{k=i;l=(j-1+L)%L;}             // dir == 3 left

  if (A[i][j]!=A[k][l])
  {
    temp = A[k][l];
    A[k][l] = A[i][j];
    A[i][j] = temp;
    return 1;}
  else return 0;

}

int MCmove(int i, int j){
  double roll;
  int k,l,dir;
  int ip,im,jp,jm;
  int kp,km,lp,lm;
  int ediff, temp;

  dir = (int)4*ran2(&seed);
  if (dir==0){k=(i-1+L)%L;l=j;}    // dir == 0 up
  else if (dir==1){k=(i+1)%L;l=j;} // dir==1 down
  else if (dir==2){k=i;l=(j+1)%L;} // dir == 2 right
  else{k=i;l=(j-1+L)%L;}           // dir == 3 left

  if (A[i][j]!=A[k][l])
  {
  // neighboor site for i,j
  ip = (i+1)%L;
  im = (i-1+L)%L;
  jp = (j+1)%L;
  jm = (j-1+L)%L;
  // neighboor site for k,l
  kp = (k+1)%L;
  km = (k-1+L)%L;
  lp = (l+1)%L;
  lm = (l-1+L)%L;
  
  // Energy difference for a swap
  ediff = ((A[ip][j] + A[i][jp] + A[im][j] + A[i][jm] - A[k][l])*A[i][j] + (A[kp][l] + A[k][lp] + A[km][l] + A[k][lm] - A[i][j])*A[k][l]);
  // Kawasaki dynamics with 2 temperature
  roll = ran2(&seed);

    if (dir<2)
    {
      if (roll<Boltzx[ediff/2]){
      temp = A[k][l];
      A[k][l] = A[i][j];
      A[i][j] = temp;
      return 1;
      }
      else return 0;
    }

    else
    {
      if (roll<Boltzy[ediff/2]){
      temp = A[k][l];
      A[k][l] = A[i][j];
      A[i][j] = temp;
      return 1;
      }
      else return 0;
    }
  }
  return 0;
}

void initA(){
  long i,j;
  float roll;
  for (i=0; i<L; ++i){
    for (j=0; j<L; ++j){ 
      if (j<pn) {A[i][j]=1; pn=pn-1;}
      else {A[i][j]=-1;}
  }
  }
  return;
  }


void plotA(){
  long i,j;
  for (i=0; i<L; ++i){
    for (j=0; j<L; ++j){
      if (A[i][j] == -1) fprintf(stderr,"0 ");
      else fprintf(stderr,"1 ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n\n");

  return;
}

void printA(){
  long i,j;
  for (i=0; i<L; ++i){
    for (j=0; j<L; ++j){
      if (A[i][j] == -1) printf("0 ");
      else printf("1 ");
    }
    printf("\n");
  }
  return;
}

// below is a NR random number generator. It generated float numbers evenly over range [0,1)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */
