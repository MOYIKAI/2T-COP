// Simple C example program that performs finite size scaling
// by Yi-Kai Mo, originally written 11/2022
//
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and uses simple graphics to show the evolution of the lattice.
//
// Parameters:
// inputfile - One coloumn with order parameter
// outfile   - output Binder clumaute

#include <stdio.h>	  // standard io library
#include <string.h>   // standard string library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	    // standard math library
#define NRANSI		    // needed for NR
#include "my_nrutil.h"    // include NR header files

int main(int argc, char *argv[]){ 
  double **data; 
  double s01, ss01, s10, ss10;
  double L, Tx, Ty, Tc, tx, ty, v;
  double y, xx,xy;
  int i, lnum, dl;

  FILE *infile;   // pointer to inputfile
  FILE *outfile;	// pointer to outputfile

  if (argc == 7){
    fprintf(stderr,"\n Good Initializtion:\n");
    L = atof(argv[1]);
    Tx = atof(argv[2]);
    Ty = atof(argv[3]);
    Tc = atof(argv[4]);
    v = atof(argv[5]);
  }
  else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: Fine1.x L Tx Ty Tc v infile\n");  // correct input syntax 
    return 1;
  }
  

  if (infile = fopen(argv[6],"r"))
  {
    i = 0;
    lnum = 0;
    data = dmatrix(0, 10000000, 0, 2);

    // Read all lines in the file and save into data array
    while (fscanf(infile, "%lf, %lf, %lf \n", &data[i][0], &data[i][1], &data[i][2]) != EOF ){i++;lnum++;}
    fclose(infile);

    ty=(Ty-Tc)/Tc;
    tx=(Tx-Tc)/Tc;
    s01=0;
    ss01=0;
    s10=0;
    ss10=0;
    dl=0;
    for (i=0; i<=lnum; i++){
      if (i%100==0){
        s10 += data[i][1];
        ss10 += data[i][1]*data[i][1];
        s01 += data[i][2];
        ss01 += data[i][2]*data[i][2];
        dl++;
      }
    }

    s10 = s10/dl;
    ss10 = ss10/dl;
    s01 = s01/dl;
    ss01 = ss01/dl;
    y = (2 - ss10/(s10*s10));
    xy = pow(L,(1/v)) * ty;
    xx = pow(L,(1/v)) * tx;

    // Print the value
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", xx, xy, s10, s01, y, L, Tx, Ty, Tc, v);




    free_dmatrix(data,0,10000000,0,2);
  }
  else{ fprintf(stderr,"No such input file \n"); return 1;}

  return 0;
}