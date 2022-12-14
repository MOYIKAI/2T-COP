// Simple C example program that performs Binder clumaute
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
  double BS01, BS10; // Binder parameter for Structure factor (0,1) and (1,0)
  float Tx, Ty;
  int L, i, lnum, dl;

  FILE *infile;   // pointer to inputfile
  FILE *outfile;	// pointer to outputfile

  if (argc == 5){
    fprintf(stderr,"\n Good Initializtion:\n");
    L = atol(argv[1]);
    Tx = atof(argv[2]);
    Ty = atof(argv[3]);
  }
  else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: Binder1.x L Tx Ty infile\n");  // correct input syntax 
    return 1;
  }
  

  if (infile = fopen(argv[4],"r"))
  {
    i = 0;
    lnum = 0;
    data = dmatrix(0, 10000000, 0, 2);

    // Read all lines in the file and save into data array
    while (fscanf(infile, "%lf, %lf, %lf \n", &data[i][0], &data[i][1], &data[i][2]) != EOF ){
      i++;
      lnum++;
    }
    fclose(infile);
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
    BS01 = 1 - ss01/(3*s01*s01);
    BS10 = 1 - ss10/(3*s10*s10);

    // Print the Binder variable
    printf("%d, %f, %f, %lf, %lf, %lf, %lf \n", L, Tx, Ty, BS10, BS01, s10, s01);



    free_dmatrix(data,0,10000000,0,2);
  }
  else{ fprintf(stderr,"No such input file \n"); return 1;}

  return 0;
}