#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double VERBOSE=0;

#include "HEMcore.c"

int main(int argc, char **argv)
{
  double inc, Ie;
  sscanf(argv[1], "%s", runname);
  sscanf(argv[2], "%lf", &a1Re);
  sscanf(argv[3], "%lf", &inc);
  sscanf(argv[4], "%lf", &Ie);
  sscanf(argv[5], "%lf", &H0H1L0);

  inc *= (M_PI/180.);
  Ie *= (M_PI/180.);

  strcpy(filename, runname);
  strcat(filename, ".bc");

  compute_constants(a1Re, inc, Ie, H0H1L0);
  printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.8le %.16le\n",
	 C1L0, C2L0, C3L0, G0L0, L1L0, H0H12L0, ww0, HL0e0, LsL0, E0);

  write_constants(filename);
}

  

