#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int VERBOSE=0;

double HL0min, HL0max;

#include "HEMcore.c"

int main(int argc, char **argv)
{
  double xp, yp, dt, nxp, nyp;
  double GL1;
  double G1L0;
  int flag;
  int npoints;

  sscanf(argv[1], "%s", runname);
  sscanf(argv[2], "%lf", &xp);
  sscanf(argv[3], "%lf", &yp);
  sscanf(argv[4], "%lf", &dt);
  if(argc == 6) {
    sscanf(argv[5], "%d", &npoints);
  } else {
    npoints = 1;
  }

  strcpy(filename, runname);
  strcat(filename, ".bc");
  read_constants(filename);

  strcpy(filenameplus, runname);
  strcat(filenameplus, ".plus");
  fileplus = fopen(filenameplus, "w");
  strcpy(filenameminus, runname);
  strcat(filenameminus, ".minus");
  fileminus = fopen(filenameminus, "w");
  strcpy(filenameerror, runname);
  strcat(filenameerror, ".error");
  fileerror = fopen(filenameerror, "w");

  /*
  printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
	 C1L0, C2L0, C3L0, G0L0, L1L0, H0H12L0, ww0, HL0e0, LsL0, E0);
   */

  GL1 = (xp*xp+yp*yp)/2.0;
  G1L0 = (1.0 - GL1)*L1L0;

  HL0min = MAX((-1.0*(G0L0 + H0H12L0)), (H0H12L0 - G1L0));
  HL0max = MIN((G0L0 - H0H12L0), (G1L0 + H0H12L0));

  if(VERBOSE) printf("minmax %.16le %.16le\n", HL0min, HL0max);

  flag = HEMmap(xp, yp, dt, npoints, &nxp, &nyp);

  if(npoints == 1) {
    printf("%.16le %.16le\n", nxp, nyp);
  }
}

  

