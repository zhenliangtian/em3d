typedef struct {
  double a1Re, H0H1L0; /* yes, redundant */
  double HL0e0, C1L0, C2L0, C3L0, G0L0, L1L0, H0H12L0, G0L1L0, ww0, LsL0, tscale, E0;
} BinaryConstants;

#include "integrator.h"
double INTEGRATION_EPSILON = 1.e-11;
double epssection = 1.e-10;

double HL0min, HL0max;

double a1Re, H0H1L0;
double HL0e0, C1L0, C2L0, C3L0, G0L0, L1L0, H0H12L0, G0L1L0, ww0, LsL0;
double tscale = 365.25*24.0*3600.;
double E0, Esimple;

/* variables used for HEMexplore */
double h_min, h_max;
double i_min, i_max;
double Ie_min, Ie_max;
int i_min_on_section, i_max_on_section;
int Ie_min_on_section, Ie_max_on_section;

char runname[100], filename[100], filenameplus[100], filenameminus[100], filenameerror[100];
FILE *fileplus, *fileminus, *fileerror;
char filenameoutput[100];
FILE *fileoutput;

void write_constants(char *runname)
{
  FILE *file;
  BinaryConstants bc;
  int flag;

  file = fopen(runname, "w");

  bc.a1Re = a1Re;
  bc.H0H1L0 = H0H1L0;

  bc.HL0e0 = HL0e0;
  bc.C1L0 = C1L0;
  bc.C2L0 = C2L0;
  bc.C3L0 = C3L0;
  bc.G0L0 = G0L0;
  bc.L1L0 = L1L0;
  bc.H0H12L0 = H0H12L0;
  bc.G0L1L0 = G0L1L0;
  bc.ww0 = ww0;
  bc.LsL0 = LsL0;
  bc.E0 = E0;
  bc.tscale = tscale;

  flag=fwrite((char *)&bc, sizeof(bc), 1, file);

  fflush(file);
  fclose(file);
}

void read_constants(char *runname)
{
  FILE *file;
  BinaryConstants bc;
  int flag;

  file = fopen(runname, "r");

  flag = fread((char *)&bc, sizeof(bc), 1, file);

  a1Re = bc.a1Re;
  H0H1L0 = bc.H0H1L0;

  HL0e0 = bc.HL0e0;
  C1L0 = bc.C1L0;
  C2L0 = bc.C2L0;
  C3L0 = bc.C3L0;
  G0L0 = bc.G0L0;
  L1L0 = bc.L1L0;
  H0H12L0 = bc.H0H12L0;
  G0L1L0 = bc.G0L1L0;
  ww0 = bc.ww0;
  LsL0 = bc.LsL0;
  E0 = bc.E0;
  tscale = bc.tscale;

  fclose(file);
}

#define NUMBER_EQUATIONS 4

double sgn(double x)
{
  if(x<0.0) {
    return(-1.0);
  } else
    return(1.0);
}      

double cubic(double a, double b, double c)
{
  double q, r, aa, bb, u, v;

  q = (a*a-3.0*b)/9.0;
  r = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  aa = -1.0*sgn(r)*pow(fabs(r)+sqrt(r*r - q*q*q), 1.0/3.0);
  bb = 0.0;
  if(aa != 0.0) bb = q/aa;
  u = aa+bb;
  v = -(1.0/3.0)*a;
  return(u+v);
}

double energy(Istate *s)
{
  double h, yp, HL0, xp;
  double H0L0, H1L0, GL1, G1L1, efactor, e2;
  double e2cos2w, cosI, cosi, sinI, sini, coseps, E;

  h = s->x[0];
  yp = s->x[1];
  HL0 = s->x[2];
  xp = s->x[3];

  H0L0 = H0H12L0 + HL0;
  H1L0 = H0H12L0 - HL0;
  GL1 = (xp*xp + yp*yp)/2.0;
  G1L1 = 1.0 - GL1;
  efactor = 1.0 - GL1;
  e2 = GL1*(2.0 - GL1);
  e2cos2w = ((xp*xp - yp*yp)/2.0)*(2.0 - GL1);
  /* printf("%.16le %.16le %.16le %.16le %.16le\n", H0H12L0, HL0, H0L0, G0L0, H0L0/G0L0); */
  cosI = H0L0/G0L0;
  cosi = H1L0/(G1L1*L1L0);
  sinI = sqrt(1.0 - cosI*cosI);
  sini = sqrt(1.0 - cosi*cosi);
  coseps = cosi*cosI + sini*sinI*cos(h);
  E = C1L0*(0.75*coseps*coseps - 0.25)/(efactor*efactor*efactor)
    + C2L0*(0.75*cosI*cosI - 0.25)
    + C3L0*((1.0+2.5*e2)*(0.75*cosi*cosi - 0.25)
	    + 1.875*(1.0 - cosi*cosi)*e2cos2w);

  /*
  printf("--\n");
  printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le\n", H0L0, H1L0, GL1, G1L1, efactor, e2, e2cos2w);
  printf("%.16le %.16le %.16le %.16le %.16le %.16le\n", cosI, cosi, sinI, sini, coseps, E);
    */

  return(E);
}

void compute_constants(double a1Re, double inc, double Ie, double H0H1L0) 
{
  double L1, c0, c1, c2, c3, C, G0, rootL1;

  double G=6.674e-11;
  double m0=5.97237e24;
  double m1=7.34612e22;
  double m2=1.9885e30;
  double Re=6378.1366e3;
  double a2=1.49597870e11;
  double J20=0.00108;
  double w0=2.0*M_PI/(24.0*3600.);
  double C0=0.33;

  double L0 = 0.3308*m0*Re*Re*sqrt(G*m0/(Re*Re*Re));
  double C10, C20, C30;
  C10 = -1.0*G*m0*m1*J20/Re;
  C20 = -1.0*G*m0*m2*J20*Re*Re/(a2*a2*a2);
  C30 = -0.5*G*m1*m2*Re*Re/(a2*a2*a2);

  L1 = m1*sqrt(G*m0*a1Re*Re);
  c0 = L1*cos(inc)/L0 -  H0H1L0;
  c1 = m0*Re*Re*C0*w0*cos(Ie)/L0;
  c2 = 0.0;
  c3= (2.0/3.0)*J20*m0*Re*Re*w0*cos(Ie)/L0;
  ww0 = cubic(c2/c3, c1/c3, c0/c3);
  /* printf("%.16le %.16le %.16le %.16le\n", c2/c3, c1/c3, c0/c3, ww0); */
  C = (C0 + J20*(2.0/3.0)*ww0*ww0)*m0*Re*Re;
  G0 = C*ww0*w0;
  rootL1 = sqrt(L1);
  if(fabs(H0H1L0 - (G0*cos(Ie) + L1*cos(inc))/L0) > 1.e-10) {
    printf("wrong H0+H1/L0 %.16le %.16le", H0H1L0, (G0*cos(Ie) + L1*cos(inc))/L0);
    exit(1);
  }
  HL0e0 = (G0*cos(Ie) - L1*cos(inc))/(2.0*L0);
  C1L0 = (C10*ww0*ww0/(a1Re*a1Re*a1Re))/L0;
  C2L0 = (C20*ww0*ww0)/L0;
  C3L0 = (C30*a1Re*a1Re)/L0;
  G0L0 = G0/L0;
  L1L0 = L1/L0;
  H0H12L0 = H0H1L0/2.0;
  G0L1L0 = (G0+L1)/L0;
  LsL0 = (G0/L0 + L1/L0);

  Istate s0;
  s0.t = 0.0;
  s0.x[0] = 0.0;
  s0.x[1] = 0.0;
  s0.x[2] = HL0e0;
  s0.x[3] = 0.0;
  E0 = energy(&s0);

  Esimple = (0.5*C*(ww0*w0)*(ww0*w0) - G*m0*m1/(2.0*a1Re*Re))/(0.5*0.3308*m0*Re*Re*w0*w0 - G*m0*m1/(2.0*60.0*Re));
  /*
  printf("Esimple: %.16le %.16le %.16le %.16le\n", 
	 0.5*C*(ww0*w0)*(ww0*w0), - G*m0*m1/(2.0*a1Re*Re), (0.5*0.3308*m0*Re*Re*w0*w0 - G*m0*m1/(2.0*60.0*Re)), Esimple);
    */
}

int equations_of_motion(Istate *state, Istate_Deriv *state_deriv)
{
  double h, yp, HL0, xp;
  double H0L0, H1L0, GL1, Gdotx, Gdoty, G1L1;
  double efactor, e2, e2cos2w, cosI, cosi;
  double cosIdotHL0, cosidotHL0, cosidotxL0, cosidotyL0;
  double e2dotxL0, e2dotyL0, efactordotxL0, efactordotyL0;
  double e2cos2wdotxL0, e2cos2wdotyL0, sinI, sini;
  double sinIdotHL0, sinidotHL0, sinidotxL0, sinidotyL0;
  double coseps, cosepsdotHL0, cosepsdotxL0, cosepsdotyL0;
  double cosepsdoth;

  h = state->x[0];
  yp = state->x[1];
  HL0 = state->x[2];
  xp = state->x[3];

  H0L0 = H0H12L0 + HL0;
  H1L0 = H0H12L0 - HL0;
  GL1 = (xp*xp + yp*yp)/2.0;
  Gdotx = xp;
  Gdoty = yp;
  G1L1 = 1.0 - GL1;
  efactor = 1.0 - GL1;
  e2 = GL1 * (2.0 - GL1);
  e2cos2w = ((xp*xp - yp*yp)/2.0)*(2.0 - GL1);
  cosI = H0L0/G0L0;
  cosi = H1L0/(G1L1*L1L0);
  if((cosI > 1.0) || (cosI < -1.0)) return(FAILURE);
  if((cosi > 1.0) || (cosi < -1.0)) return(FAILURE);
  cosIdotHL0 = 1.0/G0L0;
  cosidotHL0 = -1.0/(G1L1*L1L0);
  cosidotxL0 = cosi*(Gdotx/(G1L1*L1L0));
  cosidotyL0 = cosi*(Gdoty/(G1L1*L1L0));
  e2dotxL0 = (Gdotx*(2.0 - GL1) - GL1*Gdotx)/L1L0;
  e2dotyL0 = (Gdoty*(2.0 - GL1) - GL1*Gdoty)/L1L0;
  efactordotxL0 = -1.0*Gdotx/L1L0;
  efactordotyL0 = -1.0*Gdoty/L1L0;
  e2cos2wdotxL0 = (xp*(2.0 - GL1) - ((xp*xp - yp*yp)/2.0)*Gdotx)/L1L0;
  e2cos2wdotyL0 = (-yp*(2.0 - GL1) - ((xp*xp - yp*yp)/2.0)*Gdoty)/L1L0;
  sinI = sqrt(1.0 - cosI*cosI);
  sini = sqrt(1.0 - cosi*cosi);
  sinIdotHL0 = -1.0*cosI*cosIdotHL0/sinI;
  sinidotHL0 = -1.0*cosi*cosidotHL0/sini;
  sinidotxL0 = -1.0*cosi*cosidotxL0/sini;
  sinidotyL0 = -1.0*cosi*cosidotyL0/sini;
  coseps = cosi*cosI + sini* sinI*cos(h);
  cosepsdotHL0 = cosi*cosIdotHL0 + cosidotHL0*cosI + sinidotHL0*sinI*cos(h) + sini*sinIdotHL0*cos(h);
  cosepsdotxL0 = cosidotxL0*cosI + sinidotxL0*sinI*cos(h);
  cosepsdotyL0 = cosidotyL0*cosI + sinidotyL0*sinI*cos(h);
  cosepsdoth = -1.0*sini*sinI*sin(h);

  state_deriv->t = 1.0;
  state_deriv->x[0] = tscale*(C1L0*(1.5*coseps*cosepsdotHL0)/(efactor*efactor*efactor)
			     + C2L0*(1.5*cosI*cosIdotHL0)
			     + C3L0*((1.0 + 2.5*e2)*(1.5*cosi*cosidotHL0)
				     - 3.75*cosi*cosidotHL0*e2cos2w));
  state_deriv->x[1] = tscale*(C1L0*((1.5*coseps*cosepsdotxL0/(efactor*efactor*efactor))
				   + (0.75*coseps*coseps - 0.25)*(-3.0*efactordotxL0)/(efactor*efactor*efactor*efactor))
			     +C3L0*((2.5*e2dotxL0*(0.75*cosi*cosi - 0.25))
				    +(1.0 + 2.5*e2)*(1.5*cosi*cosidotxL0)
				    +(-3.75*cosi*cosidotxL0*e2cos2w)
				    +(1.875*(1.0 - cosi*cosi)*e2cos2wdotxL0)));
  state_deriv->x[2] = tscale*(-1.0*C1L0*1.5*coseps*cosepsdoth/(efactor*efactor*efactor));
  state_deriv->x[3] = -1.0*tscale*(C1L0*((1.5*coseps*cosepsdotyL0/(efactor*efactor*efactor))
					+ (0.75*coseps*coseps - 0.25)*(-3.0*efactordotyL0)/(efactor*efactor*efactor*efactor))
				  +C3L0*((2.5*e2dotyL0*(0.75*cosi*cosi - 0.25))
					 +(1.0 + 2.5*e2)*(1.5*cosi*cosidotyL0)
					 +(-3.75*cosi*cosidotyL0*e2cos2w)
					 +(1.875*(1.0 - cosi*cosi)*e2cos2wdotyL0)));
  return(SUCCESS);
}

double rooteps = 1.e-12;

double find_root(double xp, double yp, double divisions)
{
  Istate s0, s, smid;
  int idivisions;

  double xlow;
  double dx;

  xlow = HL0min;
  dx = (HL0max - HL0min)/divisions;

  /*  printf("xlow %.16le %.16le\n", xlow, dx); */

  if(VERBOSE) {
    double x;
    for(x=HL0min; x<HL0max; x+=(HL0max-HL0min)/100.) {
      s0.t=0.0;
      s0.x[0] = 0.;
      s0.x[1] = yp;
      s0.x[2] = x;
      s0.x[3] = xp;
      printf("%.16le %.16le\n", x, energy(&s0));
    }
    s0.t=0.0;
    s0.x[0] = 0.;
    s0.x[1] = yp;
    s0.x[2] = HL0max;
    s0.x[3] = xp;
    printf("%.16le %.16le\n", x, energy(&s0));
    /* there's an extra root near HL0max sometimes, why? */
  }

  for(idivisions=0; idivisions<((int) divisions); idivisions++) {

    double f, f0, fmid;
    double xhigh, xmid;

    s0.t=0.0;
    s0.x[0] = 0.;
    s0.x[1] = yp;
    s0.x[2] = xlow;
    s0.x[3] = xp;
    f0 = E0 - energy(&s0);
    s = s0;
    xhigh = xlow + dx;
    s.x[2] = xhigh;
    f = E0 - energy(&s);
    if(VERBOSE) printf("%.16le %.16le %.16le %.16le %.16le %.16le\n", xlow, xhigh, E0, energy(&s0), energy(&s), HL0max-xhigh);
    if(f0*f < 0.0) {
      /* refine */
      smid = s0;
      do {
	/* printf("%.16le %.16le %.16le\n", xlow, xhigh, xhigh-xlow); */
	xmid = (xlow + xhigh)/2.0;
	smid.x[2] = xmid;
	fmid = E0 - energy(&smid);
	if(f0*fmid < 0.0) {
	  f = fmid;
	  xhigh = xmid; 
	} else {
	  f0 = fmid;
	  xlow = xmid;
	}
      } while(fabs(xhigh - xlow) > rooteps);
      return((xhigh+xlow)/2.0);
    } else {
      xlow = xhigh;
    }
  }
  if(VERBOSE) printf("divisions %.16le\n", divisions);
  return(find_root(xp, yp, divisions*10.));
}

#include "integrator.c"

double pv0(double x)
{
  double TWOPI=(2.0*M_PI);
  double n;
  n = floor(x/(TWOPI));
  x -= TWOPI*n;
  if(x>M_PI) x -= TWOPI;
  return(x);
}

#define DTMAX 1.0

int integrate_until_success(Istate *s0, double *dt) 
{
  Istate s;
  Istate sb;
  double dts;

  s = *s0;
  if(fabs(*dt) > DTMAX) *dt = DTMAX*sgn(*dt);
  if(VERBOSE) printf("success dt: %.16le\n", *dt);
  dts = *dt;
  if(integrate(&s, dt) == SUCCESS) {
    if(VERBOSE) printf("success state: %.16le %.16le %.16le %.16le %.16le %.16le\n", s.t, pv0(s.x[0]), s.x[1], s.x[2], s.x[3], *dt);
    *s0 = s;
    /* the recommended dt is sometimes too small */
    if(fabs(*dt) < fabs(dts)) {
      *dt = dts;
    }
    return(SUCCESS);
  } else {
    if(VERBOSE) printf("until dt: %.16le\n", *dt);
    *dt /= 2.0;
    return(integrate_until_success(s0, dt));
  }
}

int refine_crossing(Istate *z0, double tmin, double hmin, double tmax, double hmax)
{
  Istate z;
  Istate_Deriv zd;
  int iter=0;
  double zdt;
  double h;

  z = *z0;
  do {
    equations_of_motion(&z, &zd);
    zdt = -1.0*pv0(z.x[0])/zd.x[0];
    if(VERBOSE) printf("refine: %.16le %.16le %.16le %.16le %.16le\n", tmin, z.t + zdt, tmax, hmin, hmax);
    if((z.t + zdt > tmax) || (z.t + zdt < tmin)) {
      zdt = (tmax+tmin)/2.0 - z.t;
    }
    if(integrate_until_success(&z, &zdt) == FAILURE) {
      exit(2);
    };
    h = pv0(z.x[0]);
    if(h*hmin > 0.0) {
      tmin = z.t;
      hmin = h;
    } else {
      tmax = z.t;
      hmax = h;
    }
    if(fabs(pv0(z.x[0])) < epssection) {
      if(VERBOSE) printf("refine success: %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n", z.t, z.x[0], z.x[1], z.x[2], z.x[3], h, zd.x[0]);
      *z0 = z;
      return(SUCCESS);
    }
  } while(iter++ < 4000);
  printf("refine iter too big %d\n", iter);
  return(FAILURE);
}

void compute_inclination_obliquity(Istate *s, double *inc, double *Ie)
{
  double h, yp, HL0, xp;
  double H0L0, H1L0, GL1, G1L1, efactor, e2;
  double e2cos2w, cosI, cosi, sinI, sini, coseps, E;

  h = s->x[0];
  yp = s->x[1];
  HL0 = s->x[2];
  xp = s->x[3];

  H0L0 = H0H12L0 + HL0;
  H1L0 = H0H12L0 - HL0;
  GL1 = (xp*xp + yp*yp)/2.0;
  G1L1 = 1.0 - GL1;
  cosI = H0L0/G0L0;
  cosi = H1L0/(G1L1*L1L0);
  sinI = sqrt(1.0 - cosI*cosI);
  sini = sqrt(1.0 - cosi*cosi);
  *inc = (180./M_PI)*asin(sini);
  *Ie = (180./M_PI)*asin(sinI);
}

int integrate_to_crossing(Istate *s0, double *dt0, double *tmin, double *hmin, double *tmax, double *hmax)
{
  Istate s;
  int iter;
  double dt;
  double h, h0;
  int flag;
  double tn, tx;
  double inc_during, Ie_during;

  dt = *dt0;
  s = *s0;

  /* do one step to make sure h0 is nonzero */
  if(VERBOSE) printf("crossing integrate: %.16le %.16le %.16le %.16le %.16le %.16le\n", s.t, pv0(s.x[0]), s.x[1], s.x[2], s.x[3], dt); 
  flag = integrate_until_success(&s, &dt);

  iter = 0;
  do {
    h0 = pv0(s.x[0]);
    if(h0 < h_min) {
      h_min = h0;
    } else if(h0 > h_max) {
      h_max = h0;
    }
    compute_inclination_obliquity(&s, &inc_during, &Ie_during);
    if(inc_during > i_max) {
      i_max = inc_during; 
    } 
    if(inc_during < i_min) {
      i_min = inc_during; 
    }
    if(Ie_during > Ie_max) {
      Ie_max = Ie_during; 
    } 
    if(Ie_during < Ie_min) {
      Ie_min = Ie_during;
    }

    if(VERBOSE) printf("crossing integrate: %.16le %.16le %.16le %.16le %.16le %.16le\n", s.t, pv0(s.x[0]), s.x[1], s.x[2], s.x[3], dt); 
    tn = s.t;
    flag = integrate_until_success(&s, &dt);
    tx = s.t;
    h = pv0(s.x[0]);
    /* printf("h %.16le %.16le\n", h0, h); */
    if((-2.0 < h*h0) && (h*h0 < 0.0)) {
      *s0 = s;
      *dt0 = dt;
      *tmin = tn;
      *tmax = tx;
      *hmin = h0;
      *hmax = h;
      return(SUCCESS);
    }
  } while(iter++ < 10000);
  if(VERBOSE) printf("crossing iter too big %d\n", iter);
  return(FAILURE);
}

int HEMmap(double xp, double yp, double dt, int npoints, double *nxp, double *nyp)
{
  Istate s0, s, z;
  Istate_Deriv sd, zd;
  double HL0;
  double direction, direction0;
  int ipoint;
  double tmin, tmax, hmin, hmax;

  HL0 = find_root(xp, yp, 5000.);
  if(VERBOSE) printf("%.16le\n", HL0);

  s0.t = 0.0;
  s0.x[0] = 0.0;
  s0.x[1] = yp;
  s0.x[2] = HL0;
  s0.x[3] = xp;
  s = s0;
  equations_of_motion(&s, &sd);
  direction0 = sgn(sd.x[0]);
  /* printf("direction0 %.15le\n", direction0); */

  for(ipoint=0; ipoint<npoints; ipoint++) {
    /* find a crossing */
    if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
      fprintf(fileerror, "-3\n");
      exit(-3);
    }
    z = s;
    if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
      exit(4);
    }
    equations_of_motion(&z, &zd);
    direction = sgn(zd.x[0]);
    if(direction > 0.0) {
      fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
      fflush(fileplus);
    } else {
      fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
      fflush(fileminus);
    }

    if(direction*direction0 > 0.0) {
      *nxp = z.x[3];
      *nyp = z.x[1];
      return(SUCCESS);
    }

    /* look for a second crossing */
    if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
      exit(5);
    }
    z = s;
    if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
      exit(6);
    }
    equations_of_motion(&z, &zd);
    direction = sgn(zd.x[0]);
    if(direction > 0.0) {
      fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
      fflush(fileplus);
    } else {
      fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
      fflush(fileminus);
    }

    if(direction*direction0 > 0.0) {
      *nxp = z.x[3];
      *nyp = z.x[1];
      return(SUCCESS);
    } else {
      if(VERBOSE) printf("failure %.16le %.16le %.16le %.16le %lf %lf\n", z.x[1], z.x[2], z.x[3], z.x[4], direction0, direction);
      HEMmap(xp, yp, dt/2.0, npoints, nxp, nyp);
    }
  }
}

/* hmin and hmax are variables used to zero in on a section crossing */
/* h_min and h_max are the minimum and maximum for a cycle between section crossings */

int HEMexplore(double xp, double yp, double dt)
{
  Istate s0, s, z;
  Istate_Deriv sd, zd;
  double HL0;
  double direction, direction0;
  int ipoint;
  double tmin, tmax, hmin, hmax;
  double inc_during, Ie_during;

  int npoints=1;

  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  HL0 = find_root(xp, yp, 5000.);
  if(VERBOSE) printf("%.16le\n", HL0);

  s0.t = 0.0;
  s0.x[0] = 0.0;
  s0.x[1] = yp;
  s0.x[2] = HL0;
  s0.x[3] = xp;
  s = s0;
  equations_of_motion(&s, &sd);
  direction0 = sgn(sd.x[0]);
  if(VERBOSE) printf("direction0 %.15le\n", direction0);

  compute_inclination_obliquity(&s, &inc_during, &Ie_during);
  fprintf(fileoutput, "%lf %.16le %.16le %.16le\n", direction0, HL0, inc_during, Ie_during);

  /* HEMmap loop over npoints begins here */
  /* find a crossing */ 
  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    fprintf(fileerror, "-3\n");
    exit(-3);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(4);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);
  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
  if(inc_during > i_max) {
    i_max = inc_during; 
    i_max_on_section = 1;
  } 
  if(inc_during < i_min) {
    i_min = inc_during; 
    i_min_on_section = 1;
  }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }

  fprintf(fileoutput, "1 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le %lf\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max, direction0);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  }

  /* look for a second crossing */
  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    exit(5);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(6);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);

  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
    if(inc_during > i_max) {
      i_max = inc_during; 
      i_max_on_section = 1;
    } 
  if(inc_during < i_min) {
      i_min = inc_during; 
      i_min_on_section = 1;
    }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }
  fprintf(fileoutput, "2 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le %lf\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max, direction0);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  } else {
    if(VERBOSE) printf("failure %.16le %.16le %.16le %.16le %lf %lf\n", z.x[1], z.x[2], z.x[3], z.x[4], direction0, direction);
    HEMexplore(xp, yp, dt/2.0);
  }
}

int HEMexplore0(double dt)
{
  double xp, yp; 
  Istate s0, s, z;
  Istate_Deriv sd, zd;
  double HL0;
  double direction, direction0;
  int ipoint;
  double tmin, tmax, hmin, hmax;
  double inc_during, Ie_during;

  int npoints=1;

  xp = 0.0;
  yp = 0.0;

  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  /* HL0 = find_root(xp, yp, 5000.); */
  HL0 = HL0e0;
  if(VERBOSE) printf("%.16le\n", HL0);

  s0.t = 0.0;
  s0.x[0] = 0.0;
  s0.x[1] = yp;
  s0.x[2] = HL0e0;
  s0.x[3] = xp;
  s = s0;
  printf("%.16le\n", energy(&s0));
  equations_of_motion(&s, &sd);
  direction0 = sgn(sd.x[0]);
  if(VERBOSE) printf("direction0 %.15le\n", direction0);

  compute_inclination_obliquity(&s, &inc_during, &Ie_during);
  fprintf(fileoutput, "%lf %.16le %.16le %.16le\n", direction0, HL0, inc_during, Ie_during);

  /* HEMmap loop over npoints begins here */
  /* find a crossing */ 
  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    fprintf(fileerror, "-3\n");
    exit(-3);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(4);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);
  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
  if(inc_during > i_max) {
    i_max = inc_during; 
    i_max_on_section = 1;
  } 
  if(inc_during < i_min) {
    i_min = inc_during; 
    i_min_on_section = 1;
  }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }

  printf("%.16le\n", energy(&z));
  fprintf(fileoutput, "1 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  }

  /* look for a second crossing */
  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    exit(5);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(6);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);

  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
    if(inc_during > i_max) {
      i_max = inc_during; 
      i_max_on_section = 1;
    } 
  if(inc_during < i_min) {
      i_min = inc_during; 
      i_min_on_section = 1;
    }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }
  printf("%.16le\n", energy(&z));
  fprintf(fileoutput, "2 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  } else {
    if(VERBOSE) printf("failure %.16le %.16le %.16le %.16le %lf %lf\n", z.x[1], z.x[2], z.x[3], z.x[4], direction0, direction);
    HEMexplore(xp, yp, dt/2.0);
  }
}

int HEMsineps(double xp, double yp, double dt)
{
  Istate s0, s, z;
  Istate_Deriv sd, zd;
  double HL0;
  double direction, direction0;
  int ipoint;
  double tmin, tmax, hmin, hmax;
  double inc_during, Ie_during;

  int npoints=1;

  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  HL0 = find_root(xp, yp, 5000.);
  if(VERBOSE) printf("%.16le\n", HL0);

  s0.t = 0.0;
  s0.x[0] = 0.0;
  s0.x[1] = yp;
  s0.x[2] = HL0;
  s0.x[3] = xp;
  s = s0;
  equations_of_motion(&s, &sd);
  direction0 = sgn(sd.x[0]);
  if(VERBOSE) printf("direction0 %.15le\n", direction0);

  compute_inclination_obliquity(&s, &inc_during, &Ie_during);
  fprintf(fileoutput, "%lf %.16le %.16le %.16le\n", direction0, HL0, inc_during, Ie_during);

  /* HEMmap loop over npoints begins here */
  /* find a crossing */ 
  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    fprintf(fileerror, "-3\n");
    exit(-3);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(4);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);
  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
  if(inc_during > i_max) {
    i_max = inc_during; 
    i_max_on_section = 1;
  } 
  if(inc_during < i_min) {
    i_min = inc_during; 
    i_min_on_section = 1;
  }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }

  fprintf(fileoutput, "1 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  }

  /* look for a second crossing */
  h_min = M_PI/2.0;
  h_max = -M_PI/2.0;
  i_min = 180.0;
  i_max = 0.0;
  Ie_min = 180.0;
  Ie_max = 0.0;
  i_min_on_section = 0;
  i_max_on_section = 0;
  Ie_min_on_section = 0;
  Ie_max_on_section = 0;

  if(integrate_to_crossing(&s, &dt, &tmin, &hmin, &tmax, &hmax) == FAILURE) {
    exit(5);
  }
  z = s;
  if(refine_crossing(&z, tmin, hmin, tmax, hmax) == FAILURE) {
    exit(6);
  }
  equations_of_motion(&z, &zd);
  direction = sgn(zd.x[0]);

  if(direction > 0.0) {
    fprintf(fileplus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileplus);
  } else {
    fprintf(fileminus, "%.16le %.16le %lf\n", z.x[3], z.x[1], direction);
    fflush(fileminus);
  }

  compute_inclination_obliquity(&z, &inc_during, &Ie_during);
    if(inc_during > i_max) {
      i_max = inc_during; 
      i_max_on_section = 1;
    } 
  if(inc_during < i_min) {
      i_min = inc_during; 
      i_min_on_section = 1;
    }
  if(Ie_during > Ie_max) {
    Ie_max = Ie_during; 
    Ie_max_on_section = 1;
  } 
  if(Ie_during < Ie_min) {
    Ie_min = Ie_during;
    Ie_min_on_section = 1;
  }
  fprintf(fileoutput, "2 %lf %.10le %.10le %.10le %.10le %d %d %.10le %.10le %.10le %d %d %.10le %.10le %.10le\n", 
	  direction,
	  z.x[2],
	  z.t,
	  h_min, h_max, 
	  i_min_on_section, i_max_on_section, inc_during, i_min, i_max, 
	  Ie_min_on_section, Ie_max_on_section, Ie_during, Ie_min, Ie_max);

  if(direction*direction0 > 0.0) {
    return(SUCCESS);
  } else {
    if(VERBOSE) printf("failure %.16le %.16le %.16le %.16le %lf %lf\n", z.x[1], z.x[2], z.x[3], z.x[4], direction0, direction);
    HEMexplore(xp, yp, dt/2.0);
  }
}
