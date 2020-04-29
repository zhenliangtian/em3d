double dampbase_m;

#define radius_e 6371.0/1.4959787e8 //AU/km;
#define radius_m 1737.0/1.4959787e8
const double r_e = radius_e;
const double r_m = radius_m;
const double Ce_tp = 0.3308 *mass_e*radius_e*radius_e;
const double J2MR2e_tp = 0.00108263 *mass_e*radius_e*radius_e;
const double w2e_tp = (M_PI+M_PI)*(M_PI+M_PI);
const double Am = 0.385559714218574 *mass_m*radius_m*radius_m;
const double Bm = 0.3964894783411223 *mass_m*radius_m*radius_m;
const double Cm = 0.4 *mass_m*radius_m*radius_m;


//---------Earth's shape: variable
void e_flatten(double omega_sq, double* A_p, double* C_p)
{
  double x = J2MR2e_tp*(omega_sq/w2e_tp - 1.0);
  *A_p = (Ce_tp -J2MR2e_tp) -x/3.;
  *C_p = Ce_tp +2.*x/3.;
}
//-------------------------


//---------------Moon and Earth's initial rotational state (Lie-Poisson approach)
void ini_Moon(double n, double I_ref_eclipt, VEC* Lp_p, SQM* O_p)
{
  v_cast(Lp_p, 0., 0., Cm*n); // Cm*n*1.01;
  *O_p = R_a(I_ref_eclipt);
}

void ini_Earth(double daylen, double Ce, double I, VEC* Lp_p, SQM* O_p)
{
  v_cast(Lp_p, 0., 0., Ce*TWOPI/daylen);
  *O_p = R_a(I);
}
//------------------------------------------------------------



//-------------------internal dissipation in Moon and Earth's rotational kinetic energy
void wobble_damp_Moon(SQM* O_p, VEC* Lp_p, double dt)
{
  VEC omegap;
  v_cast(&omegap, Lp_p->R[0]/Am, Lp_p->R[1]/Bm, Lp_p->R[2]/Cm);
  double spin = vlen(omegap);
  double damprate = dampbase_m *spin*spin*spin *WOBBLE_DAMP;
  double fac = exp(-damprate*dt);

  VEC Lp0=*Lp_p, Lp1;
  double Lplen = vlen(Lp0);
  v_cast(&Lp1, Lp0.R[0]*fac, Lp0.R[1]*fac, Lp0.R[2]);
  Lp1 = cv( Lplen/vlen(Lp1), Lp1 );
  *Lp_p = Lp1;

  VEC kp = cr_p(Lp0,Lp1); //axis of E-R rotation;
  double s = vlen(kp)/Lplen/Lplen;
  if(s!=0.){
    double c = sqrt(1.-s*s);
    *O_p = MN(*O_p, ezR_v(kp, c, -s));
  }
}
//------------------------------------------------------------



void spin_orbit_kick(double A, double B, double C, double m, VEC xp, double dt, VEC* Lp_p, VEC* Dpp_p)
{
  double r2 = in_p(xp, xp);
  double foo = 3.0*G_const*m/r2/r2/sqrt(r2);
  VEC Irp;
  v_cast(&Irp, A*xp.R[0], B*xp.R[1], C*xp.R[2]);
  vgrow(Lp_p, cv(dt*foo, cr_p(xp, Irp)));

  *Dpp_p = cv(-0.5*foo *(A+B+C), xp);
  vgrow(Dpp_p, cv(-foo, Irp));
  vgrow(Dpp_p, cv(2.5*foo/r2 *in_p(xp, Irp), xp));
}



//----------------free Moon and Earth rotations-------------
void axisymmetric(double B, double C, SQM* O_p, VEC* Lp_p, double dt)
{
    VEC Lp0 = *Lp_p;

    double alpha = Lp0.R[2]/B -Lp0.R[2]/C;
    SQM trR = R_c(-alpha*dt);
    *Lp_p = trMv(trR, Lp0);

    double beta = -vlen(Lp0)/B;
    *O_p = MN(MN(*O_p, R_v(Lp0,-beta*dt)), trR);
}

void triaxial_Moon(SQM* O_p, VEC* Lp_p, double dt)
{
    double gamma = Lp_p->R[0]/Bm -Lp_p->R[0]/Am;
    SQM trR = R_a(-gamma*dt);

    *Lp_p = trMv(trR, *Lp_p);
    *O_p = MN(*O_p, trR);
}

void free_Moon(SQM* O_p, VEC* Lp_p, double dt)
{
  axisymmetric(Bm, Cm, O_p, Lp_p, 0.5*dt);
  triaxial_Moon(O_p, Lp_p, dt);
  axisymmetric(Bm, Cm, O_p, Lp_p, 0.5*dt);
}
//------------------------------------------------------------

