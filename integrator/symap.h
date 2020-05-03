typedef struct{
  double t;
  VEC x[NJ], v[NJ]; //Jacobi vectors;
  SQM Oe, Om; //orien.;
  VEC Lpe, Lpm;
}SS; //sys. state;

typedef struct{
  double Ae, Ce;
  double mptE[3], eptM[3];
  double dt;
}AUX;

void H1_c(SS*, AUX*, double);
void j2b(const VEC*, VEC*);
void b2j(const VEC*, VEC*);
void ab_grav(VEC*, VEC*);
void aj_compen(VEC*, SS*, double);
void espin_orbs(const VEC*, VEC*, SQM, VEC*, double, AUX*);
void mspin_orb(SS*,double);
void tides(SS*, AUX*, double);

void map_array(SS* sp, AUX* auxp, int n)
{
  int i;
  for (i=1; i<=n; i++) {
    /* H1 */
    H1_c(sp, auxp, 1.0);
    e_flatten(sp->Lpe.R[2]*sp->Lpe.R[2]/auxp->Ce/auxp->Ce, &auxp->Ae, &auxp->Ce);
    if(WOBBLE_DAMP!=0){
      wobble_damp_Moon(&sp->Om, &sp->Lpm, (double)SPEEDUP*auxp->dt);
    }
    tides(sp, auxp, (double)SPEEDUP);
    /* H0 */
    k_evolve(sp->x, sp->v, auxp->dt);
    free_Moon(&sp->Om, &sp->Lpm, auxp->dt);
    axisymmetric(auxp->Ae, auxp->Ce, &sp->Oe, &sp->Lpe, auxp->dt);
  }

  sp->t += auxp->dt *n *(double)SPEEDUP;
}


void H1_c(SS* sp, AUX* auxp, double step_ratio)
{
  VEC xb[NJ+1], ab[NJ+1], aj[NJ];

  j2b(sp->x, xb);
  ab_grav(xb, ab);
  espin_orbs(xb, ab, sp->Oe, &sp->Lpe, step_ratio*auxp->dt, auxp);
  b2j(ab, aj);
  aj_compen(aj, sp, step_ratio*auxp->dt);
  mspin_orb(sp, step_ratio*auxp->dt);
}


void map_period(double a, AUX* auxp)
{
  double a_ref=60.0*r_e, T_ref=27.3;
  auxp->dt = T_ref *pow(a/a_ref, 1.5) /(double)MAP_FREQ;
}


void tides(SS* sp, AUX* auxp, double step_ratio)
{
    /* public part for both tides */
    double raMf[4], e;
    xv2elem_invar(kconst[0], sp->x+0, sp->v+0, raMf+0, raMf+1, &e, raMf+2, raMf+3);
    double G[TIDE_L+1][TDOMAX+1];
    int order=TDOMAX;
    if(e<0.0015)
      order=2;
    else if(e<0.04)
      order=4;
    else if(e<0.115)
      order=6;
    else if(e<0.2)
      order=8; //error control: 3e-6;
    update_hansen(e, G, order);
    VEC nhat = vhat( cr_p(sp->x[0], sp->v[0]) );
    
    /* Earth's deformation's effect */
    VEC ab_m, aj0_0;

    /* Moon's deformation's effect */
    VEC ab_e, aj0_1;

      {
	VEC w_, what;
	double wlen;
	v_cast(&w_, sp->Lpe.R[0]/auxp->Ae, sp->Lpe.R[1]/auxp->Ae, sp->Lpe.R[2]/auxp->Ce);
	w_ = Mv(sp->Oe, w_);
	wlen = vlen(w_);
	what = cv(1./wlen, w_);
	ab_m = tidal_accl(E_DEFOR, what, wlen, sp->x[0], nhat, G, order, auxp->mptE,auxp->dt, raMf); // change: Jul 26
	aj0_0 = cv(mlis[1]/mr, ab_m);
      }

      {
	VEC w_, what;
	double wlen;
	v_cast(&w_, sp->Lpm.R[0]/Am, sp->Lpm.R[1]/Bm, sp->Lpm.R[2]/Cm);
	w_ = Mv(sp->Om, w_);
	wlen = vlen(w_);
	what = cv(1./wlen, w_);
	ab_e = tidal_accl(M_DEFOR, what, wlen, cv(-1.0, sp->x[0]), nhat, G, order, auxp->eptM,auxp->dt, raMf); // change: Jul 26
	aj0_1 = cv(-mlis[0]/mr, ab_e);
      }

    /* total effect on first Jacobi vector, and on Earth's rotation */
    double dt = step_ratio*auxp->dt;
    VEC aj0 = vsum(aj0_0, aj0_1);
    vgrow(sp->v+0, cv(dt, aj0));

    VEC Dt_Lpm = trMv(sp->Om, cv(mr, cr_p(aj0_1, sp->x[0])));
    vgrow(&sp->Lpm, cv(dt, Dt_Lpm));

    VEC Dt_Lpe = trMv(sp->Oe, cv(mr, cr_p(aj0_0, sp->x[0])));
    vgrow(&sp->Lpe, cv(dt, Dt_Lpe));
}



void espin_orbs(const VEC* xb_p, VEC* ab_p, SQM Orien, VEC* Lp_p, double dt, AUX* auxp)
{
  VEC xp, Dpp, Dp;
  int i;
  for(i=1; i<NJ+1; i++){
    xp = trMv(Orien, vdif(xb_p[i], xb_p[0]));
    spin_orbit_kick(auxp->Ae, auxp->Ae, auxp->Ce, mlis[i], xp, dt, Lp_p, &Dpp);
    Dp = Mv(Orien, Dpp);
    vgrow(ab_p+i, cv(1./mlis[i], Dp));
    vgrow(ab_p+0, cv(-1./mlis[0], Dp));
  }
}
void mspin_orb(SS* sp, double dt)
{
  VEC xp = trMv(sp->Om, cv(-1., sp->x[0])), Dpp;
  spin_orbit_kick(Am, Bm, Cm, mlis[0], xp, dt, &sp->Lpm, &Dpp);
  vgrow(sp->v+0, cv(-dt/mr, Mv(sp->Om, Dpp)));
}


void ab_grav(VEC* xb_p, VEC* ab_p)
{
  int i,j;
    
  v_cast(ab_p+0, 0., 0., 0.);
  v_cast(ab_p+1, 0., 0., 0.);
  v_cast(ab_p+2, 0., 0., 0.);
    
  i=0;   // ------ omitting the "0"-"1" ordinary vector part of N-body Hamiltonian;
  j=2;
  {
    VEC xis=vdif(xb_p[j], xb_p[i]);
    double r2=in_p(xis,xis);
    double c=-G_const/r2/sqrt(r2);
    vgrow(ab_p+i, cv(mlis[j]*c, vdif(xb_p[i], xb_p[j])));
    vgrow(ab_p+j, cv(mlis[i]*c, vdif(xb_p[j], xb_p[i])));
  }

  i=1;
  j=2;
  {
    VEC xis=vdif(xb_p[j], xb_p[i]);
    double r2=in_p(xis,xis);
    double c=-G_const/r2/sqrt(r2);
    vgrow(ab_p+i, cv(mlis[j]*c, vdif(xb_p[i], xb_p[j])));
    vgrow(ab_p+j, cv(mlis[i]*c, vdif(xb_p[j], xb_p[i])));
  }
}


void aj_compen(VEC* aj_p, SS* sp, double dt) //no shape/rotation;
{
    vgrow(sp->v+0, cv(dt, aj_p[0]));
    vgrow(sp->v+1, cv(dt, aj_p[1]));
    int i=1;
    {
      double r2 = in_p(sp->x[i],sp->x[i]);
      vgrow(sp->v+i, cv( dt*kconst[i]/r2/sqrt(r2), sp->x[i]));
    }
}


void b2j(const VEC b[], VEC j[])
{
    j[0] = vdif(b[1], b[0]);
    j[1] = vdif(b[1], cv(wtI[0], j[0]));
    j[1] = vdif(b[2], j[1]);
}


void j2b(const VEC j[], VEC b[])
{
    b[2] = cv(wtI[1], j[1]);
    b[1] = vdif(b[2], j[1]);
    vgrow(b+1, cv(wtI[0], j[0]));
    b[0] = vdif(b[1], j[0]);
}


void set_angbank(SS* sp, AUX* auxp)
{
  double raMf[4], n, foo;
  xv2elem_invar(kconst[0], sp->x+0, sp->v+0, raMf+0, raMf+1, &foo, raMf+2, raMf+3);
  n = sqrt(kconst[0]/raMf[1]/raMf[1]/raMf[1]);
  VEC nhat = vhat( cr_p(sp->x[0], sp->v[0]) );

  // as tidal_accl(----, chat, -----, r_, nhat, auxp->mptE, ----);
  {
    SQM ins2xyz;
    VEC chat = sp->Oe.C[2];
    ins2xyz.C[2] = chat;  // "zhat"

    if(chat.R[2] > 1.0-1.0e-12 || chat.R[2] < -1.0+1.0e-12){
      double ang_cONxy_xINER = atan2( chat.R[1], chat.R[0] );
      v_cast(ins2xyz.C+0, cos(ang_cONxy_xINER + M_PI_2), sin(ang_cONxy_xINER + M_PI_2), 0.);
    } else{
      ins2xyz.C[0] = vhat( cr_p(zhat(), chat) );
    }

    ins2xyz.C[1] = cr_p(chat, ins2xyz.C[0]); // "yhat"

    VEC rp, rp_hat, np_hat;
    VEC r_ = sp->x[0];
    rp = trMv(ins2xyz,r_);
    rp_hat = cv(1.0/raMf[0], rp);
    np_hat = trMv(ins2xyz,nhat);
    get_incli_node_w(rp_hat, np_hat, raMf[3], &foo, auxp->mptE+NODE, auxp->mptE+PERI);
    *(auxp->mptE + MeanAno) = raMf[2] - n*auxp->dt;
  }

  // as tidal_accl(----, chat, ----, r_, nhat, auxp->eptM, ----);
  {
    SQM ins2xyz;
    VEC chat = sp->Om.C[2];
    ins2xyz.C[2] = chat;  // "zhat"

    if(chat.R[2] > 1.0-1.0e-12 || chat.R[2] < -1.0+1.0e-12){
      double ang_cONxy_xINER = atan2( chat.R[1], chat.R[0] );
      v_cast(ins2xyz.C+0, cos(ang_cONxy_xINER + M_PI_2), sin(ang_cONxy_xINER + M_PI_2), 0.);
    } else{
      ins2xyz.C[0] = vhat( cr_p(zhat(), chat) );
    }

    ins2xyz.C[1] = cr_p(chat, ins2xyz.C[0]); // "yhat"

    VEC rp, rp_hat, np_hat;
    VEC r_ = cv(-1.0, sp->x[0]);
    rp = trMv(ins2xyz,r_);
    rp_hat = cv(1.0/raMf[0], rp);
    np_hat = trMv(ins2xyz,nhat);
    get_incli_node_w(rp_hat, np_hat, raMf[3], &foo, auxp->eptM+NODE, auxp->eptM+PERI);
    *(auxp->eptM + MeanAno) = raMf[2] - n*auxp->dt;
  }
}


void map_phase_adjust(SS* sp, AUX* auxp)
{
  sp->t = 0.0;
  k_evolve(sp->x, sp->v, -0.5*auxp->dt);
  free_Moon(&sp->Om, &sp->Lpm, -0.5*auxp->dt);
  axisymmetric(auxp->Ae, auxp->Ce, &sp->Oe, &sp->Lpe, -0.5*auxp->dt);

  //e_flatten() called before ini_Earth();
  H1_c(sp, auxp, -1.);
  set_angbank(sp,auxp);
  tides(sp, auxp, -(double)SPEEDUP);

  SS s_pre_leap = *sp;
  k_evolve(sp->x, sp->v, -1.0*auxp->dt);
  //free_Moon or axisymmetric-for-Earth;
  H1_c(sp, auxp, -1.0);
  //e_flatten, wobble_damp;
  tides(sp, auxp, -(double)SPEEDUP); //angbank updated;
  *sp = s_pre_leap; //ready for maps.
}


double angmom(const SS* sp)
{
  VEC AM = cv(mr, cr_p(sp->x[0],sp->v[0]));
  vgrow(&AM, Mv(sp->Oe, sp->Lpe));
  vgrow(&AM, Mv(sp->Om, sp->Lpm));
  return vlen(AM);
}

