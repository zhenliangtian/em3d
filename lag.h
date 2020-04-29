#define TIDE_L 2

#define E_DEFOR 0
#define M_DEFOR 1
double cos_oQ[2], sin_oQ[2], c_k2[2];

#define PERI 0
#define NODE 1
#define MeanAno 2

void update_kaula(double, double[TIDE_L+1][TIDE_L+1]);
double P__2m(int, double);
double DP_2m(int, double);
double factorial(int);
double kronecker(int, int);
int even_judge(int);  //returns 1 on even;
double pv_pi(double); //to [-pi, pi)

VEC tidal_accl(unsigned char id, VEC chat, double w_c, VEC r_, VEC nhat, double G[TIDE_L+1][TDOMAX+1], int order, double* angbank, double diff_dt, double* raMf)
{
  SQM ins2xyz; //ins2xyz: vector's instantaneous-inertia-frame-components -> space-xyz-components;
  ins2xyz.C[2] = chat;  //instantaneous-inertia-coor-sys's "zhat";

  if(chat.R[2] > 1.0-1.0e-12 || chat.R[2] < -1.0+1.0e-12){
    double ang_cONxy_xINER = atan2( chat.R[1], chat.R[0] );
    v_cast(ins2xyz.C+0, cos(ang_cONxy_xINER + M_PI_2), sin(ang_cONxy_xINER + M_PI_2), 0.);
  } else{
    ins2xyz.C[0] = vhat( cr_p(zhat(), chat) );
  }

  ins2xyz.C[1] = cr_p(chat, ins2xyz.C[0]); //instantaneous-inertia-coor-sys's "yhat";

  VEC rp, rp_hat, np_hat;
  double rlen = raMf[0];
  rp = trMv(ins2xyz, r_);
  rp_hat = cv(1.0/rlen, rp);
  np_hat = trMv(ins2xyz, nhat);
  double M=raMf[2], peri, node,incli;
  get_incli_node_w(rp_hat, np_hat, raMf[3], &incli, &node, &peri);

  double peri_dot, M_dot, node_dot;
  peri_dot = pv_pi(peri-angbank[PERI])/diff_dt;
  node_dot = pv_pi(node-angbank[NODE])/diff_dt;
  M_dot = pv_pi(M-angbank[MeanAno])/diff_dt;

  *(angbank + PERI) = peri;
  *(angbank + NODE) = node;
  *(angbank + MeanAno) = M;

  double kaulaF[TIDE_L+1][TIDE_L+1];
  update_kaula(incli, kaulaF);
  
  double theta = 0.0;
  double lam = atan2_pvTWOPI(rp.R[1], rp.R[0]);
  int m,p,q;
  double T_levelm[TIDE_L+1], dTdlam_levelm[TIDE_L+1];
  double f_ie, combi_A, term_A, term_B;

  double cos_incre_levq, sin_incre_levq;
  {
    cos_incre_levq = cos(M);
    if(cos_incre_levq > 1.0 -1.0e-3  ||  cos_incre_levq < -1.0 +1.0e-3)
      sin_incre_levq = sin(M);
    else
      if (pv_pi(M)>0)
	sin_incre_levq = sqrt(1.0 - cos_incre_levq*cos_incre_levq);
      else
	sin_incre_levq = -1.0 *sqrt(1.0 - cos_incre_levq*cos_incre_levq);
  }
    
  double cos_levq_help, sin_levq_help;
  {
    double angle_levq_help = (-order)*M -M;
    cos_levq_help = cos(angle_levq_help);
    if(cos_levq_help > 1.0 -1.0e-3  ||  cos_levq_help < -1.0 +1.0e-3)
      sin_levq_help = sin(angle_levq_help);
    else
      if (pv_pi(angle_levq_help)>0)
    	sin_levq_help = sqrt(1.0 - cos_levq_help*cos_levq_help);
      else
    	sin_levq_help = -1.0 *sqrt(1.0 - cos_levq_help*cos_levq_help);
  }

  double cos_incre_levP, sin_incre_levP;
  {
    double angle_incre_levP = -2.0*(peri+M);
    cos_incre_levP = cos( angle_incre_levP );
    if(cos_incre_levP > 1.0 -1.0e-3  ||  cos_incre_levP < -1.0 +1.0e-3)
      sin_incre_levP = sin(angle_incre_levP);
    else
      if (pv_pi(angle_incre_levP)>0)
	sin_incre_levP = sqrt(1.0 - cos_incre_levP*cos_incre_levP);
      else
	sin_incre_levP = -1.0 *sqrt(1.0 - cos_incre_levP*cos_incre_levP);
  }


  double cos_incre_LEVM, sin_incre_LEVM;
  {
    double angle_incre_LEVM = (node - lam - theta);
    cos_incre_LEVM = cos( angle_incre_LEVM );
    if(cos_incre_LEVM > 1.0 -1.0e-3  ||  cos_incre_LEVM < -1.0 +1.0e-3)
      sin_incre_LEVM = sin(angle_incre_LEVM);
    else
      if (pv_pi(angle_incre_LEVM)>0)
	sin_incre_LEVM = sqrt(1.0 - cos_incre_LEVM*cos_incre_LEVM);
      else
	sin_incre_LEVM = -1.0 *sqrt(1.0 - cos_incre_LEVM*cos_incre_LEVM);
  }


  double cos_prev_LEVM, sin_prev_LEVM;
  {
    double angle_prev_LEVM = (TIDE_L-2*0)*(peri+M)  +2*(peri+M) + -1.0 *(node - lam - theta);
    cos_prev_LEVM = cos( angle_prev_LEVM );
    if(cos_prev_LEVM > 1.0 -1.0e-3  ||  cos_prev_LEVM < -1.0 +1.0e-3)
      sin_prev_LEVM = sin(angle_prev_LEVM);
    else
      if (pv_pi(angle_prev_LEVM)>0)
	sin_prev_LEVM = sqrt(1.0 - cos_prev_LEVM*cos_prev_LEVM);
      else
	sin_prev_LEVM = -1.0 *sqrt(1.0 - cos_prev_LEVM*cos_prev_LEVM);
  }



  for (m=0; m<=TIDE_L; m++) {
    T_levelm[m] = 0.0;
    dTdlam_levelm[m] = 0.0;


    double cos_prev_levP, sin_prev_levP;
    cos_prev_levP = cos_prev_LEVM *cos_incre_LEVM - sin_prev_LEVM *sin_incre_LEVM;
    sin_prev_levP = sin_prev_LEVM *cos_incre_LEVM + cos_prev_LEVM *sin_incre_LEVM;
    cos_prev_LEVM = cos_prev_levP;
    sin_prev_LEVM = sin_prev_levP;


    for (p=0; p<=TIDE_L; p++) {
      combi_A = (TIDE_L-2*p)*(peri_dot+M_dot) + m*(node_dot - w_c);
      // combi_B is: (TIDE_L-2*p)*(peri+M) + m* (node - lam - theta);


      double cos_combi_B, sin_combi_B;
      cos_combi_B = cos_prev_levP *cos_incre_levP - sin_prev_levP *sin_incre_levP;
      sin_combi_B = sin_prev_levP *cos_incre_levP + cos_prev_levP *sin_incre_levP;
      cos_prev_levP = cos_combi_B;
      sin_prev_levP = sin_combi_B;

      double cos_prev_levq = cos_combi_B *cos_levq_help - sin_combi_B *sin_levq_help;
      double sin_prev_levq = sin_combi_B *cos_levq_help + cos_combi_B *sin_levq_help;
      // angle_levq_help means: (-order)*M -M;


      for (q=-order; q<=order; q++) {
	f_ie = kaulaF[m][p] *(q>=0 ? G[p][q] : G[TIDE_L-p][-q]);

	double virtual_cos = cos_prev_levq *cos_incre_levq - sin_prev_levq *sin_incre_levq;
	double virtual_sin = sin_prev_levq *cos_incre_levq + cos_prev_levq *sin_incre_levq;
	cos_prev_levq = virtual_cos;
	sin_prev_levq = virtual_sin;

	double foruse_cos, foruse_sin;
	double checksign = combi_A + q*M_dot;

	if (checksign>0.0){
	  foruse_cos = virtual_cos *cos_oQ[id] + virtual_sin *sin_oQ[id];
	  foruse_sin = virtual_sin *cos_oQ[id] - virtual_cos *sin_oQ[id];
	}
	else if (checksign<0.0){
	  foruse_cos = virtual_cos *cos_oQ[id] - virtual_sin *sin_oQ[id];
	  foruse_sin = virtual_sin *cos_oQ[id] + virtual_cos *sin_oQ[id];
	}else{ //checksign is "0.0";
	  foruse_cos = virtual_cos;
	  foruse_sin = virtual_sin;
	}



	if (even_judge(TIDE_L-m)) {
	  term_A = f_ie *foruse_cos;
	  term_B = f_ie*(-1.0) *foruse_sin *(-1.0)*m;
	}
	else {
	  term_A = f_ie *foruse_sin;
	  term_B = f_ie *foruse_cos *(-1.0)*m;
	}
	T_levelm[m] += term_A;
	dTdlam_levelm[m] += term_B;
      }
    }
  }


  double r_xy, s_phi, c_phi, s_lam, c_lam;
  r_xy = sqrt(rp.R[0]*rp.R[0] + rp.R[1]*rp.R[1]);
  s_phi = rp_hat.R[2];
  c_phi = r_xy/rlen;
  s_lam = rp.R[1]/r_xy;
  c_lam = rp.R[0]/r_xy;
  double term_1, term_P;  // used in the "summation with running m index" process;
  double T_part = 0.0, dTdlam_part = 0.0, dTdphi_part = 0.0, fac0;
  for (m=0; m<=TIDE_L; m++) {
    term_1 = factorial(TIDE_L-m)/factorial(TIDE_L+m)*(2.0 - kronecker(0, m));
    term_P = P__2m(m, s_phi);
    T_part += term_1 *term_P* T_levelm[m];
    dTdlam_part += term_1 *term_P *dTdlam_levelm[m];
    dTdphi_part += term_1 *DP_2m(m, s_phi) *c_phi *T_levelm[m];
  }
  fac0 = c_k2[id] * pow(raMf[1]*rlen, -TIDE_L-1);

  double dTdlam, dTdphi, dTdr;
  dTdlam = fac0 * dTdlam_part;
  dTdphi = fac0 *dTdphi_part;
  dTdr = -(TIDE_L+1)/rlen *fac0 *T_part;
  
  VEC rhat, phihat, lamhat, phihatp, lamhatp, a_ast;
  // "phihatp" is the "phihat" components in equator. coor.
  v_cast(&phihatp, -s_phi*c_lam, -s_phi*s_lam, c_phi);
  v_cast(&lamhatp, -s_lam, c_lam, 0.);
  rhat = Mv(ins2xyz, rp_hat);
  phihat = Mv(ins2xyz, phihatp);
  lamhat = Mv(ins2xyz, lamhatp);

  a_ast = vsum( cv(1.0/rlen*dTdphi, phihat), cv(1.0/r_xy *dTdlam, lamhat) );
  a_ast = vsum( cv(dTdr, rhat), a_ast );
  
  return a_ast;
}


//----------------------functions directly related with "tidal_accl()"
void update_kaula(double i, double F[TIDE_L+1][TIDE_L+1])
{
    double sini = sin(i), cosi = cos(i);
    double sini_sq = sini*sini;
    double cosi_sq = cosi*cosi;
    double sini_cosi = sini*cosi;
    
    F[0][0] = -0.375*sini_sq ;
    F[0][1] = 0.75*sini_sq - 0.5 ;
    F[0][2] = -0.375*sini_sq ;
    F[1][0] = 0.75*(sini_cosi + sini) ;
    F[1][1] = -1.5*sini_cosi ;
    F[1][2] = 0.75*(sini_cosi - sini) ;
    F[2][0] = 0.75*(1.0 + 2.0*cosi + cosi_sq) ;
    F[2][1] = 1.5*sini_sq ;
    F[2][2] = 0.75*(1.0 - 2.0*cosi + cosi_sq) ;
}

double P__2m(int m, double x)  // degree l==2;
{
  if (m==0) return(1.5*x*x - 0.5);
  else if (m==1) return( 3.0*x*sqrt(1.0 - x*x) );
  else if (m==2) return(3.0 - 3.0*x*x);
  else {
    printf("\n error in calling the 'P__2m' function.");
    exit(-1);
  }
}

double DP_2m(int m, double x)  // degree l==2;
{
  if(m==0) return(3.0*x);
  else if(m==1) return( (3.0 - 6.0*x*x)/sqrt(1.0-x*x) );
  else if(m==2) return(-6.0*x);
  else {
    printf("\n error in calling the 'DP_2m' function.");
    exit(-1);
  }
}

double factorial(int n) // n in limited region: 0~4;
{
    if(n==0) return(1);
    else if(n==1) return(1);
    else if(n==2) return(2);
    else if(n==3) return(6);
    else if(n==4) return(24);
    else {
        printf("\n error in input number for computing its factorial\n");
	exit(-1);
    }
}

double kronecker(int a, int b)
{
    if(a==b) return(1.0);
    else return(0.0);
}

int even_judge(int n)   // returns 1 if n is even;
{
    if(n%2 == 0) return(1);
    else return(0);
}

double pv_pi(double x)
{
  // to [-pi,pi);
  double pv = x;
  while (pv >= M_PI)
    pv -= TWOPI;
  while (pv < -M_PI)
    pv += TWOPI;
  return pv;
}
//--------------------------------------------------------------------------------



//----------------section for update_hansen-------
void update_hansen(double e, double G[TIDE_L+1][TDOMAX+1], int order)
{
  static double newc_[TIDE_L+1][(TDOMAX/2+1)*((TDOMAX+1)/2+1)]={{0}};
  //            newcomb^{(-L-1),(L-2p)}_{q+i,i};
  void ini_newcomb(double[TIDE_L+1][(TDOMAX/2+1)*((TDOMAX+1)/2+1)]);
  static char init = 'n';
  if(init=='n'){
    ini_newcomb(newc_);
    init = 'y';
  }


  short p,q,i;
  int index;
  double base, sum, e2=e*e, e2s[TDOMAX/2+1]={1};
  for(i=1; i<=order/2; i++)
    e2s[i] = e2*e2s[i-1];

  for (p=0; p<=TIDE_L; p++) {
    index = 0;
    base = 1;

    for(q=0; q<=order; q++) {
      sum = 0;
      for(i=0; q+2*i<=order; i++){
	sum += newc_[p][index] * e2s[i];
	index++;
      }
      G[p][q] = base*sum;
      base *= e;
    }
  }
}

void ini_newcomb(double a[TIDE_L+1][(TDOMAX/2+1)*((TDOMAX+1)/2+1)])
{
    double f_newc(int, int, int, char[][TDOMAX/2+1][TDOMAX+1], double[][TDOMAX/2+1][TDOMAX+1]);
    char flags[TDOMAX+2*TIDE_L+1][TDOMAX/2+1][TDOMAX+1] = {{{0}}};
    double vals[TDOMAX+2*TIDE_L+1][TDOMAX/2+1][TDOMAX+1] = {{{0}}};

    short p,q,i;
    for(p=0; p<=TIDE_L; p++){
      int index = 0;
      for(q=0; q<=TDOMAX; q++) {
	// q domain [0,TDOMAX] is enough since G_l,p,-|q| = G_l,l-p,|q|;
	for(i=0; q+2*i<=TDOMAX; i++) {
	  a[p][index] = f_newc(TIDE_L -2*p, q+i, i, flags, vals);
	  index++;
	}
      }
    }
}


double f_newc(int b, int c, int d, char flags[][TDOMAX/2+1][TDOMAX+1], double vals[][TDOMAX/2+1][TDOMAX+1])
{
    int a = -TIDE_L -1;
    
    double spe_combi(int);
    double val, sum, c_sign;
    int j;

    if (d<0 || c<0)
      return(0.0);
    
    if(flags[TIDE_L+b+c-d][d][c] != 1) {

      if (d==0) {
	if (c==0)
	  val = 1.0;
	else if (c==1)
	  val = (double)b - a/2.0;
	else{
	  val = ((double)b-a/2.0)/c*f_newc(b+1, c-1, 0, flags, vals);
	  val += ((double)b-a)/(4.0*c)*f_newc(b+2, c-2, 0, flags, vals);
	}
      }

      else if (c<d) // i.e.,  when d > c >= 0 , 
	val = f_newc(-b, d, c, flags, vals);

      else {  // i.e. ,   when c >= d > 0;
	val = -((double)b+a/2.0)/d*f_newc(b-1, c, d-1, flags, vals);
	val += -((double)b+a)/(4.0*d)*f_newc(b-2, c, d-2, flags, vals);
	val += -( ((double)c+4.0+4.0*b+a)/(4.0*d) - 1.25 )*f_newc(b, c-1, d-1, flags, vals);
	sum = 0.0;  c_sign = -1.0;
	for (j=2; j<=d; j++) {
	  c_sign *= (-1.0);
	  sum += c_sign*spe_combi(j)*f_newc(b, c-j, d-j, flags, vals);
	}
	val += ( ((double)c+b)/(2.0*d) - 0.5 )*sum;
      }

      vals[TIDE_L+b+c-d][d][c] = val;
      flags[TIDE_L+b+c-d][d][c] = 1;
    }

    return(vals[TIDE_L+b+c-d][d][c]);
}


double spe_combi(int j) // combi( 3.0/2.0, j);   only called when j >= 2 ;
{
    if (j==2)
      return(0.375);
    else
      return( (1.5-((double)j-1.0))/(double)j*spe_combi(j-1) );
}
//------------------------------------------------
