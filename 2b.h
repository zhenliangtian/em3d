char flag_unbound_orb = 0;

#define NJ 2
#define mass_e 3.003439889929972e-6
#define mass_m 3.6949711649751186e-8
const double mlis[NJ+1] = {mass_e, mass_m, 1.0};
const double G_const = 2.9591220828559115e-4;
const double mr = mass_e*mass_m/(mass_e+mass_m);

#define K_PRECISION 1.0e-12  // used in "M_to_E" and "getDeltaE";

typedef struct{
  double a, e, M, node, peri, incli;
}k_elem;

double kconst[NJ], wtI[NJ];

double atan2_pvTWOPI(double, double);


void set_kepler_para(void)
{
  double k_r = G_const*mlis[0];
  double eta[NJ+1];
  eta[0] = mlis[0];
  int i;
  for(i=0; i<NJ; i++){
    eta[i+1] = eta[i] + mlis[i+1];
    wtI[i] = eta[i]/eta[i+1];
    kconst[i] = k_r/wtI[i];
  }
  //----- k_i(=mu_i/mi'), or mu_i, reflects splitting of Hamiltonian;
}


void ini_nb(VEC* xp, VEC* vp, double a, double e, double I, double eps)
{
    void k_elem_to_x_v(k_elem, double, VEC*, VEC*);
    k_elem planet[NJ];
    
    //----------initial orbital elements
    planet[0].a = a;
    planet[0].e = e;
    planet[0].incli = I+eps;
    planet[0].node = 0.0;
    planet[0].peri = 0.0;
    planet[0].M = 0.0;
    
    planet[1].a = 1.0 ;
    planet[1].e = 0.0;
    planet[1].incli = 0.0 ;
    planet[1].node = 0.0;
    planet[1].peri = 0.0;
    planet[1].M = 0.0;
    //----------------------------------

    k_elem_to_x_v(planet[0], kconst[0], xp+0, vp+0);
    double eta1=mlis[0]+mlis[1];
    k_elem_to_x_v(planet[1], kconst[1]*eta1/mlis[0], xp+1, vp+1);
}


void xv2elem_invar(double k, const VEC* xp, const VEC* vp, double* rp, double* ap, double* ep, double* Mp, double* fp)
{
  double r = vlen(*xp);
  *rp = r;
  
  double den = in_p(*vp,*vp)/k + (-2.0)/r;
  if (den>=0.0) {
    printf("\nxv2elem_invar: UNBOUND ORBIT.\n");
    flag_unbound_orb = 1;
    raise(SIGINT);
  }
  double a = -1.0/den;
  *ap= a;
    
  double esinE=in_p(*xp,*vp)/(sqrt(k*a));
  double ecosE=1.0-r/a;

  double e2 = ecosE*ecosE+esinE*esinE;
  *ep = sqrt(e2);
  
  double E = atan2_pvTWOPI(esinE, ecosE);
  *Mp= E - esinE;
  
  *fp = atan2_pvTWOPI(esinE*sqrt(1.0-e2), ecosE - e2);
}


void get_incli_node_w(VEC rhat, VEC nhat, double f, double* ip, double* nodep, double* wp)
{
  double i,node;
    
  if ((nhat.R[2] < 1.0-1.0e-6) && (nhat.R[2] > -1.0+1.0e-6))
    i=acos(nhat.R[2]);
  else {
    i = asin(sqrt(nhat.R[0]*nhat.R[0] + nhat.R[1]*nhat.R[1]));
    if (nhat.R[2] < 0.0)
      i = M_PI - i;
  }
  *ip=i;
  
  node = atan2_pvTWOPI(nhat.R[0], nhat.R[1]*(-1.0) );
  *nodep=node;
  
  VEC nodeline;
  if( nhat.R[2] < (1.0-5.0e-5) && nhat.R[2] > (-1.0+5.0e-5) ){
    double augfac = 1.0/sqrt(1.0-nhat.R[2]*nhat.R[2]);
    v_cast(&nodeline, -nhat.R[1]*augfac, nhat.R[0]*augfac, 0.);
  }
  else
    v_cast(&nodeline, cos(node), sin(node), 0.);

  double totalargu, c, s;
  c = in_p(nodeline, rhat);
  s = in_p(cr_p(nodeline, rhat), nhat);
  totalargu = atan2_pvTWOPI(s, c);
  *wp = totalargu - f;
}


//---------------------k elem to x,v vectors----------------------
void k_elem_to_x_v (k_elem elem, double k, VEC* xp, VEC* vp)
{
    double M_to_E(double,double);
    double a=elem.a,e=elem.e,M=elem.M,node=elem.node,peri=elem.peri,incli=elem.incli;
    double E;
    VEC x,v;
    
    E=M_to_E(e, M);
    double cE=cos(E), sE=sin(E), cphi=sqrt(1.0-e*e);
    
    x.R[0]=a*cE-a*e;
    x.R[1]=a*cphi*sE;
    x.R[2]=0.0;
    
    SQM R;
    R = MN(R_a(incli), R_c(peri));
    R = MN(R_c(node), R);

    *xp = Mv(R, x);

    double Edot = sqrt(k/(a*a*a)) /(1.0-e*cos(E));
    
    v.R[0]=-a*sE*Edot;
    v.R[1]=a*cphi*cE*Edot;
    v.R[2]=0.0;
    
    *vp = Mv(R, v);
}

double M_to_E (double e, double M)
{
    double guess=M;
    double s=sin(guess);
    double kEq=guess-e*s-M;
    
    while (kEq>K_PRECISION || kEq<-K_PRECISION){
      double c=cos(guess);
      double refineN = -kEq/(1.0-e*c);
      double refineH = -kEq/(1.0-e*c+0.5*s*refineN);
      guess += refineH;
      s=sin(guess);
      kEq=guess-e*s-M;
    }
    return guess;
}
//----------------------------------------------------------------


//---------------------- kepler evolver---------------------------
void k_evolve(VEC* xp, VEC* vp, double dt)
{   
    double getDeltaE(double,double,double);
    VEC r_i,v_i;
    double DeltaM;
    double rnorm,vsquare,u;
    double a,n,ecosE,esinE;
    double DeltaE;
    double s2,c2,s,c,fp;
    double f,g,fdot,gdot;
    double flag_minusenergy;
    int j;
    
    for (j=0; j<NJ; j++) {
      // STEP 0:
      r_i = *(xp+j);
      v_i = *(vp+j);
        
      rnorm = vlen(r_i);
      vsquare = in_p(v_i, v_i);
      u = in_p(r_i, v_i);
      
      flag_minusenergy=2.0/rnorm-vsquare/kconst[j];
      if (flag_minusenergy<=0.0){
	printf("\nk_evolve: UNBOUND ORBIT.\n");
	flag_unbound_orb = 1;
	raise(SIGINT);
      }
      
      a=1.0/flag_minusenergy;
      n=sqrt(kconst[j]/(a*a*a));
      ecosE=1.0-rnorm/a;
      esinE=u/(n*a*a);
      DeltaM=n*dt;
      
      // STEP 1:
      DeltaE=getDeltaE(ecosE, esinE, DeltaM);
      
      // STEP 2:
      // to get the f,g, fdot, gdot values  for this given initial state and time increase.
      //       -- 0. get the s2,c2,s,c,fp values.   (fp means f^prime )
      //       -- 1. get the f,g,fdot,gdot values.
      s2=sin(0.5*DeltaE);
      c2=cos(0.5*DeltaE);
      c=c2*c2-s2*s2;
      s=2.0*s2*c2;
      fp=1.0-c*ecosE+s*esinE;
      
      f=1.0-2.0*s2*s2*a/rnorm;
      g=2.0*s2*(s2*esinE+c2*rnorm/a)/n;
      fdot=-1.0*n*(a/rnorm)*s/fp;
      gdot=1.0-2.0*s2*s2/fp;
      
      // STEP 3:
      *(xp+j) = vsum( cv(f,r_i), cv(g,v_i) );
      *(vp+j) = vsum( cv(fdot,r_i), cv(gdot,v_i) );
    }
}

int sign(double x)
{
  if (x>0.0)
    return(1);
  else if(x< 0.0)
    return(-1);
  else
    return(0);
}

double increKeq_refineD(double ecosE,double esinE, double DeltaM, double guess)
{
    double c=cos(guess), s=sin(guess);
    double d1,dd,d3;
    double RN,RH,RD;
    d1=1+s*esinE-c*ecosE;
    dd=c*esinE+s*ecosE;
    d3=-s*esinE+c*ecosE;
    
    double Eq=guess+(1.0-c)*esinE-s*ecosE-DeltaM;
    RN= -1.0*Eq/d1;
    RH= -1.0*Eq/(d1+0.5*dd*RN);
    RD= -1.0*Eq/(d1+0.5*dd*RH+1.0/6.0*d3*RH*RH);
    return RD;
}

double getDeltaE(double ecosE, double esinE, double DeltaM)
{
    double guess=DeltaM/(1.0-ecosE);    // JW's guess.
    double refineD=increKeq_refineD(ecosE, esinE, DeltaM, guess);
    int count=1;

    while ( fabs(refineD)>K_PRECISION && count<5) {
	if (count==1)
	  guess += refineD;
        refineD=increKeq_refineD(ecosE, esinE, DeltaM, guess);
        guess += refineD;
        count++;
    }
    
    if (count>=5) {
        double y_inter,e_inter;
        int sig_inter;
        y_inter=DeltaM-esinE;
        sig_inter=sign(esinE*cos(y_inter)+ecosE*sin(y_inter));
        e_inter=sqrt(ecosE*ecosE+esinE*esinE);
        guess=y_inter+0.85*sig_inter*e_inter;  // danby_guess.
        
        while (fabs(refineD)>K_PRECISION) {
	  refineD=increKeq_refineD(ecosE, esinE, DeltaM, guess);
	  guess += refineD;
        }
    }
    return guess;    
}
//----------------------------------------------------------------

const double TWOPI = M_PI +M_PI;

double atan2_pvTWOPI(double y, double x)
{
  // to [0,2pi);
  double pv = atan2(y,x);
  if (pv < 0.)
    pv += TWOPI;
  return pv;
}
