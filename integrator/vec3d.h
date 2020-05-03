//dimension is 3;

typedef struct{
  double R[3];
}COL; //column;

typedef struct{
  COL C[3];
}SQM;

typedef COL VEC;

SQM MN(SQM a, SQM b) // matrix product of two square matrices;
{
  SQM c; //c.C[j].R[i] = a.C[0].R[i]*b.C[j].R[0] +a.C[1].R[i]*b.C[j].R[1] +a.C[2].R[i]*b.C[j].R[2];
  c.C[0].R[0] = a.C[0].R[0]*b.C[0].R[0] +a.C[1].R[0]*b.C[0].R[1] +a.C[2].R[0]*b.C[0].R[2];
  c.C[0].R[1] = a.C[0].R[1]*b.C[0].R[0] +a.C[1].R[1]*b.C[0].R[1] +a.C[2].R[1]*b.C[0].R[2];
  c.C[0].R[2] = a.C[0].R[2]*b.C[0].R[0] +a.C[1].R[2]*b.C[0].R[1] +a.C[2].R[2]*b.C[0].R[2];

  c.C[1].R[0] = a.C[0].R[0]*b.C[1].R[0] +a.C[1].R[0]*b.C[1].R[1] +a.C[2].R[0]*b.C[1].R[2];
  c.C[1].R[1] = a.C[0].R[1]*b.C[1].R[0] +a.C[1].R[1]*b.C[1].R[1] +a.C[2].R[1]*b.C[1].R[2];
  c.C[1].R[2] = a.C[0].R[2]*b.C[1].R[0] +a.C[1].R[2]*b.C[1].R[1] +a.C[2].R[2]*b.C[1].R[2];

  c.C[2].R[0] = a.C[0].R[0]*b.C[2].R[0] +a.C[1].R[0]*b.C[2].R[1] +a.C[2].R[0]*b.C[2].R[2];
  c.C[2].R[1] = a.C[0].R[1]*b.C[2].R[0] +a.C[1].R[1]*b.C[2].R[1] +a.C[2].R[1]*b.C[2].R[2];
  c.C[2].R[2] = a.C[0].R[2]*b.C[2].R[0] +a.C[1].R[2]*b.C[2].R[1] +a.C[2].R[2]*b.C[2].R[2];

  return c;
}

COL Mv(SQM m, COL v0){
  COL v1; //v1.R[i] = m.C[0].R[i]*v0.R[0] +m.C[1].R[i]*v0.R[1] +m.C[2].R[i]*v0.R[2];
  v1.R[0] = m.C[0].R[0]*v0.R[0] +m.C[1].R[0]*v0.R[1] +m.C[2].R[0]*v0.R[2];
  v1.R[1] = m.C[0].R[1]*v0.R[0] +m.C[1].R[1]*v0.R[1] +m.C[2].R[1]*v0.R[2];
  v1.R[2] = m.C[0].R[2]*v0.R[0] +m.C[1].R[2]*v0.R[1] +m.C[2].R[2]*v0.R[2];
  return v1;
}
COL trMv(SQM m, COL v0){
  COL v1; //v1.R[i] = m.C[i].R[0]*v0.R[0] +m.C[i].R[1]*v0.R[1] +m.C[i].R[2]*v0.R[2];
  v1.R[0] = m.C[0].R[0]*v0.R[0] +m.C[0].R[1]*v0.R[1] +m.C[0].R[2]*v0.R[2];
  v1.R[1] = m.C[1].R[0]*v0.R[0] +m.C[1].R[1]*v0.R[1] +m.C[1].R[2]*v0.R[2];
  v1.R[2] = m.C[2].R[0]*v0.R[0] +m.C[2].R[1]*v0.R[1] +m.C[2].R[2]*v0.R[2];
  return v1;
}

void v_cast(VEC* p, double x, double y, double z)
{
  p->R[0] = x;
  p->R[1] = y;
  p->R[2] = z;
}

COL cv(double c, COL input)
{
  COL Final;
  Final.R[0] = c*input.R[0];
  Final.R[1] = c*input.R[1];
  Final.R[2] = c*input.R[2];
  return Final;
}

COL vsum(COL a, COL b)
{
  COL v;
  v.R[0] = a.R[0] + b.R[0];
  v.R[1] = a.R[1] + b.R[1];
  v.R[2] = a.R[2] + b.R[2];
  return v;
}

COL vdif(COL end, COL start)
{
  COL v;
  v.R[0] = end.R[0] -start.R[0];
  v.R[1] = end.R[1] -start.R[1];
  v.R[2] = end.R[2] -start.R[2];
  return v;
}

void vgrow(COL* ap, COL b)
{
  ap->R[0] += b.R[0];
  ap->R[1] += b.R[1];
  ap->R[2] += b.R[2];
}

SQM R_a(double g)  //rotation about x (space ref.) or a (body ref.);
{
  SQM m;
  double c,s;
  c = cos(g);
  s = sin(g);

  m.C[0].R[0] = 1.0;
  m.C[0].R[1] = 0.0;
  m.C[0].R[2] = 0.0;

  m.C[1].R[0] = 0.0;
  m.C[1].R[1] = c;
  m.C[1].R[2] = s;

  m.C[2].R[0] = 0.0;
  m.C[2].R[1] = -s;
  m.C[2].R[2] = c;

  return m;
}


SQM R_b(double g)
{
  SQM m;
  double c,s;
  c = cos(g);
  s = sin(g);

  m.C[0].R[0] = c;
  m.C[0].R[1] = 0.0;
  m.C[0].R[2] = s;

  m.C[1].R[0] = 0.0;
  m.C[1].R[1] = 1.0;
  m.C[1].R[2] = 0.0;

  m.C[2].R[0] = -s;
  m.C[2].R[1] = 0.0;
  m.C[2].R[2] = c;

  return m;
}


SQM R_c(double g)  //rotation about z (space ref.) or c (body ref.);
{
  SQM m;
  double c,s;
  c = cos(g);
  s = sin(g);

  m.C[0].R[0] = c;
  m.C[0].R[1] = s;
  m.C[0].R[2] = 0.0;

  m.C[1].R[0] = -s;
  m.C[1].R[1] = c;
  m.C[1].R[2] = 0.0;

  m.C[2].R[0] = 0.0;
  m.C[2].R[1] = 0.0;
  m.C[2].R[2] = 1.0;

  return m;
}


SQM R_v(COL k, double g) //rotation about vector k;
{
  SQM sqm_sum_aA_bB_I(double, SQM, double, SQM), MN(SQM,SQM);

  double ksq = k.R[0]*k.R[0] + k.R[1]*k.R[1] + k.R[2]*k.R[2];
  double c0 = cos(g);
  double s0 = sin(g);

  SQM S; //antisym form of k;
  S.C[0].R[0] = 0.;
  S.C[0].R[1] = k.R[2];
  S.C[0].R[2] = -k.R[1];

  S.C[1].R[0] = -k.R[2];
  S.C[1].R[1] = 0.;
  S.C[1].R[2] = k.R[0];

  S.C[2].R[0] = k.R[1];
  S.C[2].R[1] = -k.R[0];
  S.C[2].R[2] = 0.;

  return sqm_sum_aA_bB_I(s0/sqrt(ksq), S, (1.-c0)/ksq, MN(S,S));
}
SQM ezR_v(COL k, double c_g, double s_g)
{
  SQM sqm_sum_aA_bB_I(double, SQM, double, SQM), MN(SQM,SQM);

  double ksq = k.R[0]*k.R[0] + k.R[1]*k.R[1] + k.R[2]*k.R[2];

  SQM S; //antisym form of k;
  S.C[0].R[0] = 0.;
  S.C[0].R[1] = k.R[2];
  S.C[0].R[2] = -k.R[1];

  S.C[1].R[0] = -k.R[2];
  S.C[1].R[1] = 0.;
  S.C[1].R[2] = k.R[0];

  S.C[2].R[0] = k.R[1];
  S.C[2].R[1] = -k.R[0];
  S.C[2].R[2] = 0.;

  return sqm_sum_aA_bB_I(s_g/sqrt(ksq), S, (1.-c_g)/ksq, MN(S,S));
}

SQM sqm_sum_aA_bB_I(double ca, SQM A, double cb, SQM B) //ca*A+cb*B+I;
{
  SQM C;
  C.C[0].R[0] = ca*A.C[0].R[0] +cb*B.C[0].R[0] +1.;
  C.C[0].R[1] = ca*A.C[0].R[1] +cb*B.C[0].R[1];
  C.C[0].R[2] = ca*A.C[0].R[2] +cb*B.C[0].R[2];

  C.C[1].R[0] = ca*A.C[1].R[0] +cb*B.C[1].R[0];
  C.C[1].R[1] = ca*A.C[1].R[1] +cb*B.C[1].R[1] +1.;
  C.C[1].R[2] = ca*A.C[1].R[2] +cb*B.C[1].R[2];

  C.C[2].R[0] = ca*A.C[2].R[0] +cb*B.C[2].R[0];
  C.C[2].R[1] = ca*A.C[2].R[1] +cb*B.C[2].R[1];
  C.C[2].R[2] = ca*A.C[2].R[2] +cb*B.C[2].R[2] +1.;
  return C;
}


void print_sqm(SQM a)
{
  int i,j;
  printf("\n-------------------------------\n");
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++)
      printf("%9.2le  ",a.C[j].R[i]);
    printf("\n");
  }
  printf("-------------------------------\n\n");
}


void print_col(COL v)
{
  int i;
  printf("\n----------\n");
  for(i=0; i<3; i++) {
    printf("%9.2le  ",v.R[i]);
    printf("\n");
  }
  printf("----------\n\n");
}



double in_p(COL a, COL v_b)
{
  return ( a.R[0]*v_b.R[0] + a.R[1]*v_b.R[1] + a.R[2]*v_b.R[2] );
}

COL cr_p(COL x, COL v)
{
  COL h;
  h.R[0] = x.R[1]*v.R[2] -v.R[1]*x.R[2];
  h.R[1] = x.R[2]*v.R[0] -v.R[2]*x.R[0];
  h.R[2] = x.R[0]*v.R[1] -v.R[0]*x.R[1];
  return h;
}

COL vhat(COL a)
{
  double len=sqrt(a.R[0]*a.R[0] +a.R[1]*a.R[1] +a.R[2]*a.R[2]);
  COL hat;
  hat.R[0] = a.R[0]/len;
  hat.R[1] = a.R[1]/len;
  hat.R[2] = a.R[2]/len;
  return hat;
}

double vlen(COL a)
{
  return sqrt(a.R[0]*a.R[0] +a.R[1]*a.R[1] +a.R[2]*a.R[2]);
}

COL xhat(void)
{
  COL a;
  a.R[0]=1.0;
  a.R[1]=0.0;
  a.R[2]=0.0;
  return(a);
}

COL yhat(void)
{
  COL b;
  b.R[0]=0.0;
  b.R[1]=1.0;
  b.R[2]=0.0;
  return(b);
}

COL zhat(void)
{
  COL c;
  c.R[0]=0.0;
  c.R[1]=0.0;
  c.R[2]=1.0;
  return(c);
}
