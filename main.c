#define TDOMAX 10
#define MAP_FREQ 40
#define SPEEDUP 1
#define WOBBLE_DAMP 1

#define _GNU_SOURCE
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include "vec3d.h"
#include "2b.h"
#include "rb.h"
#include "lag.h"
#include "symap.h"

FILE * dataf;
FILE * resu_info;

char run_id[300];

SS state;

void check_print(SS* sp)
{
  FILE* file_l=fopen("check_modify.txt","w");
  fprintf(file_l,"xJ0 %.15le %.15le %.15le\n",sp->x[0].R[0],sp->x[0].R[1],sp->x[0].R[2]);
  fprintf(file_l,"xJ1 %.15le %.15le %.15le\n",sp->x[1].R[0],sp->x[1].R[1],sp->x[1].R[2]);
  fprintf(file_l,"vJ0 %.15le %.15le %.15le\n",sp->v[0].R[0],sp->v[0].R[1],sp->v[0].R[2]);
  fprintf(file_l,"vJ1 %.15le %.15le %.15le\n",sp->v[1].R[0],sp->v[1].R[1],sp->v[1].R[2]);
  fprintf(file_l,"earth-Lp %.15le %.15le %.15le\n",sp->Lpe.R[0],sp->Lpe.R[1],sp->Lpe.R[2]);
  fprintf(file_l,"moon--Lp %.15le %.15le %.15le\n",sp->Lpm.R[0],sp->Lpm.R[1],sp->Lpm.R[2]);
  fclose(file_l);
  printf("continue checking ? -- (\"modify\")\n");
  int go_on;
  scanf("%d",&go_on);
  if(go_on==0)
    exit(-1);
}

void read_ic(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
void save_ic(double, double, double, double, double, double, double, double, double, double);
void read_meta(double, double*, double*, int*, int*);
void save_meta(double, double, int, int);
void rewrite_meta(int);

void sig_handler(int signum)
{
  rewrite_meta(signum);

  fclose(dataf);
  fclose(resu_info);

  signal(signum,SIG_DFL);
  raise(signum);
}

void set_up(SS* sp, AUX* auxp, double* a_p, double* amin_p, double* amax_p, int* arraylen_p, int* arraynum_p)
{
  double ini_a, ini_daylength, ini_Earth_I, ini_eps, ini_Moon_I, ini_e_lunar, Qe, k2e, Qm, k2m;

  set_kepler_para();
  read_ic(&ini_a, &ini_daylength, &ini_Earth_I, &ini_eps, &ini_Moon_I, &ini_e_lunar, &Qe, &k2e, &Qm, &k2m);
  
  cos_oQ[E_DEFOR] = cos(1./Qe);
  sin_oQ[E_DEFOR] = sin(1./Qe);
  cos_oQ[M_DEFOR] = cos(1./Qm);
  sin_oQ[M_DEFOR] = sin(1./Qm);
  c_k2[E_DEFOR] = k2e*pow(r_e, 2*TIDE_L+1)*G_const*mass_m;
  c_k2[M_DEFOR] = k2m*pow(r_m, 2*TIDE_L+1)*G_const*mass_e;
  dampbase_m = k2m *pow(r_m, 5) /(3.0*G_const*Cm) /Qm;

  ini_nb(sp->x, sp->v, ini_a, ini_e_lunar, ini_Earth_I, ini_eps);

  e_flatten(TWOPI/ini_daylength *TWOPI/ini_daylength, &auxp->Ae, &auxp->Ce);
  ini_Earth(ini_daylength, auxp->Ce, ini_Earth_I, &sp->Lpe, &sp->Oe);
  double a__, foo;
  xv2elem_invar(kconst[0], sp->x+0, sp->v+0, &foo, &a__, &foo, &foo, &foo);
  ini_Moon(sqrt(kconst[0]/a__/a__/a__), ini_Earth_I+ini_eps+ini_Moon_I, &sp->Lpm, &sp->Om);

  map_period(a__, auxp);
  *a_p = a__;

  double amin_re, amax_re;
  read_meta(ini_a, &amin_re, &amax_re, arraylen_p, arraynum_p);
  *amin_p = amin_re *r_e;
  *amax_p = amax_re *r_e;

  char path[300];
  sprintf(path,"./data/evol_%s.out",run_id);
  if(access(path, F_OK)==0) {printf("%s already exists.\n\n", path); exit(-1);}
  dataf = fopen(path,"wb");
  sprintf(path,"./data/resuinfo_%s.out",run_id);
  resu_info = fopen(path,"wb");

  save_ic(ini_a, ini_daylength, ini_Earth_I, ini_eps, ini_Moon_I, ini_e_lunar, Qe, k2e, Qm, k2m);
  save_meta(amin_re, amax_re, *arraylen_p, *arraynum_p);
}
void read_ic(double* a_p, double* daylen_p, double* Ie_p, double* eps_p, double* Im_p, double* e_p, double* Qe_p, double* k2e_p, double* Qm_p, double* k2m_p)
{
  *Qe_p = 100.0;
  *Qm_p = 49.934064676;
  *k2e_p = 0.299;
  *k2m_p = 0.03; // Jim Williams, Dicky;

  *a_p = 3.5*r_e;
  *daylen_p = 2.5/24.0; //(2.5/24.0)[day], or 2.5[hr];
  *Ie_p = 0.0;
  *eps_p = 0.0;
  *Im_p = 0.0;
  *e_p = 0.001;

  int index;
  double val, Apara, Ae, Ce, L_, Ls = Ce_tp *sqrt(kconst[0]/r_e/r_e/r_e);
  int set_up_correct;

  printf("\nthe default run set-up is");
  printf("\n1         2         3         4         5         6         7         8         9         10        ");
  printf("\nQ_Moon    k2m       Q_Earth   k2e       a_in_RE   e         dayl(hr)  I_e(deg)  epsi(deg) I_m(deg)");
  printf("\n%-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf",
	 *Qm_p, *k2m_p, *Qe_p, *k2e_p,
	 *a_p/r_e, *e_p, *daylen_p*24.0, *Ie_p*180.0/M_PI, *eps_p*180.0/M_PI, *Im_p*180.0/M_PI);
  double me_mm_sq = mlis[0]*mlis[0]/mlis[1]/mlis[1];
  Apara = *k2m_p/ *Qm_p /(*k2e_p/ *Qe_p) *me_mm_sq *pow(r_m/r_e,5.0);
  L_ = mr *sqrt(kconst[0]**a_p *(1.-*e_p**e_p));
  e_flatten(TWOPI/ *daylen_p *TWOPI/ *daylen_p, &Ae, &Ce);
  L_ += Ce *TWOPI/ *daylen_p;
  L_ += Cm *sqrt(kconst[0]/ *a_p/ *a_p/ *a_p);
  printf("\n(this gives A=%6.3lf, scalar sum L/Ls=%6.3lf.)",Apara,L_/Ls);

 modify_set_up:
  printf("\n\nwhich parameter should be reset for this run? (para-order, new-value) (end input: 0, 0.0)\n");

  scanf("%d %lf",&index,&val);
  while(index != 0){
    if(index==1) *Qm_p = val;
    else if(index==2) *k2m_p = val;
    else if(index==3) *Qe_p = val;
    else if(index==4) *k2e_p = val;
    else if(index==5) *a_p = val*r_e; // "val" in [RE]; transform;
    else if(index==6) *e_p = val;
    else if(index==7) *daylen_p = val/24.0; // "val" in [hr]; transform;
    else if(index==8) *Ie_p = val/180.0*M_PI; // "val" in [deg]; transform;
    else if(index==9) *eps_p = val/180.0*M_PI; // "val" in [deg]; transform;
    else if(index==10) *Im_p = val/180.0*M_PI; // "val" in [deg]; transform;
    else printf("\npara-order not known\n\n");
    scanf("%d %lf",&index,&val);
  }

  printf("\nthe new run set-up is");
  printf("\n1         2         3         4         5         6         7         8         9         10        ");
  printf("\nQ_Moon    k2m       Q_Earth   k2e       a_in_RE   e         dayl(hr)  I_e(deg)  epsi(deg) I_m(deg)");
  printf("\n%-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf %-9.3lf",
	 *Qm_p, *k2m_p, *Qe_p, *k2e_p,
	 *a_p/r_e, *e_p, *daylen_p*24.0, *Ie_p*180.0/M_PI, *eps_p*180.0/M_PI, *Im_p*180.0/M_PI);
  Apara = *k2m_p/ *Qm_p /(*k2e_p/ *Qe_p) *me_mm_sq *pow(r_m/r_e,5.0);
  L_ = mr *sqrt(kconst[0]**a_p *(1.-*e_p**e_p));
  e_flatten(TWOPI/ *daylen_p *TWOPI/ *daylen_p, &Ae, &Ce);
  L_ += Ce *TWOPI/ *daylen_p;
  L_ += Cm *sqrt(kconst[0]/ *a_p/ *a_p/ *a_p);
  printf("\n(this gives A=%6.3lf, scalar sum L/Ls=%6.3lf.)",Apara,L_/Ls);

  printf("\ncorret? (1-yes, 0-no)\n");
  scanf("%d",&set_up_correct);
  if( set_up_correct != 1 ) goto modify_set_up;
}
void read_meta(double ini_a, double* amin_re_p, double* amax_re_p, int* arraylen_p, int* arraynum_p)
{
  *amax_re_p = 15.0;
  *amin_re_p = 1.0;
  *arraylen_p = 10000;
  *arraynum_p = 20000000;

  int index;
  double t_kyr;
  size_t datasize_MB;
  AUX aux;
  map_period(ini_a,&aux);
  char valstr[50];
  int ctrl_para_correct;

  t_kyr = aux.dt*(double)*arraylen_p*(double)*arraynum_p /365.25/1000.0 *(double)SPEEDUP;
  datasize_MB = sizeof(SS)**arraynum_p/1024/1024;
  printf("\ncomputation control:");
  printf("\n1         2         3         4          ");
  printf("\na_min(re) a_max(re) leap_l    leap_num   ");
  printf("\n%-9.3lf %-9.3lf %-9d %-9d",*amin_re_p,*amax_re_p,*arraylen_p,*arraynum_p);
  //estimated final time [kyr]
  //estimated data file size.
  printf("\nestimated final time %.0le kyr, datafile size %zu MB.\n",t_kyr,datasize_MB);

 comp_ctrl_change:
  printf("\n\nwhich ctrl parameter should be reset for this run? (para-order, new-value) (end input: 0, 0.0)\n");
  scanf("%d %s",&index,valstr);
  while(index != 0){
    if(index==1) sscanf(valstr,"%lf",amin_re_p);
    else if(index==2) sscanf(valstr,"%lf",amax_re_p);
    else if(index==3) sscanf(valstr,"%d",arraylen_p);
    else if(index==4) sscanf(valstr,"%d",arraynum_p);
    else printf("\npara-order not known\n\n");
    scanf("%d %s",&index,valstr);
  }

  t_kyr = aux.dt*(double)*arraylen_p*(double)*arraynum_p /365.25/1000.0 *(double)SPEEDUP;
  datasize_MB = sizeof(SS)**arraynum_p/1024/1024;
  printf("\nnew computation control:");
  printf("\n1         2         3         4          ");
  printf("\na_min(re) a_max(re) leap_l    leap_num   ");
  printf("\n%-9.3lf %-9.3lf %-9d %-9d",*amin_re_p,*amax_re_p,*arraylen_p,*arraynum_p);
  printf("\nestimated final time %.0le kyr, datafile size %zu MB.\n",t_kyr,datasize_MB);

  printf("\ncorrect? (1-yes, 0-no)\n");
  scanf("%d",&ctrl_para_correct);
  if( ctrl_para_correct != 1 ) goto comp_ctrl_change;
}

int main(int argc, const char * argv[])
{
  if(argc==2){
    sprintf(run_id,"%s",argv[1]);
  } else{
    printf("\nargument list: run index\n\n");
    exit(-1);
  }

  /*ic, para, and meta (ctrl)*/
  AUX aux;
  double amin, amax, a__;
  int arraylen, arraynum;
  set_up(&state, &aux, &a__, &amin, &amax, &arraylen, &arraynum);
  signal(SIGINT,sig_handler);
  signal(SIGTERM,sig_handler);

  /*pre-leap preparation*/
  map_phase_adjust(&state, &aux);

  /*beginning leap*/
  int i = 0;
  double a_mark = a__;
  while (i<arraynum && a__<=amax && a__>=amin) {
    if(a__ >= a_mark*1.59 || a__ <= a_mark/1.59){
      map_period(a__, &aux);
      a_mark = a__;
    }
    fwrite(&state, sizeof(SS), 1, dataf);
    if(i%100==0){
      fwrite(&state, sizeof(SS), 1, resu_info);
      fwrite(&aux, sizeof(AUX), 1, resu_info);
      fwrite(&a_mark, sizeof(a_mark), 1, resu_info);
    }
    map_array(&state, &aux, arraylen);
    double foo;
    xv2elem_invar(kconst[0], state.x+0, state.v+0, &foo, &a__, &foo, &foo, &foo);
    i++;
  }

  fclose(dataf);
  fclose(resu_info);

  if(i>=arraynum)
    rewrite_meta(0);
  else if(a__>amax)
    rewrite_meta(-1);
  else
    rewrite_meta(-2);

  return 0;
}


void save_ic(double ini_a, double ini_daylength, double ini_Earth_I, double ini_eps, double ini_Moon_I, double ini_e_lunar,
	     double Qe, double k2e, double Qm, double k2m)
{
  char path[500];
  sprintf(path,"./data/ic_%s.txt", run_id);
  FILE* fic=fopen(path,"w");

  fprintf(fic,"\nrun index: %s\n", run_id);

  fprintf(fic,"\nQ_Earth, k2e: %9.5lf %9.5lf\n",Qe,k2e);
  fprintf(fic,"Q_Moon, k2m: %9.5lf %9.5lf\n",Qm,k2m);
  double me_mm_sq = mlis[0]*mlis[0]/mlis[1]/mlis[1];
  fprintf(fic,"A parameter %9.5lf\n", k2m/Qm /(k2e/Qe) *me_mm_sq *pow(r_m/r_e,5.0) );

  fprintf(fic,"\ninitial configuration\na %6.2lf (Earth radii)\n",ini_a/r_e);
  fprintf(fic,"e %.3le\n",ini_e_lunar);
  fprintf(fic,"daylength %6.2lf (hours)\n",ini_daylength*24.0);
  fprintf(fic,"Earth's obliquity %6.2lf (deg)\n",ini_Earth_I*180.0/M_PI);
  fprintf(fic,"mutual obliquity (moon orbit to earth equator) %6.2lf (deg)\n",ini_eps*180.0/M_PI);
  fprintf(fic,"Moon's obliquity %6.2lf (deg)\n",ini_Moon_I*180.0/M_PI);

  fprintf(fic,"\n");
  fclose(fic);
}
void save_meta(double amin_re, double amax_re, int arraylen, int arraynum)
{
  char path[500];
  sprintf(path,"./data/meta_%s.txt",run_id);
  FILE* fmeta=fopen(path,"w");

  fprintf(fmeta,"\nrun index: %s\n",run_id);

  fprintf(fmeta,"\ncomputation control parameters\n");
  fprintf(fmeta,"a_min(R_earth) %lf   a_max(R_earth) %lf\n",amin_re,amax_re);
  fprintf(fmeta,"leap length %d   leap num %d\n",arraylen,arraynum);

  fprintf(fmeta,"\nL/Ls\n");
  double Ls = Ce_tp *sqrt(kconst[0]/r_e/r_e/r_e);
  double Lem0 = angmom(&state);
  fprintf(fmeta,"when run starts %9.2lf\n",Lem0/Ls);

  fprintf(fmeta,"\nrun time\n");
  time_t time0;
  time(&time0);
  struct tm *tm_p = gmtime(&time0);
  fprintf(fmeta,"UTC    %4d-%02d-%02d-%02d:%02d:%02d\n",
	  1900+tm_p->tm_year, 1+tm_p->tm_mon, tm_p->tm_mday, tm_p->tm_hour, tm_p->tm_min, tm_p->tm_sec);

  fclose(fmeta);
}
void rewrite_meta(int signum)
{
  char str0[500] = "when run starts"; //L/Ls at start;
  char str1[500] = "UTC";             //time at start;

  char oldpath[500], newpath[500];
  sprintf(oldpath,"./data/meta_%s.txt",run_id);
  sprintf(newpath,"./data/tmpmeta_%s.txt",run_id);

  FILE* oldmeta=fopen(oldpath,"r"), * newmeta=fopen(newpath,"w");
  {
    char* line = NULL;
    size_t len = 0;
    ssize_t nread;
    while((nread = getline(&line, &len, oldmeta)) != -1){
      if(strncmp(line, str0, strlen(str0))==0){
	//line of L/Ls at start;
	fprintf(newmeta, "%s", line);
	double Lem1 = angmom(&state);
	double Ls = Ce_tp *sqrt(kconst[0]/r_e/r_e/r_e);
	fprintf(newmeta, "when run ends   %9.2lf\n",Lem1/Ls);

      } else if(strncmp(line, str1, strlen(str1))==0){
	//line of time at start;
	time_t time0;
	struct tm tm;
	strptime(line, "UTC    %F-%T", &tm);
	time0 = timegm(&tm);
	char* local_t0 = ctime(&time0);
	fprintf(newmeta, "starts %s %s",tzname[0],local_t0);
	fprintf(newmeta, "%s", line);

	time_t time1;
	time(&time1);
	fprintf(newmeta, "ends   %s %s",tzname[0],ctime(&time1));
	int elapse = difftime(time1,time0);
	fprintf(newmeta, "time spent in the run %d seconds\n",elapse);
	int elapsH = elapse/3600;
	elapse -= 3600*elapsH;
	int elapsM = elapse/60;
	elapse -= 60*elapsM;
	fprintf(newmeta, "(%2d hours %2d minutes %2d seconds)\n",elapsH,elapsM,elapse);
	if(flag_unbound_orb==1)
	  fprintf(newmeta, "UNBOUND ORBIT\n");
	else if(signum==SIGINT)
	  fprintf(newmeta, "SIGINT\n");
	else if(signum==SIGTERM)
	  fprintf(newmeta, "SIGTERM\n");
	else if(signum==0)
	  fprintf(newmeta, "COMPLETION-ARRAYNUM\n");
	else if(signum==-1)
	  fprintf(newmeta, "COMPLETION-SEMIMAJOR-AXIS-MAX\n");
	else if(signum==-2)
	  fprintf(newmeta, "COMPLETION-SEMIMAJOR-AXIS-MIN\n");
	fprintf(newmeta, "\n");

	char str_host[500];
	gethostname(str_host, sizeof(str_host));
	fprintf(newmeta, "theater %s\n", str_host);
	fprintf(newmeta, "\n");

      } else{
	fprintf(newmeta, "%s", line);
      }
    }
  }
  fclose(oldmeta);
  fclose(newmeta);

  remove(oldpath);
  rename(newpath, oldpath);
}

