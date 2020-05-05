#define MAX_NUMBER_EQUATIONS 5

typedef struct {
  double t;
  double x[MAX_NUMBER_EQUATIONS];
} Istate;

typedef struct {
  double t;
  double x[MAX_NUMBER_EQUATIONS];
} Istate_Deriv;

#define SUCCESS 1
#define FAILURE 2

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
