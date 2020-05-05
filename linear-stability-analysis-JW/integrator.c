/*  Bulirsch-Stoer integrator */

int integrate(Istate *state, double *official_delta_t)
{
  Istate local_state;
  Istate_Deriv state_deriv;
  static int lt[10] = {1,2,3,4,6,8,12,16,24,32};
  static double temp[MAX_NUMBER_EQUATIONS][12];
  double d[6], delta_t;
  int equations_of_motion();

  int i, m, m1, k, i1max, i1;
  double xa, xb, varm, h, hd, eta2, dta;
  double yb, c, b1, den, dtn, b, var, varma;
  double foo, flt;

  d[0] = 0.0;
  d[1] = 0.0;
  d[2] = 0.0;
  d[3] = 0.0;
  d[4] = 0.0;
  d[5] = 0.0;

 start:

  local_state.t = state->t;
  for(i=0;i<NUMBER_EQUATIONS;i++) {
    local_state.x[i] = state->x[i];
  }
  delta_t = *official_delta_t;

  for(i=0;i<NUMBER_EQUATIONS;i++) {
    for(m=0; m<12; m++) {
      temp[i][m] = 0.;
    }
  }

  xa = local_state.t;

  if(equations_of_motion(&local_state, &state_deriv) == FAILURE) {
    return(FAILURE);
  }
 
  for(i=0; i<NUMBER_EQUATIONS; i++)
    {
      temp[i][1] = ABS(local_state.x[i]);
      if(temp[i][1] < INTEGRATION_EPSILON)
	temp[i][1] = INTEGRATION_EPSILON;
      temp[i][4] = state_deriv.x[i];
      temp[i][0] = local_state.x[i];
    }

 outside:
  xb = delta_t + xa;
  for(m=0; m<10; m++)
    {
      flt = lt[m];
      varm = 0.0;
      m1 = MIN(m, 6);
      if(m1 != 0)
	  for(k=0; k<m1; k++)
	    {
	      d[k] = flt/lt[m-k-1];
	      d[k] *= d[k];
	    }
      h = delta_t / flt;
      hd = 0.5 * h;
      for(i=0; i<NUMBER_EQUATIONS; i++)
	{
	  temp[i][3] = temp[i][0];
	  local_state.x[i] = temp[i][0] + hd*temp[i][4];
	}
      i1max = 2*flt - 1;
      local_state.t = xa;
      for(i1=0; i1<i1max; i1++)
	{
	  local_state.t += hd;

  	  if(equations_of_motion(&local_state, &state_deriv) == FAILURE) {
	    return(FAILURE);
	  }

	  for(i=0; i<NUMBER_EQUATIONS; i++)
	    {
	      foo = local_state.x[i];
	      temp[i][1] = MAX(temp[i][1], ABS(foo));
	      if((temp[i][1]) == (0.0))
		{
		  printf("difsys: temp[i][1] = 0\n");
		  return(FAILURE);
		}
	      eta2 = temp[i][3] + h*state_deriv.x[i];
	      temp[i][3] = local_state.x[i];
	      local_state.x[i] = eta2;
	    }
	}

      local_state.t = xb;

      if(equations_of_motion(&local_state, &state_deriv) == FAILURE) {
	return(FAILURE);
      }

      for(i=0; i<NUMBER_EQUATIONS; i++)
	{
	  dta = temp[i][11];
	  yb = (temp[i][3] + local_state.x[i] 
		   + hd*state_deriv.x[i]) / 2.0;
	  c = yb;
	  temp[i][11] = yb;
	  if(m1 != 0)
	    {
	      for(k=0; k<m1; k++)
		{
		  b1 = d[k] * dta;
		  den = b1 - c;
		  dtn = dta;
		  if(den != 0)
		    {
		      b = (c - dta) / den;
		      dtn = c * b;
		      c = b1 * b;
		    }
		  dta = temp[i][11-k-1];
		  temp[i][11-k-1] = dtn;
		  yb += dtn;
		}
	      var = ABS(temp[i][2] - yb) / temp[i][1];
	      varm = MAX(varm, var);
	    }
	  temp[i][2] = yb;
	}
      
      if(m < 3)
	varma = varm; 
      else if(varm <= INTEGRATION_EPSILON) /* check for convergence or divergence */
	goto end;
      else if(varm >= varma)
	break;
    }
  
  delta_t /= 2.0;

  goto outside;
  
 end: local_state.t = xb;
  for(i=0; i<NUMBER_EQUATIONS; i++)
    local_state.x[i] = temp[i][2];

/* reset state variables */
  state->t = local_state.t; 
  for(i=0;i<NUMBER_EQUATIONS;i++)
    state->x[i] = local_state.x[i];
  *official_delta_t = delta_t * 1.5 * pow(0.6, (double)(m-m1));

  return(SUCCESS);
}


