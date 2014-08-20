#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP window_mean(SEXP X, SEXP times_start, SEXP times_end)
{
	int i,j,k,l,L,P,N;
	double temp;
	SEXP res;
	SEXP Rdim;
	
	Rdim = getAttrib(X, R_DimSymbol);
  N = INTEGER(Rdim)[0];
  P = INTEGER(Rdim)[1];
	
	
	X = coerceVector(X,REALSXP);
	times_start = coerceVector(times_start,INTSXP);
	times_end = coerceVector(times_end,INTSXP);
	
	L = length(times_start);
	
	protect(res = allocMatrix(REALSXP,P,L));
	
	for (l = 0; l < L; l++) {
	  for (j = 0; j < P; j++){
	  	temp = 0;
	  	for (i = INTEGER(times_start)[l] - 1; i < INTEGER(times_end)[l]; i++){
	  		temp += REAL(X)[i + j * N];
	  	}
	  	
	  	REAL(res)[ j + l * P] = temp / (INTEGER(times_end)[l] - INTEGER(times_start)[l] + 1);
	  }
	}
	
	UNPROTECT(1);
	return res;	
}


SEXP log_norm(SEXP Rvals,SEXP Rmax, SEXP RID)
{
  
  int i,K,ID;
  double temp = 1;
  double log_temp;
  double temp_sum =0;
  double max;
  SEXP Rval;
  SEXP z_norm;
  SEXP res;
  SEXP log_Rval;

  
  
  Rvals = coerceVector(Rvals,REALSXP);
  /* RID = coerceVector(RID,INTSXP);*/
  K=length(Rvals);
  max = REAL(Rmax)[0];
  ID = INTEGER(RID)[0];

  protect(Rval=allocVector(REALSXP,K));
  protect(log_Rval=allocVector(REALSXP,K));
  protect(z_norm=allocVector(REALSXP,1));
  protect(res=allocVector(VECSXP,3));

  
  for (i = 0; i < K; i++) {
    if ( i != ID) {
      temp += exp(REAL(Rvals)[i] - max);
    }

  }
  temp = log(temp) + max;
  log_temp = temp;
  temp = exp(temp);
  
  for (i = 0; i < K; i++) {
    REAL(log_Rval)[i] = REAL(Rvals)[i] - log_temp;
    REAL(Rval)[i] = exp(REAL(log_Rval)[i]);

  }
  SET_VECTOR_ELT(res,0,Rval);
  REAL(z_norm)[0] = log_temp;
  SET_VECTOR_ELT(res,1,log_Rval);
  SET_VECTOR_ELT(res,2,z_norm);
  
  UNPROTECT(4);

  return res;  
}

SEXP eval_weights_risk(SEXP simu_exp, SEXP data){
  SEXP Rdim1;
  SEXP Rdim2;
  SEXP res;
  double temp_sum;
  double var_2;
  double argument;
  int i, j, k, l, N, P, K;
  
  Rdim1 = getAttrib(simu_exp, R_DimSymbol);
  Rdim2 = getAttrib(data, R_DimSymbol);
  N = INTEGER(Rdim2)[0];
  P = INTEGER(Rdim2)[1];
  K = INTEGER(Rdim1)[2];
  
  simu_exp = coerceVector(simu_exp,REALSXP);
  data = coerceVector(data,REALSXP);
  
  protect(res=allocVector(REALSXP, K ) );
  
  for (k = 0; k < K; k++){
    temp_sum = 0;
    for (j = 1; j < P; j++){
      for (i = 0; i < N; i++){
        var_2 = 0.01 + 0.04 * REAL(simu_exp)[i + N*j + P*N*k] * REAL(simu_exp)[i + N*j + P*N*k];
        argument = REAL(simu_exp)[i + N*j + P*N*k] - REAL(data)[i + N*j];
        temp_sum += (- argument * argument / var_2 - log(var_2) ) /2;
      }
    }
    REAL(res)[k] = temp_sum;
  }
  
  UNPROTECT(1);
  return res;
}

SEXP sum_over_time(SEXP simu_exp, SEXP data_arr){
  SEXP Rdim1;
  SEXP Rdim2;
  SEXP res;
  SEXP dim_res;
  double temp_sum;
  double var_2;
  double argument;
  int i, t, T, j, k, l, n, N, P, K;
  
  Rdim1 = getAttrib(simu_exp, R_DimSymbol);
  Rdim2 = getAttrib(data_arr, R_DimSymbol);
  T = INTEGER(Rdim1)[0];
  P = INTEGER(Rdim1)[1];
  K = INTEGER(Rdim1)[2];
  
  N = INTEGER(Rdim2)[2];
  
  simu_exp = coerceVector(simu_exp,REALSXP);
  data_arr = coerceVector(data_arr,REALSXP);
  
  protect(res=allocVector(REALSXP, (P-1) * K * N ) );
  PROTECT(dim_res = allocVector(INTSXP, 3));
  INTEGER(dim_res)[0] = (P-1);
  INTEGER(dim_res)[1] = K;
  INTEGER(dim_res)[2] = N;
  
  setAttrib(res, R_DimSymbol, dim_res);
  
  for (n = 0; n < N; n++){
  /*printf("%i\n", n);*/
    for (k = 0; k < K; k++){
      for (j = 1; j < P; j++){
        temp_sum = 0;
        for (t = 0; t < T; t++){
          var_2 = 0.01 + 0.04 * REAL(simu_exp)[t + T*j + P*T*k] * REAL(simu_exp)[t + T*j + P*T*k];
          argument = REAL(simu_exp)[t + T*j + P*T*k] - REAL(data_arr)[t + T*j + P*T*n];
          temp_sum += (- argument * argument / var_2 - log(var_2) ) /2;
        }
      	REAL(res)[j-1 + k * (P-1) + K * (P-1) * n] = temp_sum;
      }
    }
  }
  UNPROTECT(2);
  return res;
}

