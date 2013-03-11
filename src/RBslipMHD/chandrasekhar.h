#include<stdlib.h>
#include<math.h>
#define _USE_MATH_DEFINES

double xfunc_stationary(double x, double Q, 
			double P1, double P2){
  // Note: P1 & P2 not used here. Value passed does not matter.
  return 2*pow(x, 3) + 3*pow(x, 2) - 1 - Q/pow(M_PI, 2);
}

double xfunc_oscillatory(double x, double Q, double P1, double P2){
  return 2*pow(x, 3) + 3*pow(x, 2) - 1 - 
    Q*pow(P1, 2)/((pow(M_PI, 2)*(1 + P1)*(P1 + P2)));
}

double solvex(double Q, double P1, double P2,
	      double left, double right, 
	      double (*xfunc)(double, double, double, double)){
  double middle = (left + right) / 2,
    fl = (*xfunc)(left, Q, P1, P2),
    fr = (*xfunc)(right, Q, P1, P2),
    fm = (*xfunc)(middle, Q, P1, P2);
  if(fr - fl < 0.0001){
    return (left + right)/2.0;
  }
  else if(fm > 0)
    return solvex(Q, P1, P2, left, middle, xfunc);
  else 
    return solvex(Q, P1, P2, middle, right, xfunc);
}

double findRcVert_stationary(double Q){
  /* Here P1 & P2 have no meaning. Used for avoiding type mismatch while
  passing fn pointer.

  At x=0, xfunc_stationary is -ve and when x = 1 + Q/ M_PI^2 xfunc_stationary
  is +ve as for x>1, 2x^3 + 3x^2 - x > 0.

  So these bounds are passed during bissection method call
  */
  double x = solvex(Q, 0, 0, 0, (1 + Q/( M_PI * M_PI)), xfunc_stationary);
  return pow(M_PI, 4) * (1 + x) * (pow(1 +x, 2) + Q / pow(M_PI, 2)) / x;
}

double findKVert_stationary(double Q){
  /* Here P1 & P2 have no meaning. Used for avoiding type mismatch while
  passing fn pointer.
  
  Bounds logic same as above.
  */
  double x = solvex(Q, 0, 0, 0, (1 + Q/( M_PI * M_PI)), xfunc_stationary);
  return M_PI * sqrt(x);
}

double findRcVert_oscillatory(double Q, double P1, double P2){
  /*
    Bounds logic similar as above.
   */
  double r = 1 + Q*pow(P1, 2)/(pow(M_PI, 2)*(1 + P1)*(P1 + P2));
  double x = solvex(Q, P1, P2, 0, r, xfunc_oscillatory);
  
  return (((pow(M_PI, 2))*(1 + P2)*(P1 + P2)*(1 + x))/(pow(P2, 2)*x))*
    (pow(M_PI, 2)*pow((1 + x), 2) + Q*pow(P1, 2)/((1 + P1)*(P1 + P2)));
}

double findKVert_oscillatory(double Q, double P1, double P2){
  /*
    Bounds logic similar as above.
  */
  double r = 1 + Q*pow(P1, 2)/(pow(M_PI, 2)*(1 + P1)*(P1 + P2));
  double x = solvex(Q, P1, P2, 0, r, xfunc_oscillatory);
  
  return M_PI * sqrt(x);
}
