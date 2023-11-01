#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

// Function to calculate Stein-Kernel Distance
NumericVector stein_thin_indices_cpp(NumericMatrix samp, int m) {
  int n = samp.nrow();
  int d = samp.ncol();
  double b = -0.5; // Values required for Stein-Kernel
  double c = 1;
  
  // Calculate KSD of the form K_p(X_i, x_i)
  NumericVector KSD_n(n);
  
  for (int i_n = 0; i_n < n; i_n++) {
    NumericVector x = samp(i_n, _);
    NumericVector y = samp(i_n, _);
    NumericVector w_sq(d);
    
    NumericVector b_x = -x;
    NumericVector b_y = -y;
    
    for (int j_kp = 0; j_kp < d; j_kp++) {
      double val = pow(c, 2) + sum(pow(x - y, 2));
      double term1 = b_x[j_kp] * b_y[j_kp] * pow(val, b);
      double term2 = b_x[j_kp] * (b * pow(val, b - 1) * (-2) * (x[j_kp] - y[j_kp]));
      double term3 = b_y[j_kp] * (b * pow(val, b - 1) * 2 * (x[j_kp] - y[j_kp]));
      double term4 = -2 * b * ((b - 1) * pow(val, b - 2) * 2 * pow(x[j_kp] - y[j_kp], 2) + pow(val, b - 1));
      
      w_sq[j_kp] = term1 + term2 + term3 + term4;
    }
    
    KSD_n[i_n] = sum(w_sq);
  }
  
  // We are intended to choose the indices
  NumericVector ind_chosen(m);
  ind_chosen[0] = which_min(KSD_n); // First value is based on KSD_n only
  
  // Loop for j = 2:m
  for (int j = 1; j < m; j++) {
    NumericVector KSD_trail_n(n); // Storing the required quantity for each n
    
    // Loop for n
    for (int i = 0; i < n; i++) {
      double f_term = KSD_n[i] / 2; // First term
      double s_term = 0;            // Second term
      NumericVector y = samp(i, _);
      
      for (int j_d = 0; j_d < j; j_d++) { // Loop for j'
        NumericVector x = samp(ind_chosen[j_d], _);
        NumericVector w_sq(d);
        NumericVector b_x = -x;
        NumericVector b_y = -y;
        
        for (int j_kp = 0; j_kp < d; j_kp++) {
          double val = pow(c, 2) + sum(pow(x - y, 2));
          double term1 = b_x[j_kp] * b_y[j_kp] * pow(val, b);
          double term2 = b_x[j_kp] * (b * pow(val, b - 1) * (-2) * (x[j_kp] - y[j_kp]));
          double term3 = b_y[j_kp] * (b * pow(val, b - 1) * 2 * (x[j_kp] - y[j_kp]));
          double term4 = -2 * b * ((b - 1) * pow(val, b - 2) * 2 * pow(x[j_kp] - y[j_kp], 2) + pow(val, b - 1));
          
          w_sq[j_kp] = term1 + term2 + term3 + term4;
        }
        
        s_term += sum(w_sq);
      }
      
      KSD_trail_n[i] = f_term + s_term;
    }
    
    ind_chosen[j] = which_min(KSD_trail_n);
  }
  
  return ind_chosen + 1;
}

