#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppEigen)]]

using namespace Rcpp;
using namespace arma;

static int num(IntegerVector x, int c);
static arma::uvec which(IntegerVector x, int c);

// [[Rcpp::export]]
Rcpp::List imoscato(arma::mat X, arma::mat Y, IntegerVector cell_type, IntegerVector z, NumericVector v,  NumericVector s, arma::mat P, bool sc_nb, bool st_nb, bool find_domain, int iter, int burn) {
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;
  int m = X.n_rows;
  int K = Rcpp::unique(cell_type).size();
  int D = Rcpp::unique(z).size();
  
  NumericVector number(K);
  NumericVector prop(K);
  for (int k = 0; k < K; k++){
    number(k) = num(cell_type, k);
    prop(k) = number(k)/m;
    // prop(k) = pow(K, -1);
  }
  
  // Markov random field (MRF) prior setting
  IntegerVector e(D);
  for (int d = 0; d < D; d++){
    e(d) = 1;
  }
  int n_neighbor = P.n_cols;
  int min_element = 10;
  double f = 1;

  // Calculate the minimum and maximum sample normalized expression levels used for add-delete algorithm
  arma::mat Mu(m, p);
  for (int i = 0; i < m; i++){
    for (int j = 0; j < p; j++){
      Mu(i, j) = X(i, j)/v(i);
    }
  }

  // Hyperparamters
  double delta_p = 0.5;
  double a_alpha = 0.001;
  double b_alpha = 0.001;
  double a_mu = 0.001;
  double b_mu = 0.001;
  double a_0 = 0.001;
  double b_0 = 0.001;
  double a_phi = 0.001;
  double b_phi = 0.001;
  double a_psi = 0.001;
  double b_psi = 0.001;
  double theta = 0.5;
  double vartheta = 0.5;
  double omega = 0.05;

  // Tuning parameters
  double tau_omega = 0.1;
  double tau_pi = 0.1;
  double tau_alpha = 0.1;
  double tau_mu = 0.1;
  double tau_0 = 0.1;
  double tau_phi = 1;
  double tau_psi = 1;
  int p_update = p*0.2;
  
  // Set temporary variables
  bool zidd = true;
  bool zinb = true;
  bool st_data = true;
  bool feature_selection = true;
  int t, i, ii, j, jj, k, kk, ee, d, num = 10, gamma_new, delta_new, Delta_temp, p_temp, H_temp, Z_temp, neighbor_index;
  double hastings, accept_gamma = 0, accept_omega = 0, accept_delta = 0, accept_alpha = 0, accept_mu = 0, accept_phi = 0, accept_psi = 0;
  double try_gamma = 0, try_omega = 0, try_delta = 0, try_alpha = 0, try_mu = 0, try_phi = 0, try_psi = 0;
  double d_temp, h_temp, mu_sum_temp, mu_sum_new, omega_temp, omega_new, alpha_temp, alpha_new, phi_new, phi_temp, psi_new, psi_temp, map_z_temp, map_z, map_gamma_temp;
  double pi = 3.141593, max_temp, sum_temp;
  double c_mu;
  IntegerVector Delta_sum_temp(n, 0);
  IntegerVector z_map, gamma_map;
  NumericMatrix A_new(D, K);
  NumericMatrix M_new(K, p);
  IntegerVector H_sum_temp(n, 0);
  IntegerVector Z_sum_temp(m, 0);
  NumericVector pi_temp(K);
  NumericVector pi_new(K);
  NumericVector Omega_new(K);
  
  
  // Create the spaces to store the results
  arma::cube Pi_store(n, K, iter);
  arma::cube Omega_store(n, K, iter);
  NumericMatrix Delta_ppi(n, K);
  IntegerVector Delta_sum(iter);
  arma::cube A_store(D, K, iter);
  IntegerMatrix z_store(iter, n);
  NumericVector map_z_store(iter);
  IntegerMatrix gamma_store(iter, p);
  IntegerVector gamma_sum(iter);
  NumericVector gamma_ppi(p);
  NumericVector gamma_BF(p);
  NumericVector map_gamma_store(iter);
  arma::cube M_store(K, p, iter);
  // arma::cube H_store(n, p, iter);
  NumericMatrix H_ppi(n, p);
  IntegerVector H_sum(iter);
  // arma::cube Z_store(m, p, iter);
  NumericMatrix Z_ppi(m, p);
  IntegerVector Z_sum(iter);
  NumericMatrix phi_store(iter, p);
  NumericMatrix psi_store(iter, p);

  // Initialization
  NumericMatrix A(D, K);
  for (k = 0; k < K; k++){
    for (d = 0; d < D; d++){
      A(d, k) = 1;
    }
  }

  NumericMatrix ct_p(n, K);
  NumericMatrix Omega(n, K);
  IntegerMatrix Delta(n, K);
  
  if (zidd) {
    for (i = 0; i < n; i++){
      for (k = 0; k < K; k++){
        Omega(i, k) = 1;
        
      }
    }
    
    for (i = 0; i < n; i++){
      for (k = 0; k < K; k++){
        ct_p(i, k) = Omega(i, k)/sum(Omega(i, _));
      }
    }
    
    for (i = 0; i < n; i++) {
      for (k = 0; k < K; k++) {
        Delta_ppi(i, k) = 0;
        
        if (Omega(i, k) != 0) {
          Delta(i, k) = 0;
        } else {
          Delta(i, k) = 0;
          Delta_sum_temp(i) = Delta_sum_temp(i) + Delta(i, k);
        }
      }
    }
    
  }else{
    for (k = 0; k < K; k++){
      for (i = 0; i < n; i++){
        ct_p(i, k) = prop(k);
        // ct_p(i, k) = pow(K, -1);
        
      }
    }
  }
  
  IntegerVector gamma(p);
  NumericMatrix M(K, p);
  for (k = 0; k < K; k++){
    for (j = 0; j < p; j++){
      M(k, j) = mean(mean(X));
    }
  }

  IntegerMatrix H(n, p);
  IntegerMatrix Z(m, p);
  double phi_start;
  double psi_start;
  
  if (sc_nb){
    psi_start = 10;
  }else{
    psi_start = 100000000;
  }
  if (st_nb){
    phi_start = 10;
  }else{
    phi_start = 100000000;
  }
  NumericVector phi(p);
  NumericVector psi(p);
  
  for (j = 0; j < p; j++) {
    phi(j) = phi_start;
    psi(j) = psi_start;

    for (i = 0; i < n; i++) {
      H_ppi(i, j) = 0;
      if (Y(i, j) != 0) {
        H(i, j) = 0;
      } else {
        H(i, j) = 0;
        H_sum_temp(i) = H_sum_temp(i) + H(i, j);
      }
    }

    for (ii = 0; ii < m; ii++){
      Z_ppi(ii, j) = 0;

      if (X(ii, j) != 0) {
        Z(ii, j) = 0;
      } else {
        Z(ii, j) = 0;
        Z_sum_temp(ii) = Z_sum_temp(ii) + Z(ii, j);

      }
    }
  }

  IntegerVector clusters(K-1);
  for (k = 0; k < (K-1); k++){
    clusters(k) = k;
  }
  
  // Spatial domains
  NumericVector prob_init(D);
  IntegerVector domains(D);

  for (d = 0; d < D; d++){
    domains(d) = d;
    prob_init(d) = pow(D, -1);
  }

  if (find_domain) {
    IntegerVector z(n);
    for (int i = 0; i < n; i++){
      z(i) = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(domains, 1, TRUE, prob_init));
    }
  }
  

  // Start MCMC algorithms
  for (t = 0; t < iter; t++){
    if (zidd){
      // Add-delete algorithm to update cell type proportion Pi, Omega, and Delta
      for (i = 0; i < n; i++) {
        for (ee = 0; ee < 5; ee++) {
          k = rand()%K;
          delta_new = 1 - Delta(i, k);

          for (kk = 0; kk < K; kk++){
            pi_temp(kk) = ct_p(i, kk);
            pi_new(kk) = ct_p(i, kk);
            Omega_new(kk) = Omega(i, kk);
          }
          
          omega_temp = Omega(i, k);

          // Delete step
          if (delta_new == 0){
            omega_new = rgamma(1, A(z(i), k), 1 )(0);
            Omega_new(k) = omega_new;
            
            for (kk = 0; kk < K; kk++){
              pi_new(kk) = Omega_new(kk)/sum(Omega_new);
            }
            
            // Calculate hastings ratio
            hastings = 0;

            // SRT data likelihood ratio
            for (j = 0; j < p; j++){
              if (gamma(j) == 1){
                mu_sum_new = sum(pi_new*M(_, j));
                mu_sum_temp = sum(pi_temp*M(_, j));

                if (H(i, j) == 0){
                  hastings = hastings + phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
                  hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
                }
              }
            }

            // // Conditional prior of omega given delta
            // hastings = hastings + (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_new) - omega_new);
            // 
            // // Proposal density ratio
            // hastings = hastings - (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_new) - omega_new);

            // Prior ratio of delta
            hastings = hastings + log(1 - delta_p);
            hastings = hastings - log(delta_p);

            // Check if accept the proposed new values
            if (t >= burn) {
              try_delta = try_delta + 1;
            }

            if (hastings >= log(double(rand()%10001)/10000)) {
              Delta(i, k) = delta_new;
              Omega(i, k) = omega_new;
              for (kk = 0; kk < K; kk++){
                ct_p(i, kk) = Omega(i, kk)/sum(Omega(i, _));
              }
              
              if (t >= burn) {
                accept_delta = accept_delta + 1;
              }
            }

          } // End of delete step

          // Add step
          else {
            omega_new = 0;
            Omega_new(k) = omega_new;
            for (kk = 0; kk < K; kk++){
              pi_new(kk) = Omega_new(kk)/sum(Omega_new);
            }

            // Calculate hastings ratio
            hastings = 0;

            // SRT data likelihood ratio
            for (j = 0; j < p; j++){
              if (gamma(j) == 1){
                mu_sum_new = sum(pi_new*M(_, j));
                mu_sum_temp = sum(pi_temp*M(_, j));

                if (H(i, j) == 0){
                  hastings = hastings + phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
                  hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
                }
              }
            }

            // // Conditional prior of omega given delta
            // hastings = hastings - (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_temp) - omega_temp);
            // 
            // // Proposal density ratio
            // hastings = hastings + (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_temp) - omega_temp);

            // Prior ratio of delta
            hastings = hastings + log(delta_p);
            hastings = hastings - log(1 - delta_p);

            // Check if accept the proposed new values
            if (t >= burn) {
              try_delta = try_delta + 1;
            }

            if (hastings >= log(double(rand()%10001)/10000)) {
              // For any i, delta_{ik} cannot be 1 for all k
              if (sum(Delta(i, _)) < K - 1){
                Delta(i, k) = delta_new;
                Omega(i, k) = omega_new;
                for (kk = 0; kk < K; kk++){
                  ct_p(i, kk) = Omega(i, kk)/sum(Omega(i, _));
                }
                
                if (t >= burn) {
                  accept_delta = accept_delta + 1;
                }
              }
              
            }

          } // End of add step

        }
      } // End of updating cell type proportion Pi, Omega, and Delta
      
      
      // RWMH to update the cell type proportion Pi and Omega ##################
      for (i = 0; i < n; i++) {
        for (k = 0; k < K; k++) {
          if (Delta(i, k) == 0){
            
            for (kk = 0; kk < K; kk++){
              pi_temp(kk) = ct_p(i, kk);
              pi_new(kk) = ct_p(i, kk);
              Omega_new(kk) = Omega(i, kk);
            }
            
            omega_temp = Omega(i, k);
            // omega_new = exp(rnorm(1, log(omega_temp), tau_omega)(0));
            omega_new = r_truncnorm(omega_new, 10*tau_omega, 1, 100);

            Omega_new(k) = omega_new;
            
            for (kk = 0; kk < K; kk++){
              pi_new(kk) = Omega_new(kk)/sum(Omega_new);
            }
            
            hastings =  0;
            
            // SRT data likelihood ratio
            for (j = 0; j < p; j++){
              if (gamma(j) == 1){
                mu_sum_new = sum(pi_new*M(_, j));
                mu_sum_temp = sum(pi_temp*M(_, j));
                
                if (H(i, j) == 0){
                  hastings = hastings + phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
                  hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
                }
              }
            }
            
            // Prior ratio
            hastings = hastings + (A(z(i), k) - 1)*log(omega_new) - omega_new;
            hastings = hastings - ((A(z(i), k) - 1)*log(omega_temp) - omega_temp);
            
            // Check if accept the proposed new values
            if (hastings >= log(double(rand()%10001)/10000)) {
              Omega(i, k) = omega_new;
              
              for (kk = 0; kk < K; kk++){
                ct_p(i, kk) = Omega(i, kk)/sum(Omega(i, _));
              }
              
              if (t >= burn){
                accept_omega = accept_omega + 1;
              }
            }
            
            if (t >= burn){
              try_omega = try_omega + 1;
            }
          }
          
        }
      } // End of updating the cell type proportion Pi and Omega ###############
      
    } else {
      // RWMH to update the cell type proportion Pi ##############################
      for (i = 0; i < n; i++) {
        // Propose new pi vector
        for (kk = 0; kk < K; kk++){
          pi_temp(kk) = ct_p(i, kk);
          pi_new(kk) = ct_p(i, kk);
        }
        
        int k_temp = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE));
        h_temp = r_truncnorm(0, tau_pi, -pi_temp(k_temp + 1), pi_temp(k_temp));
        
        pi_new(k_temp) = pi_temp(k_temp) - h_temp;
        pi_new(k_temp + 1) = pi_temp(k_temp + 1) + h_temp;
        
        // Calculate hastings ratio
        hastings = 0;
        
        // SRT data likelihood ratio
        for (j = 0; j < p; j++){
          if (gamma(j) == 1){
            mu_sum_new = sum(pi_new*M(_, j));
            mu_sum_temp = sum(pi_temp*M(_, j));
            if (H(i, j) == 0){
              hastings = hastings + phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
              hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
            }
          }
        }
        
        // Prior ratio
        for (k = 0; k < K; k++) {
          hastings = hastings + (A(z(i), k) - 1)*log(pi_new(k));
          hastings = hastings - (A(z(i), k) - 1)*log(pi_temp(k));
        }
        
        // check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)){
          for (k = 0; k < K; k++){
            ct_p(i, k) = pi_new(k);
          }
        }
      } // End of updating the cell type proportion Pi #########################
    }
    
    
    // RWMH to update the concentration parameter A  ###########################
    for (d = 0; d < D; d++){
      for (k = 0; k < K; k++){
        alpha_temp = A(d, k);
        
        for (kk = 0; kk < K; kk++){
          A_new(d, kk) = A(d, kk);
        }
        
        // alpha_new = exp(rnorm(1, log(alpha_temp), tau_alpha)(0));
        alpha_new = r_truncnorm(alpha_temp, 10*tau_alpha, 1, 100);
        A_new(d, k) = alpha_new;
        
        hastings = 0;
        
        if (zidd){ // zero-inflated gamma distribution
          // Conditional prior ratio
          for (i = 0; i < n; i++){
            if (z(i) == d && Delta(i, k) == 0){
              hastings = hastings + (-lgamma(alpha_new) + (alpha_new - 1)*log(Omega(i, k)));
              hastings = hastings - (-lgamma(alpha_temp) + (alpha_temp - 1)*log(Omega(i, k)));
            }
          }
          
        } else {
          // Conditional prior ratio
          for (i = 0; i < n; i++){
            if (z(i) == d){
              hastings = hastings + lgamma(sum(A_new(d, _)));
              hastings = hastings - lgamma(sum(A(d, _)));
              
              for (kk = 0; kk < K; kk++) {
                hastings = hastings + (-lgamma(A_new(d, kk)) + (A_new(d, kk) - 1)*log(ct_p(i, kk)));
                hastings = hastings - (- lgamma(A(d, kk)) + (A(d, kk) - 1)*log(ct_p(i, kk)));
              }
              
            }
          }
          
        }
        
        // Prior ratio
        hastings = hastings + (a_alpha - 1)*log(alpha_new) - b_alpha*alpha_new;
        hastings = hastings - ((a_alpha - 1)*log(alpha_temp) - b_alpha*alpha_temp);
        
        // Check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)) {
          A(d, k) = alpha_new;
          
          if (t >= burn){
            accept_alpha = accept_alpha + 1;
          }
        }
        
        if (t >= burn){
          try_alpha = try_alpha + 1;
        }
      }
    } // End of updating concentration parameter A #############################
    

    // Gibbs Sampling to update the spatial domain allocation z ################
    if (find_domain) {
      for (i = 0; i < n; i++) {
        
        // Compute I(z_i = z_i')
        IntegerVector neighbor_temp(D);
        
        for (int ii = 0; ii < n_neighbor; ii++){
          if (P(i, ii) != 0) {
            neighbor_index = P(i, ii) - 1;
            neighbor_temp(z(neighbor_index)) = neighbor_temp(z(neighbor_index)) + 1;
          }
        }
        
        // Calculate posterior probability of z
        NumericVector loglklh(D);
        for (d = 0; d < D; d++){
          if (zidd) { // zero-inflated gamma distribution
            // Conditional prior ratio of omega
            for (k = 0; k < K; k++){
              if (Delta(i, k) == 0) {
                loglklh(d) = loglklh(d) -lgamma(A(d, k)) + (A(d, k) - 1)*log(Omega(i, k)) - Omega(i, k);
                
              }
              
            }
            
          } else{
            // Conditional prior ratio of pi
            loglklh(d) = loglklh(d) + lgamma(sum(A(d, _)));
            
            for (k = 0; k < K; k++) {
              loglklh(d) = loglklh(d) - lgamma(A(d, k)) + (A(d, k) - 1)*log(ct_p(i, k));
            }
            
          }
          
          // Prior ratio of z
          loglklh(d) = loglklh(d) + e(d) +  f*(neighbor_temp(d));
          
        }
        
        
        // Posterior probability of each spatial domain
        NumericVector prob(D);
        for (d = 0; d < D; d++) {
          double diff = 0;
          for (int m = 0; m < D; m++) {
            if (m != d) {
              diff = diff + exp(loglklh(m) - loglklh(d));
            }
          }
          prob(d) = 1/(1 + diff);
        }
        
        
        // Generate new sample
        int z_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(domains, 1, TRUE, prob));
        
        // Make sure each group has at least multiple observations
        d = z(i);
        if (sum(z == d) == min_element && z_new != d) {
          z(i) = d;
        }else {
          z(i) = z_new;
        }
      } 
      
    } // End of updating spatial domain allocation parameter z ###############


    // Add-delete algorithm to update feature selection parameter gamma ########
    if (feature_selection) {
      for (ee = 0; ee < p_update; ee++) {
        j = rand()%p;
        gamma_new = 1 - gamma(j);
        // double c_mu = mean(M(_, j));
        double c_mu = mean(Mu.col(j));
        
        // Delete step
        if (gamma_new == 0){
          hastings = 0;
          
          NumericVector mu_temp = M(_, j);
          double mu_new = exp(rnorm(1, log(c_mu), tau_0)(0));
          
          // SRT data likelihood ratio
          if (st_data){
            for (i = 0; i < n; i++) {
              mu_sum_temp = sum(ct_p(i, _)*M(_, j));
              mu_sum_new = mu_new;
              if (H(i, j) == 0){
                hastings = hastings + (phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j))));
                hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
              }
            }
            
          }
          
          // scRNA-seq data likelihood ratio
          for (i = 0; i < m; i++) {
            if (Z(i, j) == 0) {
              hastings = hastings + psi(j)*( - log(mu_new*v(i) + psi(j))) + X(i, j)*(log(mu_new*v(i)) - log(mu_new*v(i) + psi(j)));
              hastings = hastings - (psi(j)*( - log(mu_temp(cell_type(i))*v(i) + psi(j))) + X(i, j)*(log(mu_temp(cell_type(i))*v(i)) - log(mu_temp(cell_type(i))*v(i) + psi(j))));
            }
          }
          
          // Conditional prior ratio
          hastings = hastings + (a_0*log(b_0) - lgamma(a_0) + (a_0 - 1)*log(mu_new) - b_0*mu_new);
          for (k = 0; k < K; k++) {
            hastings = hastings - (a_mu*log(b_mu) - lgamma(a_mu) + (a_mu - 1)*log(mu_temp(k)) - b_mu*mu_temp(k));
          }
          
          // Prior ratio
          hastings = hastings - log(omega) + log(1 - omega);
          
          // Conditional proposal ratio
          for (k = 0; k < K; k++) {
            hastings = hastings + (-log(2*pi)/2 - log(10*tau_mu) - pow((log(mu_temp(k)) - log(c_mu)), 2)/2/pow(10*tau_mu, 2));
          }
          hastings = hastings - (-log(2*pi)/2 - log(tau_0) - pow((log(mu_new) - log(c_mu)), 2)/2/pow(tau_0, 2));
          
          // Check if accept the proposed new values
          if (t >= burn) {
            try_gamma = try_gamma + 1;
          }
          if (hastings >= log(double(rand()%10001)/10000)) {
            gamma(j) = gamma_new;
            for (k = 0; k < K; k++) {
              M(k, j) = mu_new;
            }
            
            if (t >= burn) {
              accept_gamma = accept_gamma + 1;
            }
          }
        } // End of delete step
        
        // Add step
        else{
          hastings = 0;
          double mu_temp = M(0, j);
          NumericVector mu_new(K);
          for (k = 0; k < K; k++){
            mu_new(k) = exp(rnorm(1, log(c_mu), 10*tau_mu)(0));
          }
          
          //SRT data likelihood ratio
          if (st_data) {
            for (i = 0; i < n; i++) {
              mu_sum_temp = mu_temp;
              mu_sum_new =  sum(ct_p(i, _)*mu_new);
              if (H(i, j) == 0){
                hastings = hastings + (phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j))));
                hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
              }
            }
            
          }
          
          
          // scRNA-seq data likelihood ratio
          for (i = 0; i < m; i++) {
            if (Z(i, j) == 0) {
              hastings = hastings + psi(j)*(- log(mu_new(cell_type(i))*v(i) + psi(j))) + X(i, j)*(log(mu_new(cell_type(i))*v(i)) - log(mu_new(cell_type(i))*v(i) + psi(j)));
              hastings = hastings - (psi(j)*(- log(mu_temp*v(i) + psi(j))) + X(i, j)*(log(mu_temp*v(i)) - log(mu_temp*v(i) + psi(j))));
            }
          }
          
          // Conditional prior ratio
          for (k = 0; k < K; k++) {
            hastings = hastings + (a_mu*log(b_mu) - lgamma(a_mu) + (a_mu - 1)*log(mu_new(k)) - b_mu*mu_new(k));
          }
          hastings = hastings - (a_0*log(b_0) - lgamma(a_0) + (a_0 - 1)*log(mu_temp) - b_0*mu_temp);
          
          // Prior ratio
          hastings = hastings + log(omega) - log(1 - omega);
          
          // Conditional proposal ratio
          hastings = hastings + (-log(2*pi)/2 - log(tau_0) - pow((log(mu_temp) - log(c_mu)), 2)/2/pow((tau_0), 2));
          for (k = 0; k < K; k++) {
            hastings = hastings - (-log(2*pi)/2 - log(10*tau_mu) - pow((log(mu_new(k)) - log(c_mu)), 2)/2/pow((10*tau_mu), 2));
          }
          
          // Check if accept the proposed new values
          if (t >= burn) {
            try_gamma = try_gamma + 1;
          }
          if (hastings >= log(double(rand()%10001)/10000)) {
            gamma(j) = gamma_new;
            for (k = 0; k < K; k++) {
              M(k, j) = mu_new(k);
            }
            if (t >= burn) {
              accept_gamma = accept_gamma + 1;
            }
          }
        } // End of add step
      } // End of updating gamma
    }

    
    for (j = 0; j < p; j++){
      gamma_store(t, j) = gamma(j);
    }
    gamma_sum(t) = sum(gamma);

    // Calculate PPI of gamma
    if (t >= burn) {
      for (j = 0; j < p; j++){
        gamma_ppi(j) = gamma_ppi(j) + gamma(j);
      }
    } // End of feature selection


    // RWMH to update the normalized expression levels M #######################
    for (j = 0; j < p; j++){
      if (gamma(j) == 0) {
        double mu_temp = M(0, j);
        double mu_new = exp(rnorm(1, log(mu_temp),  tau_mu)(0));

        hastings = 0;

        // SRT data likelihood ratio
        if (st_data) {
          for (i = 0; i < n; i++) {
            mu_sum_temp = mu_temp;
            mu_sum_new =  mu_new;
            if (H(i, j) == 0){
              hastings = hastings + (phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j))));
              hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
            }
          }
        }

        //scRNA-seq data likelihood ratio
        for (i = 0; i < m; i++) {
          if (Z(i, j) == 0) {
            hastings = hastings + psi(j)*(log(psi(j)) - log(mu_new*v(i) + psi(j))) + X(i, j)*(log(mu_new*v(i)) - log(mu_new*v(i) + psi(j)));
            hastings = hastings - (psi(j)*(log(psi(j)) - log(mu_temp*v(i) + psi(j))) + X(i, j)*(log(mu_temp*v(i)) - log(mu_temp*v(i) + psi(j))));

          }
        }

        // Prior ratio
        hastings = hastings + (a_0 - 1)*log(mu_new) - b_0*mu_new;
        hastings = hastings - ((a_0 - 1)*log(mu_temp) - b_0*mu_temp);

        // check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)) {
          for (k = 0; k < K; k++) {
            M(k, j) = mu_new;
          }

          if (t >= burn){
            accept_mu = accept_mu + 1;
          }
        }

        if (t >= burn){
          try_mu = try_mu + 1;
        }
      } // End of updating M for non-discriminating genes

    else{
      for (k = 0; k < K; k++) {
        double mu_temp = M(k, j);
        double mu_new = exp(rnorm(1, log(mu_temp), tau_mu)(0));

        M_new(_, j) = M(_, j);
        M_new(k, j) = mu_new;

        hastings = 0;

        // SRT data likelihood ratio
        if (st_data) {
          for (i = 0; i < n; i++) {
            mu_sum_temp =  sum(ct_p(i, _)*M(_, j));
            mu_sum_new =  sum(ct_p(i, _)*M_new(_, j));
            if (H(i, j) == 0){
              hastings = hastings + (phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j))));
              hastings = hastings - (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
            }
          }

        }


        // scRNA-seq data likelihood ratio
        for (i = 0; i < m; i++) {
          if (Z(i, j) == 0 && cell_type(i) == k) {
            hastings = hastings + psi(j)*(log(psi(j)) - log(mu_new*v(i) + psi(j))) + X(i, j)*(log(mu_new*v(i)) - log(mu_new*v(i) + psi(j)));
            hastings = hastings - (psi(j)*(log(psi(j)) - log(mu_temp*v(i) + psi(j))) + X(i, j)*(log(mu_temp*v(i)) - log(mu_temp*v(i) + psi(j))));
          }
        }

        // Prior ratio
        hastings = hastings + (a_mu - 1)*log(mu_new) - b_mu*mu_new;
        hastings = hastings - ((a_mu - 1)*log(mu_temp) - b_mu*mu_temp);

        // Check if accept the proposed new values
        if (hastings >= log(double(rand()%10001)/10000)) {
          M(k, j) = mu_new;

          if (t >= burn){
            accept_mu = accept_mu + 1;
          }
        }

        if (t >= burn){
          try_mu = try_mu + 1;
        }
      }

    } // End of updating M for discriminating genes
  } // End of updating the normalized expression levels M ######################


    // Gibbs sampler to update the zero-inflation parameter Z in scRNA-seq data
    if (zinb) {
      for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
          NumericVector prob_temp(2);
          double mu_temp =  M(cell_type(i), j);
          
          // Only update for x_{ij} = 0
          if (X(i, j) == 0) {
            prob_temp(0) = psi(j)*(log(psi(j)) - log(v(i)*mu_temp + psi(j))) + log(1 - vartheta);
            prob_temp(1) = log(vartheta);
            prob_temp(1) = 1/(1 + exp(prob_temp(0) - prob_temp(1)));
            Z_temp = Z(i, j);
            Z(i, j) = rbinom(1, 1, prob_temp(1))(0);
            Z_sum_temp(i) = Z_sum_temp(i) - Z_temp + Z(i, j);
          }
        }
      }
    }  // End of update the zero-inflation parameter Z in scRNA-seq data
    

    // Gipps sampler to update the zero-inflation parameter H in SRT data ###
    if (zinb) {
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          NumericVector prob_temp(2);
          mu_sum_temp =  sum(ct_p(i, _)*M(_, j));

          // only update for y_{ij} = 0
          if (Y(i, j) == 0) {
            prob_temp(0) = phi(j)*(log(phi(j)) - log(s(i)*mu_sum_temp + phi(j))) + log(1 - theta);
            prob_temp(1) = log(theta);
            prob_temp(1) = 1/(1 + exp(prob_temp(0) - prob_temp(1)));
            H_temp = H(i, j);
            H(i, j) = rbinom(1, 1, prob_temp(1))(0);
            H_sum_temp(i) = H_sum_temp(i) - H_temp + H(i, j);
          }
        }
      }
    }  // End of update of the zero-inflation parameter H in SRT data
    

    // RWMH to update of the dispersion parameter psi in scRNA-seq data ########
    if(sc_nb) {
      for (j = 0; j < p; j++) {
        psi_temp = psi(j);
        psi_new = exp(r_truncnorm(log(psi(j)), tau_psi, log(1), log(100)));
        // psi_new = exp(rnorm(1, log(psi_temp), tau_psi)(0));
        
        hastings = 0;
        
        for (i = 0; i < m; i++) {
          double mu_temp =  M(cell_type(i), j);
          
          // SRT data likelihood ratio
          if (Z(i, j) == 0){
            hastings = hastings + lgamma(X(i, j) + psi_new) - lgamma(psi_new) +
              psi_new*(log(psi_new) - log(mu_temp*v(i) + psi_new)) + X(i, j)*( - log(mu_temp*v(i) + psi_new));
            hastings = hastings - (lgamma(X(i, j) + psi_temp) - lgamma(psi_temp) +
              psi_temp*(log(psi_temp) - log(mu_temp*v(i) + psi_temp)) + X(i, j)*( - log(mu_temp*v(i) + psi_temp)));
          }
        }
        
        // Prior ratio
        hastings = hastings + (a_psi - 1)*log(psi_new) - b_psi*psi_new;
        hastings = hastings - ((a_psi - 1)*log(psi_temp) - b_psi*psi_temp);
        
        // Check if accept the proposed new value
        if (t >= burn){
          try_psi = try_psi + 1;
        }
        if (hastings >= log(double(rand()%10001)/10000)){
          psi(j) = psi_new;
          
          if (t >= burn){
            accept_psi = accept_psi + 1;
          }
        }
        
      }
    } // End of updating the dispersion parameter psi in scRNA-seq data
    
    
    // RWMH to update of the dispersion parameter Phi in SRT data ############
    if (st_nb) {
      for (j = 0; j < p; j++) {
        phi_temp = phi(j);
        phi_new = exp(r_truncnorm(log(phi(j)), tau_phi, log(1), log(100)));
        // phi_new = exp(rnorm(1, log(phi_temp),  tau_phi)(0));
        hastings = 0;

        for (i = 0; i < n; i++) {
          mu_sum_temp =  sum(ct_p(i, _)*M(_, j));

          // SRT data likelihood ratio
          if (H(i, j) == 0){
            hastings = hastings + lgamma(Y(i, j) + phi_new) - lgamma(phi_new) +
              phi_new*(log(phi_new) - log(mu_sum_temp*s(i) + phi_new)) + Y(i, j)*( - log(mu_sum_temp*s(i) + phi_new));
            hastings = hastings - (lgamma(Y(i, j) + phi_temp) - lgamma(phi_temp) +
              phi_temp*(log(phi_temp) - log(mu_sum_temp*s(i) + phi_temp)) + Y(i, j)*( - log(mu_sum_temp*s(i) + phi_temp)));
          }
        }

        // Prior ratio
        hastings = hastings + (a_phi - 1)*log(phi_new) - b_phi*phi_new;
        hastings = hastings - ((a_phi - 1)*log(phi_temp) - b_phi*phi_temp);

        // Check if accept the proposed new value
        if (t >= burn){
          try_phi = try_phi + 1;
        }
        if (hastings >= log(double(rand()%10001)/10000)){
          phi(j) = phi_new;

          if (t >= burn){
            accept_phi = accept_phi + 1;
          }
        }

      }
    } // End of updating the dispersion parameter Phi in SRT data


    // Monitor the process
    if((t*100/(iter-1)) == num) {
      Rcout<<num<< "% has been done\n";
      num = num + 10;
    }

    // Calculate the marginal posterior probability to obtain MAP ############
    double map_z = 0;
    
    for (i = 0; i < n; i++) {
      z_store(t, i) = z(i);
      
      // Prior of z
      int n_neighbor_temp = 0;
      
      // Compute I(z_i = z_i')
      for (int ii = 0; ii < n_neighbor; ii++){
        if (P(i, ii) != 0) {
          if(z(i) == z(P(i, ii) - 1)){
            n_neighbor_temp = n_neighbor_temp + 1;
          }
        }
      }
      
      map_z = map_z + e(z(i)) + f*n_neighbor_temp;
      
      if (zidd) {
        for (k = 0; k < K; k++) {
          map_z = map_z - lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(Omega(i, k)) - Omega(i, k);
        }
        
      } else {
        // Conditional prior of pi
        map_z = map_z + lgamma(sum(A(z(i), _)));
        
        for (k = 0; k < K; k++) {
          map_z = map_z - lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(ct_p(i, k));
        }
      }
      
    }
    
    // Calculate the marginal posterior probability of gamma to obtain MAP =====
    double map_gamma = 0;
    
    for (j = 0; j < p; j++) {
      map_gamma = map_gamma + gamma(j)*(log(omega) - log(1 - omega));
      
      // SRT data likelihood ratio
      if (st_data) {
        for (i = 0; i < n; i++) {
          mu_sum_temp = sum(ct_p(i, _)*M(_, j));
          if (H(i, j) == 0){
            map_gamma = map_gamma + lgamma(Y(i, j) + phi(j)) - lgamma(phi(j))  - lgamma(Y(i, j) + 1) +
              phi(j)*(log(phi(j)) - log(mu_sum_temp*s(i) + phi(j))) + Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j)));
          }
        }
      }
      
      
      // scRNA-seq data likelihood ratio
      for (i = 0; i < m; i++) {
        double mu_temp = M(cell_type(i), j);
        if (Z(i, j) == 0){
          map_gamma = map_gamma + lgamma(X(i, j) + psi(j)) - lgamma(psi(j))  - lgamma(X(i, j) + 1) +
            psi(j)*(log(psi(j)) - log(mu_temp*v(i) + psi(j))) + X(i, j)*(log(mu_temp*v(i)) - log(mu_temp*v(i) + psi(j)));
        }
      }
    }
    
    map_gamma_store(t) = map_gamma;
    
    
    // Obtain the MAP estimates of z and gamma #################################
    map_z_store(t) = map_z;
    
    if (t == burn) {
      z_map = z;
      map_z_temp = map_z;
      
      gamma_map = gamma;
      map_gamma_temp = map_gamma;
      
    }else if (t >= burn) {
      if (map_z > map_z_temp) {
        z_map = z;
        map_z_temp = map_z;
      }
      if (map_gamma > map_gamma_temp) {
        gamma_map = gamma;
        map_gamma_temp = map_gamma;
      }
    }
    
    // Store the results
    Pi_store.slice(t) = as<arma::mat>(ct_p);
    Omega_store.slice(t) = as<arma::mat>(Omega);
    A_store.slice(t) = as<arma::mat>(A);
    M_store.slice(t) = as<arma::mat>(M);
    
    for (i = 0; i < n; i++){
      if (t >= burn){
        for(k = 0; k < K; k++) {
          Delta_ppi(i, k) = Delta_ppi(i, k) + Delta(i, k);
        }
      }
    }

    H_sum(t) = 0;
    for (i = 0; i < n; i++){
      H_sum(t) = H_sum(t) + H_sum_temp(i);
    }

    Z_sum(t) = 0;
    for (i = 0; i < m; i++){
      Z_sum(t) = Z_sum(t) + Z_sum_temp(i);

    }

    for (j = 0; j < p; j++){
      phi_store(t, j) = phi(j);
      psi_store(t, j) = psi(j);

      // for (i = 0; i < n; i++){
      //   H_store(i, j, t) = H(i, j);
      // }

      // for (ii = 0; ii < m; ii++){
      //   Z_store(ii, j, t) = Z(ii, j);
      // }
    }

    for (j = 0; j < p; j++){
     if (t >= burn){
       for(i = 0; i < n; i++) {
         H_ppi(i, j) = H_ppi(i, j) + H(i, j);
       }
       for(ii = 0; ii < m; ii++) {
         Z_ppi(ii, j) = Z_ppi(ii, j) + Z(ii, j);
       }
     }
    }

  }// End of iterations


  // Calculate acceptance rate, PPI and BF of gamma
  if (feature_selection){
    accept_gamma = accept_gamma/try_gamma;
    for (j = 0; j < p; j++) {
      gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
      gamma_BF(j) = gamma_ppi(j)/(1 - gamma_ppi(j))*(1 - omega)/omega;
    }
  }

  for(i = 0; i < n; i++) {
    for(k = 0; k < K; k++) {
      Delta_ppi(i, k) = Delta_ppi(i, k)/(iter - burn);
    }
  }
  
  for(j = 0; j < p; j++) {
    for(i = 0; i < n; i++) {
      H_ppi(i, j) = H_ppi(i, j)/(iter - burn);
    }
    for(ii = 0; ii < m; ii++) {
      Z_ppi(ii, j) = Z_ppi(ii, j)/(iter - burn);
    }
  }

  accept_alpha = accept_alpha/try_alpha;
  accept_omega = accept_omega/try_omega;
  accept_delta = accept_delta/try_delta;
  accept_mu = accept_mu/try_mu;
  accept_phi = accept_phi/try_phi;
  accept_psi = accept_psi/try_psi;

  return Rcpp::List::create(Rcpp::Named("Pi_store") = Pi_store, Rcpp::Named("A_store") = A_store, Rcpp::Named("Omega_store") = Omega_store, Rcpp::Named("Delta_ppi") = Delta_ppi,
                            Rcpp::Named("z_store") = z_store, Rcpp::Named("z_map") = z_map, 
                            Rcpp::Named("gamma_store") = gamma_store, Rcpp::Named("gamma_sum") = gamma_sum, Rcpp::Named("gamma_ppi") = gamma_ppi, Rcpp::Named("gamma_BF") = gamma_BF,
                            Rcpp::Named("M_store") = M_store,
                            Rcpp::Named("H_ppi") = H_ppi, Rcpp::Named("Z_ppi") = Z_ppi,
                            Rcpp::Named("phi_store") = phi_store, Rcpp::Named("psi_store") = psi_store,
                            Rcpp::Named("accept_gamma") = accept_gamma, Rcpp::Named("accept_alpha") = accept_alpha,
                            Rcpp::Named("accept_mu") = accept_mu, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("accept_psi") = accept_psi);
  
} // End of function


// Used functions ##############################################################
int num(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      count++;
    }
  }
  return count;
}

arma::uvec which(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  int m = num(x, c);
  arma::uvec index(m);
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      index(count) = i;
      count++;
    }
  }
  return index;
}

Rcpp::Environment base("package:base");
Function do_unique = base["unique"];
