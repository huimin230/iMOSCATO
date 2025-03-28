#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <RcppEigen.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppEigen)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List imoscato(arma::mat Y, arma::mat B, IntegerVector z, NumericVector s, arma::mat P,  int iter, int burn) {
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;
  int D = Rcpp::unique(z).size();
  int K = B.n_rows;
  
  bool st_nb = true;
  bool zero_inflation = true;
  bool find_domain = true;
  
  // Markov random field (MRF) prior setting
  IntegerVector e(D);
  for (int d = 0; d < D; d++){
    e(d) = 1;
  }
  int n_neighbor = P.n_cols;
  int min_element = 2;
  double f = 1;
  bool domain_data_likelihood = false;

  // Hyperparamters
  double delta_p = 0.5;
  double a_alpha = 0.001;
  double b_alpha = 0.001;
  double a_phi = 0.001;
  double b_phi = 0.001;
  double theta = 0.5;

  // Tuning parameters
  double tau_omega = 1;
  double tau_alpha = 1;
  double tau_phi = 1;
  
  // Set temporary variables
  int t, i, ii, j, jj, k, kk, ee, d, dd, num = 10, delta_new,H_temp, neighbor_index;
  double epsilon, hastings, accept_omega = 0, accept_delta = 0, accept_alpha = 0, accept_phi = 0;
  double try_omega = 0, try_delta = 0, try_alpha = 0, try_phi = 0;
  double h_temp, mu_sum_temp, mu_sum_new, omega_temp, omega_new, alpha_temp, alpha_new, phi_new, phi_temp;
  double map_pi_temp, map_alpha_temp, map_z_temp, map_phi_temp;
  double max_temp, sum_temp;
  IntegerVector z_map;
  NumericMatrix Pi_map(n, K);
  NumericMatrix Omega_map(n, K);
  NumericMatrix A_map(D, K);
  NumericVector phi_map(p);
  NumericMatrix A_new(D, K);
  IntegerVector H_sum_temp(n, 0);
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
  // arma::cube H_store(n, p, iter);
  NumericMatrix H_ppi(n, p);
  IntegerVector H_sum(iter);
  NumericMatrix phi_store(iter, p);
  
  NumericVector map_pi_store(iter);
  NumericVector map_alpha_store(iter);
  NumericVector map_z_store(iter);
  NumericVector map_phi_store(iter);
  
  // Initialization
  NumericMatrix A(D, K);
  for (k = 0; k < K; k++){
    for (d = 0; d < D; d++){
      A(d, k) = 10;
    }
  }
  
  NumericMatrix ct_p(n, K);
  NumericMatrix Omega(n, K);
  IntegerMatrix Delta(n, K);
  
  for (i = 0; i < n; i++){
    double sum_omega = 0;
    
    for (k = 0; k < K; k++){
      Omega(i, k) = 10;
      sum_omega += Omega(i, k);
    }
    
    for (k = 0; k < K; k++){
      ct_p(i, k) = Omega(i, k)/sum_omega;
      
    }
  }

  
  for (i = 0; i < n; i++) {
    for (k = 0; k < K; k++) {
      Delta_ppi(i, k) = 0;
      
      if (Omega(i, k) != 0) {
        Delta(i, k) = 0;
      } else {
        Delta(i, k) = 0;
      }
    }
  }
  

  // Spatial domains
  NumericVector prob_init(D);
  IntegerVector domains(D);
  
  for (d = 0; d < D; d++){
    domains(d) = d;
    prob_init(d) = pow(D, -1);
  }
  
  NumericMatrix M(K, p);
  for (j = 0; j < p; j++){
    for (k = 0; k < K; k++){
      M(k, j) = B(k,j);
    }
  }
  
  IntegerMatrix H(n, p);
  NumericVector phi(p);
  double phi_start;

  if (st_nb){
    phi_start = 10;
  }else{
    phi_start = 100000000;
  }
  
  for (j = 0; j < p; j++) {
  phi(j) = phi_start;

    for (i = 0; i < n; i++) {
      H_ppi(i, j) = 0;
      if (Y(i, j) != 0) {
        H(i, j) = 0;
      } else {
        H(i, j) = 0;
        H_sum_temp(i) += H(i, j);
      }
    }
  }

  IntegerVector clusters(K-1);
  for (k = 0; k < (K-1); k++){
    clusters(k) = k;
  }
  
  
  
  
  // Start MCMC algorithms
  for (t = 0; t < iter; t++){
    
    // Gipps sampler to update the zero-inflation parameter H in SRT data ###
    if (zero_inflation) {
      for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
          // only update for y_{ij} = 0
          if (Y(i, j) == 0) {
            NumericVector prob_temp(2);
            mu_sum_temp =  sum(ct_p(i, _)*M(_, j));
            
            prob_temp(0) = phi(j)*(log(phi(j)) - log(s(i)*mu_sum_temp + phi(j))) + log(1 - theta);
            prob_temp(1) = log(theta);
            
            // Numerical stability improvement
            max_temp = max(prob_temp);
            prob_temp(0) = exp(prob_temp(0) - max_temp);
            prob_temp(1) = exp(prob_temp(1) - max_temp);
            
            sum_temp = prob_temp(0) + prob_temp(1);
            prob_temp(0) = prob_temp(0)/sum_temp;
            prob_temp(1) = prob_temp(1)/sum_temp;
            
            H_temp = H(i, j);
            H(i, j) = rbinom(1, 1, prob_temp(1))(0);
            H_sum_temp(i) = H_sum_temp(i) - H_temp + H(i, j);
          }
        }
      }
    }  // End of update of the zero-inflation parameter H in SRT data


    // RWMH to update of the dispersion parameter Phi in SRT data ############
    if (st_nb) {
      for (j = 0; j < p; j++) {
        phi_temp = phi(j);
        phi_new = r_truncnorm(phi_temp, tau_phi, 1, 1000);
        // phi_new = exp(rnorm(1, log(phi_temp), tau_phi)(0));
        
        hastings = 0;

        for (i = 0; i < n; i++) {
          mu_sum_temp =  sum(ct_p(i, _)*M(_, j));

          // SRT data likelihood ratio
          if (H(i, j) == 0){
            hastings += lgamma(Y(i, j) + phi_new) - lgamma(phi_new) +
              phi_new*(log(phi_new) - log(mu_sum_temp*s(i) + phi_new)) + 
              Y(i, j)*( - log(mu_sum_temp*s(i) + phi_new));
            hastings -= (lgamma(Y(i, j) + phi_temp) - lgamma(phi_temp) +
              phi_temp*(log(phi_temp) - log(mu_sum_temp*s(i) + phi_temp)) + 
              Y(i, j)*( - log(mu_sum_temp*s(i) + phi_temp)));
          }
        }

        // Prior ratio
        hastings += (a_phi - 1)*log(phi_new) - b_phi*phi_new;
        hastings -= ((a_phi - 1)*log(phi_temp) - b_phi*phi_temp);

        // Check if accept the proposed new value
        if (t >= burn){
          try_phi++;
        }
        
        if (hastings >= log(double(rand()%10001)/10000)){
          phi(j) = phi_new;

          if (t >= burn){
            accept_phi++;
          }
        }

      }
    } // End of updating the dispersion parameter Phi in SRT data


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
        NumericVector loglklh(D, 0.0);
        epsilon = 1e-20;

        // if (domain_data_likelihood) {
        //   NumericMatrix st_M(n, p);
        //   for (ii = 0; ii < n; ii++) {
        //     for (j = 0; j < p; j++){
        //       st_M(ii, j) = sum(ct_p(ii, _)*M(_, j));
        //     }
        //   }
        //   
        //   NumericMatrix st_B(D, p);
        //   for (dd = 0; dd < D; dd++) {
        //     for (j = 0; j < p; j++){
        //       // Need to properly compute mean for spots in domain d
        //       double sum_domain = 0.0;
        //       int count_domain = 0;
        //       
        //       // Calculate mean of normalized expression for spots in domain d
        //       for (ii = 0; ii < n; ii++) {
        //         if (z(ii) == dd) {
        //           sum_domain += st_M(ii, j);
        //           count_domain++;
        //         }
        //       }
        //       
        //       // Avoid division by zero
        //       if (count_domain > 0) {
        //         st_B(dd, j) = sum_domain / count_domain;
        //       } else {
        //         st_B(dd, j) = 0.0;
        //       }
        //     }
        //   }
        // }
        

        for (d = 0; d < D; d++){
          // if (domain_data_likelihood) {
          //   // Evaluate data likelihood under this domain assignment
          //   for (j = 0; j < p; j++){
          //     if (H(i, j) == 0){
          //       double mu = st_B(z(i), j)+epsilon;
          //       
          //       // Add data likelihood to log probability
          //       loglklh(d) += phi(j)*(log(phi(j)) - log(mu*s(i) + phi(j))) +
          //         Y(i, j)*(log(mu*s(i)) - log(mu*s(i) + phi(j)));
          //     }
          //   }
          // }
          
          for (k = 0; k < K; k++){
            if (Delta(i, k) == 0) {
              
              // Conditional prior ratio of Omega_i
              loglklh(d) += -lgamma(A(d, k)) + (A(d, k) - 1)*log(Omega(i, k)+epsilon) - Omega(i, k);
              
            }
            
          }
          
          // Spatial MRF prior
          loglklh(d) += e(d) + f*(neighbor_temp(d));

        }
        
      
        // Normalize to avoid numerical issues
        double max_loglklh = loglklh(0);
        for (d = 1; d < D; d++) {
          if (loglklh(d) > max_loglklh) {
            max_loglklh = loglklh(d);
          }
        }
        
        NumericVector prob(D);
        double sum_prob = 0;
        for (d = 0; d < D; d++) {
          prob(d) = exp(loglklh(d) - max_loglklh);
          sum_prob += prob(d);
        }
        
        // Normalize probabilities
        for (d = 0; d < D; d++) {
          prob(d) /= sum_prob;
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


    // RWMH to update the concentration parameter A  ###########################
    for (d = 0; d < D; d++){
      for (k = 0; k < K; k++){
        alpha_temp = A(d, k);

        for (kk = 0; kk < K; kk++){
          A_new(d, kk) = A(d, kk);
        }

        alpha_new = r_truncnorm(alpha_temp, tau_alpha, 1, 1000);
        // alpha_new = exp(rnorm(1, log(alpha_temp), tau_alpha)(0));
        
        A_new(d, k) = alpha_new;

        hastings = 0;

        // Conditional prior ratio
        for (i = 0; i < n; i++){
          if (z(i) == d && Delta(i, k) == 0){
            hastings += (-lgamma(alpha_new) + (alpha_new - 1)*log(Omega(i, k)+epsilon));
            hastings -= (-lgamma(alpha_temp) + (alpha_temp - 1)*log(Omega(i, k)+epsilon));
          }
        }

        // Prior ratio
        hastings += (a_alpha - 1)*log(alpha_new) - b_alpha*alpha_new;
        hastings -= ((a_alpha - 1)*log(alpha_temp) - b_alpha*alpha_temp);

        if (t >= burn){
          try_alpha++;
        }
        
        if (hastings >= log(double(rand()%10001)/10000)) {
          A(d, k) = alpha_new;

          if (t >= burn){
            accept_alpha++;
          }
        }


      }


    } // End of updating concentration parameter A #############################
    

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
        double a_omega =  A(z(i), k);

        if (delta_new == 0){
          omega_new = rgamma(1, a_omega, 1)(0);
          Omega_new(k) = omega_new;
          
          // Recalculate proportions
          double sum_omega_new = 0;
          for (kk = 0; kk < K; kk++){
            sum_omega_new += Omega_new(kk);
          }
          
          for (kk = 0; kk < K; kk++){
            pi_new(kk) = Omega_new(kk) / sum_omega_new;
          }
          
          // Calculate hastings ratio
          hastings = 0;
          
          // SRT data likelihood ratio
          for (j = 0; j < p; j++){
            if (H(i, j) == 0){
                mu_sum_new = sum(pi_new*M(_, j));
                mu_sum_temp = sum(pi_temp*M(_, j));
                
                hastings += phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + 
                  Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
                hastings -= (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + 
                  Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
              }
          }
          
          // Conditional prior of omega given delta
          hastings += (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_new) - omega_new);
          
          // Proposal density ratio
          hastings -= (-lgamma(a_omega) + (a_omega - 1)*log(omega_new) - omega_new);

          // Prior ratio of delta
          hastings += log(1 - delta_p);
          hastings -= log(delta_p);
          
          // Check if accept the proposed new values
          if (t >= burn) {
            try_delta++;
          }
          
          if (hastings >= log(double(rand()%10001)/10000)) {
            Delta(i, k) = delta_new;
            Omega(i, k) = omega_new;
            
            // Recalculate all proportions
            double sum_omega = 0;
            for (kk = 0; kk < K; kk++){
              sum_omega += Omega(i, kk);
            }
            
            for (kk = 0; kk < K; kk++){
              ct_p(i, kk) = Omega(i, kk) / sum_omega;
            }
            
            if (t >= burn) {
              accept_delta++;
            }
          }
          
        } // End of delete step
        
        // Add step
        else {
          omega_new = 0;
          Omega_new(k) = omega_new;
          
          // Recalculate proportions
          double sum_omega_new = 0;
          for (kk = 0; kk < K; kk++){
            sum_omega_new += Omega_new(kk);
          }
          
          for (kk = 0; kk < K; kk++){
            pi_new(kk) = Omega_new(kk) / sum_omega_new;
          }
          
          // Calculate hastings ratio
          hastings = 0;
          
          // SRT data likelihood ratio
          for (j = 0; j < p; j++){
            if (H(i, j) == 0){
              mu_sum_new = sum(pi_new*M(_, j));
              mu_sum_temp = sum(pi_temp*M(_, j));
              
              hastings += phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + 
                Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
              hastings -= (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + 
                Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
            }
          }
          
          // Conditional prior of omega given delta
          hastings -= (-lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(omega_temp) - omega_temp);
          
          // Proposal density ratio
          hastings += (-lgamma(a_omega) + (a_omega - 1)*log(omega_temp) - omega_temp);

          // Prior ratio of delta
          hastings += log(delta_p);
          hastings -= log(1 - delta_p);
          
          // Check if accept the proposed new values
          if (t >= burn) {
            try_delta++;
          }
          
          if (hastings >= log(double(rand()%10001)/10000)) {
            // For any i, delta_{ik} cannot be 1 for all k
            if (sum(Delta(i, _)) < K){// Ensure at least one active component
              Delta(i, k) = delta_new;
              Omega(i, k) = omega_new;
              
              // Recalculate all proportions
              double sum_omega = 0;
              for (kk = 0; kk < K; kk++){
                sum_omega += Omega(i, kk);
              }
              
              for (kk = 0; kk < K; kk++){
                ct_p(i, kk) = Omega(i, kk) / sum_omega;
              }
              
              if (t >= burn) {
                accept_delta++;
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
          omega_new = r_truncnorm(omega_temp, tau_omega, 1, 1000);
          // omega_new = exp(rnorm(1, log(omega_temp), tau_omega)(0));
          
          Omega_new(k) = omega_new;
          
          // Recalculate proportions
          double sum_omega_new = 0;
          for (kk = 0; kk < K; kk++){
            sum_omega_new += Omega_new(kk);
          }
          
          for (kk = 0; kk < K; kk++){
            pi_new(kk) = Omega_new(kk) / sum_omega_new;
          }
          
          hastings =  0;
          
          // SRT data likelihood ratio
          for (j = 0; j < p; j++){
            if (H(i, j) == 0){
              mu_sum_new = sum(pi_new*M(_, j));
              mu_sum_temp = sum(pi_temp*M(_, j));
              
              hastings += phi(j)*(- log(mu_sum_new*s(i) + phi(j))) + 
                Y(i, j)*(log(mu_sum_new*s(i)) - log(mu_sum_new*s(i) + phi(j)));
              hastings -= (phi(j)*(- log(mu_sum_temp*s(i) + phi(j))) + 
                Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j))));
            }
          }
          
          // Prior ratio
          hastings += (A(z(i), k) - 1)*log(omega_new) - omega_new;
          hastings -= ((A(z(i), k) - 1)*log(omega_temp) - omega_temp);
          
          if (t >= burn){
            try_omega++;
          }
          
          // Check if accept the proposed new values
          if (hastings >= log(double(rand()%10001)/10000)) {
            Omega(i, k) = omega_new;
            
            // Recalculate all proportions
            double sum_omega = 0;
            for (kk = 0; kk < K; kk++){
              sum_omega += Omega(i, kk);
            }
            
            for (kk = 0; kk < K; kk++){
              ct_p(i, kk) = Omega(i, kk) / sum_omega;
            }
            
            if (t >= burn){
              accept_omega++;
            }
          }
          
          
        }
        
      }
    } // End of updating the cell type proportion Pi and Omega ###############
    

    // Monitor the process
    if((t*100/(iter-1)) == num) {
      Rcout<<num<< "% has been done\n";
      num = num + 10;
    }

    // Calculate the marginal posterior probability to obtain Pi and Omega MAP #
    double map_pi = 0;
    double map_data = 0;
    
    // SR data likelihood
    for (i = 0; i < n; i++){
      for (j = 0; j < p; j++){
        mu_sum_temp = sum(ct_p(i, _)*M(_, j));
        
        if (H(i, j) == 0){
          map_data += lgamma(Y(i, j) + phi(j)) - lgamma(phi(j))  - lgamma(Y(i, j) + 1) +
            phi(j)*(log(phi(j)) - log(mu_sum_temp*s(i) + phi(j))) + 
            Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j)));
        }
      }
    }
    
    
    for (i = 0; i < n; i++){
      // SR data likelihood
      for (j = 0; j < p; j++){
        mu_sum_temp = sum(ct_p(i, _)*M(_, j));
        if (H(i, j) == 0) {
          map_pi += lgamma(Y(i, j) + phi(j)) - lgamma(phi(j))  - lgamma(Y(i, j) + 1) +
            phi(j)*(log(phi(j)) - log(mu_sum_temp*s(i) + phi(j))) + 
            Y(i, j)*(log(mu_sum_temp*s(i)) - log(mu_sum_temp*s(i) + phi(j)));
        }
        for (k = 0; k < K; k++){
          // Prior of Omega
          if (Delta(i, k) == 0) {
            map_pi += - lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(Omega(i, k)+epsilon) - Omega(i, k);
          }
        }

      }
    }

    map_pi_store(t) = map_pi;
    
    
    // Calculate the marginal posterior probability of A to obtain domain MAP ###
    double map_alpha = 0;
    
    for (d = 0; d < D; d++){
      for (k = 0; k < K; k++){
        // Conditional prior
        for (i = 0; i < n; i++){
          if (z(i) == d && Delta(i, k) == 0){
            map_alpha += (-lgamma(A(d, k)) + (A(d,k) - 1)*log(Omega(i, k)+epsilon));
          }
        }
        
        // Prior
        map_alpha += (a_alpha - 1)*log(A(d, k)) - b_alpha*A(d, k);

      }
    } 
    
    map_alpha_store(t) = map_alpha;

    // Calculate the marginal posterior probability of z to obtain domain MAP ###
    double map_z = 0;

    for (i = 0; i < n; i++) {
      z_store(t, i) = z(i);

      // Prior of z
      int n_neighbor_temp = 0;

      // Compute I(z_i = z_i')
      for (int ii = 0; ii < n_neighbor; ii++){
        if (P(i, ii) != 0) {
          if(z(i) == z(P(i, ii) - 1)){
            n_neighbor_temp ++;

          }
        }
      }

      map_z += e(z(i)) + f*n_neighbor_temp;

      for (k = 0; k < K; k++) {
        if (Delta(i, k) == 0) {
          map_z += - lgamma(A(z(i), k)) + (A(z(i), k) - 1)*log(Omega(i, k)+epsilon) - Omega(i, k);
        }
      }
    }
    
    // Calculate the marginal posterior probability of phi to obtain MAP #######
    double map_phi = 0;
    
    for (j = 0; j < p; j++) {
      map_phi += (a_phi - 1)*log(phi(j)) - b_phi*phi(j);    
    }
    
    map_phi += map_data;
    map_phi_store(t) = map_phi;
    
    // Obtain the MAP estimates of pi/omega, and z #############################
    map_z_store(t) = map_z;

    if (t == (burn-1)) {
      Pi_map = clone(ct_p);
      Omega_map = clone(Omega);
      A_map = clone(A);
      z_map = clone(z);
      phi_map = clone(phi);

      map_pi_temp = map_pi;
      map_alpha_temp = map_alpha;
      map_z_temp = map_z;
      map_phi_temp = map_phi;

    }else {
      if (map_pi > map_pi_temp) {
        Pi_map = clone(ct_p);
        Omega_map = clone(Omega);
        map_pi_temp = map_pi;
      }
      
      if (map_alpha > map_alpha_temp) {
        A_map = clone(A);
        map_alpha_temp = map_alpha;
      }
      
      if (map_z > map_z_temp) {
        z_map = clone(z);
        map_z_temp = map_z;

      }

      if (map_phi > map_phi_temp) {
        phi_map = clone(phi);
        map_phi_temp = map_phi;
      }

      
    }

    // Store the results
    Delta_sum(t) = sum(Delta);
    Pi_store.slice(t) = as<arma::mat>(ct_p);
    Omega_store.slice(t) = as<arma::mat>(Omega);
    A_store.slice(t) = as<arma::mat>(A);

    for (i = 0; i < n; i++){
      if (t >= burn){
        for(k = 0; k < K; k++) {
          Delta_ppi(i, k) += Delta(i, k);
        }
      }
    }

    for (i = 0; i < n; i++){
      H_sum(t) += H_sum_temp(i);
    }

    for (j = 0; j < p; j++){
      phi_store(t, j) = phi(j);
      
     if (t >= burn){
       for(i = 0; i < n; i++) {
         H_ppi(i, j) += H(i, j);
       }
     }
    }
    
  }// End of iterations

  
  // Calculate acceptance rate
  for(i = 0; i < n; i++) {
    for(k = 0; k < K; k++) {
      Delta_ppi(i, k) = Delta_ppi(i, k)/(iter - burn);
    }
  }

  
  for(j = 0; j < p; j++) {
    for(i = 0; i < n; i++) {
      H_ppi(i, j) = H_ppi(i, j)/(iter - burn);
    }
  }

  accept_alpha = accept_alpha/try_alpha;
  accept_omega = accept_omega/try_omega;
  accept_delta = accept_delta/try_delta;
  accept_phi = accept_phi/try_phi;
  
  return Rcpp::List::create(Rcpp::Named("Pi_map") = Pi_map, Rcpp::Named("Omega_map") = Omega_map, Rcpp::Named("Delta_ppi") = Delta_ppi, Rcpp::Named("A_map") = A_map, Rcpp::Named("z_map") = z_map,Rcpp::Named("phi_map") = phi_map, Rcpp::Named("H_ppi") = H_ppi,
                            Rcpp::Named("Pi_store") = Pi_store, Rcpp::Named("Omega_store") = Omega_store, Rcpp::Named("A_store") = A_store, Rcpp::Named("z_store") = z_store,Rcpp::Named("phi_store") = phi_store,
                            Rcpp::Named("map_pi_store") = map_pi_store, Rcpp::Named("map_alpha_store") = map_alpha_store, Rcpp::Named("map_z_store") = map_z_store,Rcpp::Named("map_phi_store") = map_phi_store,
                            Rcpp::Named("Delta_sum") = Delta_sum, Rcpp::Named("H_sum") = H_sum,
                            Rcpp::Named("accept_omega") = accept_omega, Rcpp::Named("accept_delta") = accept_delta, Rcpp::Named("accept_alpha") = accept_alpha, Rcpp::Named("accept_phi") = accept_phi);
  
} // End of function

