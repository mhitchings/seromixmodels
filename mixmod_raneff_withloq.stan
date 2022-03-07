functions {
  
  real expit(real x) {
    real y;
    y = 1/(1+exp(-x));
    return y;
  }
  
}

// Input data.
data {
  int <lower=0> N; //the number of serosamples
  int <lower=0> r_b; //the number of covariates on response probability
  int <lower=0> r_g; //the number of covariates on response probability
  real y[N]; //antibody levels
  real reg_vars_b[N,r_b+1]; //covariates on seroresponse probability
  real reg_vars_g[N,r_g]; //covariates on seroresponse level
  int <lower = 0,upper=1 > model_nonresponse; //distribution of non-responders (0=normal, 1=loglogistic)
  int <lower = 0,upper=2 > model_response; //distribution of responders (0=normal,1=lognormal,2=gamma)
  int <lower = 0,upper=1 > negprior; //whether to use pre-pandemic negatives to inform priors
  real negprior_params[3]; //priors from pre-pandemic negatives
  int<lower=0> num_hh; //number of households
  int hhids[N]; //list of household ids for each household (1:num_hh)
  real<lower=0> loq; //lower limit of quantitation (if there is no LOQ, set it to 0)
  real gammadiff_prior;
}

// The parameters accepted by the model
parameters {
  // For prior on pre-pandemic negatives, these are uniform with the limits below
  real<lower=negprior_params[2]-negprior_params[3],upper=negprior_params[2]+negprior_params[3]> sigma0; //variance of non-response distribution
  real<lower=negprior_params[1]-negprior_params[3],upper=negprior_params[1]+negprior_params[3]> mean_nonresponse; //mean among non-responders
  real<lower=-5,upper=(1-model_response)> sigma;   //variance of seroresponse distribution
  real<lower=-4,upper=4> beta_covs[r_b+1]; // beta coefficients
  real<lower = 0,upper=10> gamma_diff; //baseline seroresponse
  real<lower=-3,upper=5> gamma_covs[r_g]; // beta coefficients
  real<lower=0,upper=2> raneff_sigma; // variance of random effect
  real hh_res[num_hh];
}

// The model to be estimated
model {
  real l_response;
  real l_nonresponse;
  real p;
  real gamma_i;
  real response_shape;
  real response_rate;
  real nonresponse_rate;
  
  
  // Prior on gamma_diff to help with convergence
  gamma_diff ~ normal(gammadiff_prior,0.5);
  
  // Household random effects
  for (i in 1:num_hh) {
    hh_res[i] ~ normal(0,raneff_sigma);
  }
  
  // Likelihood for each serosample
  for (i in 1:N) {
    
    // Probability of seropositivity for the individual (including random effect
    // and covariates)
    p = expit(hh_res[hhids[i]]+dot_product(beta_covs,reg_vars_b[i,]));
    
    // Mean seroresponse for a seronegative and seropositive individual including
    // covariates
    if (model_nonresponse==0) {
      // Normal distribution
      nonresponse_rate = exp(sigma0);
      if (r_g>0) {
        gamma_i = mean_nonresponse + gamma_diff + dot_product(gamma_covs,reg_vars_g[i,]);
      } else {
        gamma_i = mean_nonresponse + gamma_diff;
      }
    } else if (model_nonresponse==1) {
      // Log-logistic distribution
      nonresponse_rate = 1+exp(sigma0);
      if (r_g>0) {
        gamma_i = pi()/nonresponse_rate / (mean_nonresponse*sin(pi()/nonresponse_rate)) + gamma_diff + dot_product(gamma_covs,reg_vars_g[i,]);
      } else {
        gamma_i = pi()/nonresponse_rate / (mean_nonresponse*sin(pi()/nonresponse_rate)) + gamma_diff;
      }
    }
    
    // Shape parameters for the seroresponse among seropositives
    if (model_response==0) {
      // Normal distribution
      response_shape = gamma_i;
      response_rate = exp(sigma);
    } else if (model_response==1) {
      // Log-normal distribution
      response_shape = log(gamma_i);
      response_rate = exp(sigma);
    } else if (model_response==2) {
      // Gamma distribution
      response_shape = exp(sigma);
      response_rate = exp(sigma)/gamma_i;
    }
    
    // Now calculate likelihood for seronegative and seropositive compartments
    // given the y value and covariates
    
    if (y[i]<loq) {
      // Values below the LOQ
    
      // Seroresponse
      if (model_response==0) { // Normal
    
        l_response = normal_cdf(loq,response_shape,response_rate)*p;

      } else if (model_response==1) { // Log-normal

        l_response = lognormal_cdf(loq,response_shape,response_rate)*p;

      } else if (model_response==2) { //Gamma
    
        l_response = gamma_cdf(loq,response_shape,response_rate)*p;

      }

    // Non-responders
      if (model_nonresponse==0) { // Normal

        l_nonresponse = normal_cdf(loq,mean_nonresponse, nonresponse_rate)*(1-p);
        
      } else if (model_nonresponse==1) { // Log-logistic
      
        l_nonresponse = (1/(1+(mean_nonresponse*loq)^(-nonresponse_rate)))*(1-p);

      }
    } else {
      // Values above the LOQ
      
      // Seroresponse
      if (model_response==0) { // Normal
    
        l_response = exp(normal_lpdf(y[i] | response_shape, response_rate))*p;
        
      } else if (model_response==1) { // Log-normal

        l_response = exp(lognormal_lpdf(y[i] | response_shape, response_rate))*p;
        
      } else if (model_response==2) { //Gamma
    
        l_response = exp(gamma_lpdf(y[i] | response_shape, response_rate))*p;
        
      }

    // Non-responders
      if (model_nonresponse==0) { // Normal

        l_nonresponse = exp(normal_lpdf(y[i] | mean_nonresponse, nonresponse_rate))*(1-p);
        
      } else if (model_nonresponse==1) { // Log-logistic

        l_nonresponse = mean_nonresponse*nonresponse_rate*
          (mean_nonresponse*y[i])^(nonresponse_rate-1)/(1+(mean_nonresponse*y[i])^(nonresponse_rate))^2 * (1-p);
        
      }
    }
    
    // Final contribution to the likelihood is the log sum of the two compartments

    target += log(l_response + l_nonresponse);

  }
}

generated quantities {

  vector[N] log_lik;
  real l_response;
  real l_nonresponse;
  real p;
  real gamma_i;
  real response_shape;
  real response_rate;
  real nonresponse_rate;
  
  // Store log_lik for each individual so we can use the loo package
  
  
  for (i in 1:N) {
    
    p = expit(hh_res[hhids[i]]+dot_product(beta_covs,reg_vars_b[i,]));
    
    if (model_nonresponse==0) {
      nonresponse_rate = exp(sigma0);
      if (r_g>0) {
        gamma_i = mean_nonresponse + gamma_diff + dot_product(gamma_covs,reg_vars_g[i,]);
      } else {
        gamma_i = mean_nonresponse + gamma_diff;
      }
    } else if (model_nonresponse==1) {
      nonresponse_rate = 1+exp(sigma0);
      if (r_g>0) {
        gamma_i = pi()/nonresponse_rate / (mean_nonresponse*sin(pi()/nonresponse_rate)) + gamma_diff + dot_product(gamma_covs,reg_vars_g[i,]);
      } else {
        gamma_i = pi()/nonresponse_rate / (mean_nonresponse*sin(pi()/nonresponse_rate)) + gamma_diff;
      }
    }
    
    if (model_response==0) {
      response_shape = gamma_i;
      response_rate = exp(sigma);
    } else if (model_response==1) {
      response_shape = log(gamma_i);
      response_rate = exp(sigma);
    } else if (model_response==2) {
      response_shape = exp(sigma);
      response_rate = exp(sigma)/gamma_i;
    }
    
    if (y[i]<loq) {
      // Values below the LOQ
    
      // Seroresponse
      if (model_response==0) { // Normal
    
        l_response = normal_cdf(loq,response_shape,response_rate)*p;

      } else if (model_response==1) { // Log-normal

        l_response = lognormal_cdf(loq,response_shape,response_rate)*p;

      } else if (model_response==2) { //Gamma
    
        l_response = gamma_cdf(loq,response_shape,response_rate)*p;

      }

    // Non-responders
      if (model_nonresponse==0) { // Normal

        l_nonresponse = normal_cdf(loq,mean_nonresponse, nonresponse_rate)*(1-p);
        
      } else if (model_nonresponse==1) { // Log-logistic
      
        l_nonresponse = (1/(1+(mean_nonresponse*loq)^(-nonresponse_rate)))*(1-p);

      }
    } else {
      // Values above the LOQ
      
      // Seroresponse
      if (model_response==0) { // Normal
    
        l_response = exp(normal_lpdf(y[i] | response_shape, response_rate))*p;
        
      } else if (model_response==1) { // Log-normal

        l_response = exp(lognormal_lpdf(y[i] | response_shape, response_rate))*p;
        
      } else if (model_response==2) { //Gamma
    
        l_response = exp(gamma_lpdf(y[i] | response_shape, response_rate))*p;
        
      }

    // Non-responders
      if (model_nonresponse==0) { // Normal

        l_nonresponse = exp(normal_lpdf(y[i] | mean_nonresponse, nonresponse_rate))*(1-p);
        
      } else if (model_nonresponse==1) { // Log-logistic

        l_nonresponse = mean_nonresponse*nonresponse_rate*
          (mean_nonresponse*y[i])^(nonresponse_rate-1)/(1+(mean_nonresponse*y[i])^(nonresponse_rate))^2 * (1-p);
        
      }
    }

    log_lik[i] = log(l_response + l_nonresponse);

  }
}
