
functions{
  real a_calc_homog(real Source, real alpha_A,  real T){
      real c = alpha_A;
      real g = Source;
      
      real a = (g+c*T)/T;  
       return a;
    
  }
  
  }
  
data{

    real kcounts_y; 
    real alpha_A_data;
    real beta_data;
    real Source_data;
    real a_y_data;
    real a_t_data;
    real b_y_data;
    real b_t_data;
    
 
}

parameters{

    // priors for parameters
    real<lower = 0.001,    upper = 1>          alpha_A; //division of fast
    real<lower = 0.2,       upper = 0.5>      beta;
    real<lower = ((kcounts_y))*0.005,       upper = ((kcounts_y))*0.1>    Source;

    real<lower = 0,    upper = 1> a_y;
    real<lower = 0,    upper = .1> b_y;
    real<lower = 0,    upper = 1> a_t;
    real<lower = 0,    upper = .1> b_t;
}

  transformed parameters{

 
    real delta_A = a_calc_homog(Source,alpha_A,  ((kcounts_y)));
}

model{
    

    //priors for parameters
    alpha_A ~ lognormal(log(alpha_A_data),  1);
    beta    ~ lognormal(log(beta_data),     1);
    Source  ~ lognormal(log(Source_data),   1);
    
    //priors for precursors
    a_y    ~ lognormal(log(abs(a_y_data)),     0.1);
    b_y    ~ lognormal(log(abs(b_y_data)),      0.1);
    a_t    ~ lognormal(log(abs(a_t_data)),      0.1);
    b_t    ~ lognormal(log(abs(b_t_data)),   0.1);
    
    
}     

