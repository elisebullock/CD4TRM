functions{
#include functions.stan
    real expdecaysol(real t, real a, real b){
        return a*exp(-b*(t));
    }
    
    real expincsol(real t, real a, real b){
        return a*(1-exp(-b*(t-4)));
    }
    
  real intplus(int j, vector y, int numofki67int, int numofbrduint, int start) 
  {
      // Calculating a jumping index (calculates indexes down a ki67 strand)
      array[numofki67int+2] int indices;
      for (i in 1:(numofki67int+2)) {
          indices[i]=start+j+(numofbrduint+2)*(i-1);
        }
      real summed = sum(y[indices]);
      
      return summed;   
  }
   
  real a_calc_homog(real Source, real alpha_A,  real T){
      real c = alpha_A;
      real g = Source;
      
      real a = (g+c*T)/T;  
       return a;
    
  }
  
  
      vector KIAB_inb_yfp(real t,  vector y,
              // unpack parameters
              real alpha_A,
              real delta_A,
              real beta,
              real Source,
              int numofki67int,
              real a_y,
              real b_y,
              real kf,
              real a_y_data) 
      { // BRDU labeling of A and B
              int numki67comp = numofki67int+1;
              int switch1 = (2)*(numofki67int+2);
              real exp_decay_value;
            if (a_y_data==0.372){
                exp_decay_value = expincsol(t, a_y,b_y);}
            else
              {exp_decay_value = expdecaysol(t, a_y,b_y);}
            //real exp_decay_value = 0.5;
                
            vector[switch1] dy;
            
            // EM Fast Yfp+- Ki67+-int
            dy[1] = Source*exp_decay_value*kf - (delta_A+numki67comp*beta+alpha_A)*y[1] + 2*alpha_A*intplus(0,y,numofki67int,0,1); 
            
            dy[2] = Source*(1-exp_decay_value)*kf -(delta_A+alpha_A+numki67comp*beta)*y[2] + 2*alpha_A*intplus(1,y,numofki67int,0,1); 
                
            for (i in 0:(numofki67int-1)){
                for (j in 0:(1)){
                    dy[(2)*(1+i)+j+1] = -(delta_A+alpha_A+numki67comp*beta)*y[(2)*(1+i)+j+1]+numki67comp*beta*y[(2)*(i)+j+1];
             }}
                
            dy[(1+numofki67int)*(2)+1] = Source*exp_decay_value*(1-kf)-(delta_A+alpha_A)*y[(1+numofki67int)*(2)+1]+numki67comp*beta*y[(numofki67int)*(2)+1];
            dy[switch1]  =Source*(1-exp_decay_value)*(1-kf) -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2)*(1+numofki67int)];

          return dy;
      }
      
      vector KIAB_tmf_yfp(real t,  vector y,
              // unpack parameters
              real alpha_A,
              real delta_A,
              real beta,
              real Source,
              real eff,
              int numofki67int,
              real a_y,
              real b_y,
              real kf,
              real a_y_data) 
      { // BRDU labeling of A and B
              int numki67comp = numofki67int+1;
              int switch1 = (2)*(numofki67int+2);
              real exp_decay_value;
              if (a_y_data==0.372){
                  //exp_decay_value = a_y *(1-exp(-b_y*(t)));
                  exp_decay_value=0;}
              else
                {exp_decay_value =expdecaysol(t, a_y,b_y);}
            //real exp_decay_value = 0.5;
                
            vector[switch1] dy;
            
            // EM Fast Yfp+- Ki67+-int
            dy[1] = eff*y[2]+Source*exp_decay_value*kf - (delta_A+numki67comp*beta+alpha_A)*y[1] + 2*alpha_A*intplus(0,y,numofki67int,0,1); 
            
            dy[2] = -eff*y[2]+Source*(1-exp_decay_value)*kf -(delta_A+alpha_A+numki67comp*beta)*y[2] + 2*alpha_A*intplus(1,y,numofki67int,0,1); 
                
            for (i in 0:(numofki67int-1)){
                for (j in 0:(1)){
                    dy[(2)*(1+i)+j+1] = ((-j*2)+1)*eff*y[(2)*(1+i)+2]-(delta_A+alpha_A+numki67comp*beta)*y[(2)*(1+i)+j+1]+numki67comp*beta*y[(2)*(i)+j+1];
             }}
                
            dy[(1+numofki67int)*(2)+1] = Source*exp_decay_value*(1-kf)-(delta_A+alpha_A)*y[(1+numofki67int)*(2)+1]+numki67comp*beta*y[(numofki67int)*(2)+1];
            dy[switch1]  =Source*(1-exp_decay_value)*(1-kf) -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2)*(1+numofki67int)];

          return dy;
      }      
    
      vector equilibrium_ode_yfp(real t,  vector y,
              // unpack parameters
              real alpha_A,
              real delta_A,
              real beta,
              real Source,
              int numofki67int,
              real kf) 
      { // Calculating the steady state values
              int numki67comp = numofki67int+1;
              int switch1 = (2)*(numofki67int+2);
              int exp_decay_value = 0;
                  
              vector[switch1] dy;
              
            // EM Fast Yfp+- Ki67+-int
            dy[1] = Source*exp_decay_value*kf - (delta_A+numki67comp*beta+alpha_A)*y[1] + 2*alpha_A*intplus(0,y,numofki67int,0,1); 
            
            dy[2] = Source*(1-exp_decay_value)*kf -(delta_A+alpha_A+numki67comp*beta)*y[2] + 2*alpha_A*intplus(1,y,numofki67int,0,1);
                
            for (i in 0:(numofki67int-1)){
                for (j in 0:(1)){
                    dy[(2)*(1+i)+j+1] = -(delta_A+alpha_A+numki67comp*beta)*y[(2)*(1+i)+j+1]+numki67comp*beta*y[(2)*(i)+j+1];
               }}
                
            dy[(1+numofki67int)*(2)+1] = Source*exp_decay_value*(1-kf)-(delta_A+alpha_A)*y[(1+numofki67int)*(2)+1]+numki67comp*beta*y[(numofki67int)*(2)+1];
            dy[switch1]  =Source*(1-exp_decay_value)*(1-kf) -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2)*(1+numofki67int)];
            

          return  dy;
      }

     vector KIAB_inb_tom(real t,  vector y,
             // unpack parameters
             real alpha_A,
             real delta_A,
             real Source,
             real a_t,
             real b_t) 
     { // tom labeling of A and B
            //real exp_decay_value = (0.04+0.01017982*pow(t,4.74886179)*exp(-0.75671596*t));
            real exp_decay_value = expdecaysol(t, a_t,b_t);
               	
 
            vector[2] dy;
            
            // EM Fast tom+- Ki67+-int
            dy[1] = Source*exp_decay_value - (delta_A+alpha_A)*y[1] + 2*alpha_A*y[1]; 
            
            dy[2] = Source*(1-exp_decay_value) -(delta_A+alpha_A)*y[2] + 2*alpha_A*y[2]; 

   
         return dy;
     }
  
  
     vector equilibrium_ode_tom(real t,  vector y,
             // unpack parameters
             real alpha_A,
             real delta_A,
             real Source,
             real efft) 
     { // Calculating the steady state values 
            real exp_decay_value = efft;
                
                  
            vector[2] dy;
            
            // EM Fast tom+- Ki67+-int
            dy[1] = Source*exp_decay_value - (delta_A+alpha_A)*y[1] + 2*alpha_A*y[1]; 
            
            dy[2] = Source*(1-exp_decay_value) -(delta_A+alpha_A)*y[2] + 2*alpha_A*y[2]; 

         return  dy;
     }

    vector effic(vector initial_conditions, int numofki67int, real effy){
        int switch1 = (2)*(numofki67int+2);
        vector[switch1]  initial_cond_eff = initial_conditions;
        
        for (i in 1:((numofki67int+1))){
            initial_cond_eff[i*2-1]= initial_conditions[i*2]*effy;
            initial_cond_eff[i*2]= initial_conditions[i*2]*(1-effy);

        }        
        return initial_cond_eff;
    }
    
    vector effic2(vector initial_conditions, int numofki67int, real effy){
        int switch1 = (2)*(numofki67int+2);
        vector[switch1]  initial_cond_eff = initial_conditions;
        
        for (i in 1:((numofki67int+1))){
            initial_cond_eff[i*2-1]= initial_conditions[i*2-1]+initial_conditions[i*2]*effy;
            initial_cond_eff[i*2]= initial_conditions[i*2]*(1-effy);

        }        
        return initial_cond_eff;
    }    
  

  array[] vector solveode_y(vector init_cond, array[] real time_index_y, array[] real time_index_equilibrium,
        real alpha_A,
        real delta_A,
        real beta,
        real Source,
        int numofki67int,
        int switch1,
        int numObs_y,
        real a_y,
        real b_y,
        real kf,
        real timestart,
        real a_y_data
        )
    {
        
        // code for either non-swtiching or for real data
        array[numObs_y] vector[switch1] k_hat_y = ode_rk45(KIAB_inb_yfp, init_cond, timestart-1, time_index_y, alpha_A, delta_A, beta, Source, numofki67int, a_y, b_y, kf, a_y_data);
                
    return k_hat_y;
  }
   
 array[] vector solveode_t(vector init_cond, array[] real time_index_t, array[] real time_index_equilibrium,
       real alpha_A,
       real delta_A,
       real Source,
       int numObs_t,
       real a_t,
       real b_t,
       real timestart
       )
   {
       
       // code for either non-swtiching or for real data
       array[numObs_t] vector[2] k_hat_t = ode_rk45(KIAB_inb_tom, init_cond, timestart-1, time_index_t, alpha_A, delta_A, Source, a_t, b_t);
               
 
   return k_hat_t;
 }    
     
}
data{
  int<lower  = 1>    numObs_y;        // number of observations for donor fractions may be different that cell counts
  int switch1;                      // Switching array point between A and B cells
  array[numObs_y] real time_index_y;    // time array corresponding to data points
  array[numObs_y] real kcounts_y;       // Counts of all cells
  /*
  array[numObs_y] int yfphi_kihi_c;
  array[numObs_y] int yfphi_kilo_c;
  array[numObs_y] int yfplo_kilo_c;
  array[numObs_y] int yfplo_kihi_c;*/
  array[100] real fulltime_y;
  
  vector[switch1] init_cond_obs;  // Initial conditions used for equilibrium setting
  vector[2] init_cond_obs_t;  // Initial conditions used for equilibrium setting

  int numofki67int;                 // Num of ki67 intermediates
  array[1] real time_index_equilibrium; // time array for steady state equili
  array[1] real time_index_yfpeff;
  int counts_t0_y;
  real scalescale_y;
  real stdcounts_y;
    int numObs_t_pre;
  int<lower  = 1>    numObs_t;        // number of observations for donor fractions may be different that cell counts
  array[numObs_t] real time_index_t;    // time array corresponding to data points
  array[numObs_t_pre] real time_index_t_pre;
  real timestart;
  /*array[numObs_t] real kcounts_t;       // Counts of all cells
  array[numObs_t] int tomhi_c;
  array[numObs_t] int tomlo_c;*/
  int counts_t0_t;
  real scalescale_t;
  real stdcounts_t; 
  array[100] real fulltime_t;
  array[numObs_y] real ki67inyfphi;
  array[numObs_y] real ki67inyfplo;
  array[numObs_y] real yfphi;
  array[numObs_t] real tomhi;
array[numObs_y] real preyfp;
array[numObs_t_pre] real pretom;
  real alpha_A_data;
  real delta_A_data;
  real beta_data;
  real Source_data;
  real effy_data;
  real efft_data;
    
  real a_y_data;
  real a_t_data;
  real b_y_data;
  real b_t_data; 
  real kf; 
  //real Source_hi;
  //real Source_lo;
   
  //real Source;                      // Naive pool
  //real fs;
  //real eff;                         // Efficiency of BRDU uptake upon division
  //real mu;                          // Degradation of BRDU during unlabeling phase
  //real alpha_A;
  //real alpha_B;
  //real delta_A;
  //real delta_B;
  //real gamma;
  //real beta;

  
  }

/*
transformed data{
  array[numObs_y,4] int data_fractional_y;

  array[numObs_y] real ki67inyfphi;
  array[numObs_y] real ki67inyfplo;
  array[numObs_y] real yfphi;
  
  
  // data transformations    
  for (i in 1:numObs_y){
      data_fractional_y[i,1] = yfphi_kilo_c[i];
      data_fractional_y[i,2] = yfphi_kihi_c[i];
      data_fractional_y[i,3] = yfplo_kilo_c[i];
      data_fractional_y[i,4] = yfplo_kihi_c[i];
      
      ki67inyfphi[i] = yfphi_kihi_c[i]* 1.0/(yfphi_kihi_c[i]+yfphi_kilo_c[i]);
      ki67inyfplo[i] = yfplo_kihi_c[i]* 1.0/(yfplo_kihi_c[i]+yfplo_kilo_c[i]) ;
      yfphi[i] = (yfphi_kihi_c[i]+yfphi_kilo_c[i])* 1.0/(yfphi_kihi_c[i]+yfphi_kilo_c[i]+yfplo_kihi_c[i]+yfplo_kilo_c[i]) ;
  }
  
  array[numObs_t,2] int data_fractional_t;
  array[numObs_t] real tomhi;
  
  // data transformations    
  for (i in 1:numObs_t){
      data_fractional_t[i,1] = tomhi_c[i];
      data_fractional_t[i,2] = tomlo_c[i];
      
      tomhi[i] = (tomhi_c[i])* 1.0/(tomhi_c[i]+tomlo_c[i]) ;
       
  }
}

*/

parameters{

    // priors for parameters  

    //real<lower = delta_A_data/10,    upper = delta_A_data*10>   delta_A; //death of fast
    real<lower = 0.001,    upper = 1>          alpha_A; //division of fast
    real<lower = 0.2,       upper = 0.5>      beta;
    real<lower = (mean(kcounts_y))*0.005,       upper = (mean(kcounts_y))*0.1>    Source;
    //real<lower = Source_lo,       upper = Source_hi>    Source;

    real<lower = 0.9,        upper = 0.999>                effy;
    real<lower = 0.6,        upper = 0.95>                efft;
    
    real<lower = 0,    upper = 1> a_y;
    real<lower = 0,    upper = .1> b_y;
    real<lower = 0,    upper = 1> a_t;
    real<lower = 0,    upper = .1> b_t;


/*
real a_y;
real b_y;
real a_t;
real b_t;

real<lower = (a_y_data/10),    upper = (a_y_data*10)> a_y;
real<lower = (b_y_data/10),    upper = (b_y_data*10)> b_y;
real<lower = (a_t_data/10),    upper = (a_t_data*10)> a_t;
real<lower = (b_t_data/10),    upper = (b_t_data*10)> b_t;


    //real<lower = delta_A_data/10,    upper = delta_A_data*10>   delta_A; //death of fast
    real<lower = alpha_A_data/10,    upper = alpha_A_data*10>          alpha_A; //division of fast
    real<lower = beta_data/10,       upper = beta_data*10>      beta;
    real<lower = Source_data/10,       upper = Source_data*10>    Source;
    //real<lower = Source_lo,       upper = Source_hi>    Source;

    real<lower = eff_data/10,        upper = 0.999>                eff;
*/

    
    // priors for data likelihoods
    
    real<lower = 0>           scale_y;
    real<lower = 0>           phi_inv_y;
    real<lower = 0>           scale_t;    
    real<lower = 0>           phi_inv_t;
    
    real<lower = 0>    sigma_tomhi;
    real<lower = 0>    sigma_khiyfphi ;
    real<lower = 0>    sigma_khiyfplo ;
    real<lower = 0>    sigma_yfphi;
    real<lower = 0>    sigma_preyfp;
    real<lower = 0>    sigma_pretom;

  }

  transformed parameters{

 
    real delta_A = a_calc_homog(Source,alpha_A,  (mean(kcounts_y)));


    vector<lower = 0>[numObs_y]   c_total_y;        // Total cell counts
    vector<lower = 0,upper=1>[numObs_y]  f_yfphi_kihi;   // fraction of BRDuHi ki67Hi
    vector<lower = 0,upper=1>[numObs_y]  f_yfphi_kilo;   // fraction of BRDuHi ki67Lo 
    vector<lower = 0,upper=1>[numObs_y]  f_kihi_yfphi;
    vector<lower = 0,upper=1>[numObs_y]  f_kihi_yfplo;    
    vector<lower = 0,upper=1>[numObs_y]  f_yfphi;        // fraction of BRDUHi
    vector<lower = 0,upper=1>[numObs_y]  f_kihi_y;         // fraction of ki67Hi 
  
    vector<lower = 0>[numObs_y]  c_kihi_y;        // fraction of ki67Hi 
    vector<lower = 0>[numObs_y]  c_kilo_y;        // fraction of ki67Hi 
    vector<lower = 0>[numObs_y]  c_yfplo;       // fraction of ki67Hi 
    vector<lower = 0>[numObs_y]  c_yfphi_kihi;  // fraction of ki67Hi
    vector<lower = 0>[numObs_y]  c_yfphi_kilo;  // fraction of ki67Hi
    vector<lower = 0>[numObs_y]  c_yfplo_kilo;  // fraction of ki67Hi  
    vector<lower = 0>[numObs_y]  c_yfplo_kihi;  // fraction of ki67Hi  
    
    array[numObs_y] vector[4] norm_c_total_y;    // Simplex to store the normalized quadrant counts

    array[1] vector[switch1] init_conditions = ode_bdf(equilibrium_ode_yfp, init_cond_obs, 0.0, time_index_equilibrium, alpha_A, delta_A, beta, Source,  numofki67int, kf);
    
    vector[switch1] init_condit_eff = effic(init_conditions[1],numofki67int,effy);
        /*
    array[1] vector[switch1] init_conditions2 = ode_bdf(KIAB_inb_yfp, init_condit_eff, 0.0, time_index_yfpeff, alpha_A, delta_A, beta, Source, numofki67int, a_y, b_y, kf);
    vector[switch1] init_condit_eff2 = effic2(init_conditions2[1],numofki67int,eff);
    
    array[numObs_y] vector[switch1] k_hat_y = solveode_y(init_condit_eff2, time_index_y, time_index_equilibrium,
      alpha_A,
         delta_A,
         beta,
         Source,
         numofki67int,
         switch1,
         numObs_y,
         a_y,
         b_y,
         kf); */
             
    array[1] vector[switch1] init_conditions2 = ode_bdf(KIAB_tmf_yfp, init_condit_eff, 0.0, time_index_yfpeff, alpha_A, delta_A, beta, Source, effy,numofki67int, a_y, b_y, kf, a_y_data);

    array[numObs_y] vector[switch1] k_hat_y = solveode_y(init_conditions2[1], time_index_y, time_index_equilibrium,
      alpha_A,
         delta_A,
         beta,
         Source,
         numofki67int,
         switch1,
         numObs_y,
         a_y,
         b_y,
         kf,
         timestart,
         a_y_data); 

    for (i in 1:numObs_y){
        // calculating
        c_yfplo_kilo[i] =  k_hat_y[i,switch1];
        
        c_kihi_y[i] = sum(k_hat_y[i,1:((2)*(1+numofki67int))]);
        c_yfplo[i] = intplus(1, k_hat_y[i,:], numofki67int,  0, 1);
        
        // total counts
        c_total_y[i] = sum(k_hat_y[i, 1:switch1]);
        
        c_kilo_y[i] = c_total_y[i] - c_kihi_y[i];
        c_yfphi_kihi[i] = c_total_y[i] - c_yfplo[i] - sum(k_hat_y[i,((2)*(1+numofki67int)+1):(switch1-1)]);

        // fractions of yfphi cells in ki67Hi
        f_yfphi_kihi[i] = (c_yfphi_kihi[i])/c_kihi_y[i];
      
        // fractions of yfphi cells in KI67lo
        c_yfphi_kilo[i] = (sum(k_hat_y[i,(1+(2)*(1+numofki67int)):(switch1-1)]));
        f_yfphi_kilo[i] = c_yfphi_kilo[i]/c_kilo_y[i];
        
        c_yfplo_kihi[i] = (c_total_y[i]-c_yfphi_kilo[i]-c_yfphi_kihi[i]-c_yfplo_kilo[i]);
        
        // fractions of ki67hi cells in yfpHi
        f_kihi_yfphi[i] = (c_yfphi_kihi[i])/(c_yfphi_kihi[i]+c_yfphi_kilo[i]);
      
        // fractions of ki6s7hi cells in yfplo
        f_kihi_yfplo[i] = c_yfplo_kihi[i]/(c_yfplo_kihi[i]+c_yfplo_kilo[i]);
                      
        // fractions of yfphi
        f_yfphi[i] = (c_total_y[i]-c_yfplo[i])/c_total_y[i];
        
        // fractions of ki67Hi
        f_kihi_y[i] = c_kihi_y[i]/c_total_y[i];
                
        norm_c_total_y[i,1] = c_yfphi_kilo[i]/c_total_y[i];
        norm_c_total_y[i,2] = c_yfphi_kihi[i]/c_total_y[i];
        norm_c_total_y[i,3] = c_yfplo_kilo[i]/c_total_y[i];
        norm_c_total_y[i,4] = c_yfplo_kihi[i]/c_total_y[i];  
    
    }

    vector<lower = 0>[numObs_t]   c_total_t;        // Total cell counts
    vector<lower = 0,upper=1>[numObs_t]  f_tomhi;        // fraction of BRDUHi

    vector<lower = 0>[numObs_t]  c_tomlo;       // fraction of ki67Hi  

    array[numObs_t] vector[2] norm_c_total_t;    // Simplex to store the normalized quadrant counts

    array[1] vector[2] init_conditions_t = ode_bdf(equilibrium_ode_tom, init_cond_obs_t, 0.0, time_index_equilibrium, alpha_A, delta_A, Source, efft);

    array[numObs_t] vector[2] k_hat_t = solveode_t(init_conditions_t[1], time_index_t, time_index_equilibrium,
      alpha_A,
         delta_A,
         Source,
         numObs_t,
         a_t,
         b_t,
         timestart); 

    for (i in 1:numObs_t){
        // calculating
        c_tomlo[i] = k_hat_t[i, 2];
        
        // total counts
        c_total_t[i] = sum(k_hat_t[i, 1:2]);
        
        // fractions of tomhi
        f_tomhi[i] = (c_total_t[i]-c_tomlo[i])/c_total_t[i];
                
        norm_c_total_t[i,1] = k_hat_t[i, 1]/c_total_t[i];
        norm_c_total_t[i,2] = k_hat_t[i, 2]/c_total_t[i]; 
    
    }
}

model{
    //priors for liklihoods
    
    scale_y ~ lognormal(log(stdcounts_y), scalescale_y);
    scale_t ~ lognormal(log(stdcounts_t), scalescale_t);
    phi_inv_y ~ exponential(5);
    phi_inv_t ~ exponential(5);
    
    sigma_tomhi~ normal(0, 1);
    sigma_khiyfphi ~ normal(0, 1);
    sigma_khiyfplo ~ normal(0, 1);
    sigma_yfphi~ normal(0, .1);
    
    sigma_preyfp ~ normal(0, 1);
    sigma_pretom ~ normal(0, 1);
        
    //priors for parameters
    alpha_A ~ lognormal(log(alpha_A_data),  .1);
    //delta_A ~ lognormal(log(delta_A_data),  .1);
    beta    ~ lognormal(log(beta_data),     0.1);
    effy    ~ lognormal(log(effy_data),      0.05);
    efft    ~ lognormal(log(efft_data),      0.05);
    Source  ~ lognormal(log(Source_data),   0.1);
    
    //priors for precursors
    a_y    ~ lognormal(log(abs(a_y_data)),     0.1);
    b_y    ~ lognormal(log(abs(b_y_data)),      0.1);
    a_t    ~ lognormal(log(abs(a_t_data)),      0.1);
    b_t    ~ lognormal(log(abs(b_t_data)),   0.1);
   
    
    // prior for total counts at t=0
    //counts_t0_y ~ lognormal(log(sum(init_conditions[1,:])), stdcounts_y*10/numObs_y);

    for (i in 1:numObs_y) {
        //data_fractional_y[i,:] ~ dirichlet_multinomial((norm_c_total_y[i,:]/phi_inv_y)); 
        logit(ki67inyfphi[i]) ~ normal(logit(f_kihi_yfphi[i]),sigma_khiyfphi);
        logit(ki67inyfplo[i]) ~ normal(logit(f_kihi_yfplo[i]),sigma_khiyfplo);
        logit(yfphi[i]) ~ normal(logit(f_yfphi[i]), sigma_yfphi);
        
        if (a_y_data==0.372){
            logit(preyfp[i]) ~ normal(logit(expincsol(time_index_y[i],inv_logit(a_y),inv_logit(b_y))), sigma_preyfp);}
        else{
            logit(preyfp[i]) ~ normal(logit(expdecaysol(time_index_y[i],inv_logit(a_y),inv_logit(b_y))), sigma_preyfp);}
   }
    
    // prior for total counts at t=0
    //counts_t0_t ~ lognormal(log(sum(init_conditions_t[1,:])), stdcounts_t*10/numObs_t);

    for (i in 1:numObs_t) {
        //data_fractional_t[i,:] ~ dirichlet_multinomial((norm_c_total_t[i,:]/phi_inv_t));  
        //logit(ki67intomhi[i]) ~ normal(logit(f_kihi_tomhi[i]),sigma_khitomhi);
        //logit(ki67intomlo[i]) ~ normal(logit(f_kihi_tomlo[i]),sigma_khitomlo);
        logit(tomhi[i]) ~ normal(logit(f_tomhi[i]), sigma_tomhi);  
        
    }
    
    for (i in 1:numObs_t_pre) {
                  
        logit(pretom[i]) ~ normal(logit(expdecaysol(time_index_t_pre[i],inv_logit(a_t),inv_logit(b_t))), sigma_pretom);      
    }
    
}


generated quantities{
    //posteriors
    real ppd_alpha_A = lognormal_rng(log(alpha_A_data), 0.5);
    real ppd_delta_A = lognormal_rng(log((delta_A_data)), 0.1);
    real ppd_beta = lognormal_rng(log(beta_data), 0.5);
    real ppd_Source = lognormal_rng(log(Source_data), 0.5);
    real ppd_effy = lognormal_rng(log(effy_data), 0.1);
    real ppd_efft = lognormal_rng(log(efft_data), 0.1);
    
    vector [numObs_y+numObs_t] log_lik;
 
    vector[numObs_y] log_lik_y;
    vector[numObs_y] data_counts_y_hat;
    array[numObs_y,4] int  norm_f_total_y_hat;
    array[numObs_y] real  f_yfphi_hat;
    array[numObs_y] real  ki67inyfphi_hat;
    array[numObs_y] real  ki67inyfplo_hat;
        
    
    for (i in 1:numObs_y) {
        log_lik_y[i] = normal_lpdf(logit(yfphi[i]) | logit(f_yfphi[i]), sigma_yfphi)+normal_lpdf(logit(ki67inyfphi[i]) | logit(f_kihi_yfphi[i]),sigma_khiyfphi)+normal_lpdf(logit(ki67inyfplo[i]) | logit(f_kihi_yfplo[i]),sigma_khiyfplo);
        norm_f_total_y_hat[i,:] = dirichlet_multinomial_rng(norm_c_total_y[i,:]/phi_inv_y,counts_t0_y); 
        data_counts_y_hat[i] = lognormal_rng(log(c_total_y[i]),scale_y);
        
        f_yfphi_hat[i] = inv_logit(normal_rng(logit(f_yfphi[i]),sigma_yfphi));
        ki67inyfphi_hat[i] = inv_logit(normal_rng(logit(ki67inyfphi[i]),sigma_khiyfphi));
        ki67inyfplo_hat[i] = inv_logit(normal_rng(logit(ki67inyfplo[i]),sigma_khiyfplo));
        
        log_lik[i] = log_lik_y[i];
        }
        
        
   /*
     array[100] vector[switch1] k_hat_y_calc = solveode_y(init_condit_eff2, fulltime_y, time_index_equilibrium, 
      alpha_A,
         delta_A,
         beta,
         Source,
         numofki67int,
         switch1,
         100,
         a_y,
         b_y,
         kf); */
         
    array[100] vector[switch1] k_hat_y_calc = solveode_y(init_conditions2[1], fulltime_y, time_index_equilibrium, 
     alpha_A,
        delta_A,
        beta,
        Source,
        numofki67int,
        switch1,
        100,
        a_y,
        b_y,
        kf,
        timestart,
        a_y_data); 
            vector<lower = 0>[100]   c_total_y_calc;        // Total cell counts
            vector<lower = 0,upper=1>[100]  f_yfphi_kihi_calc;   // fraction of BRDuHi ki67Hi
            vector<lower = 0,upper=1>[100]  f_yfphi_kilo_calc;   // fraction of BRDuHi ki67Lo 
            vector<lower = 0,upper=1>[100]  f_yfphi_calc;        // fraction of BRDUHi
            vector<lower = 0,upper=1>[100]  f_kihi_y_calc;         // fraction of ki67Hi 
            vector<lower = 0,upper=1>[100]  f_kihi_yfphi_calc;
            vector<lower = 0,upper=1>[100]  f_kihi_yfplo_calc;
              
            vector<lower = 0>[100]  c_kihi_y_calc;        // fraction of ki67Hi 
            vector<lower = 0>[100]  c_kilo_y_calc;        // fraction of ki67Hi 
            vector<lower = 0>[100]  c_yfplo_calc;       // fraction of ki67Hi 
            vector<lower = 0>[100]  c_yfphi_kihi_calc;  // fraction of ki67Hi
            vector<lower = 0>[100]  c_yfphi_kilo_calc;  // fraction of ki67Hi
            vector<lower = 0>[100]  c_yfplo_kilo_calc;  // fraction of ki67Hi  
            vector<lower = 0>[100]  c_yfplo_kihi_calc;
            
            vector<lower = 0>[100]  f_pre_yfphi_calc;       // fraction of ki67Hi 
            
            array[100] vector[4] norm_c_total_y_calc;    // Simplex to store the normalized quadrant counts
            //vector[100] clonal;
            //array[100] vector[2] lamdas = ode_rk45(lamda,init_lamdas, 0.0, linspaced_array(100, 1, 1000), alpha_A, alpha_B, delta_A, delta_B, gamma);

            for (i in 1:100){
                // calculating
                c_yfplo_kilo_calc[i] =  k_hat_y_calc[i,switch1];
                
                c_kihi_y_calc[i] = sum(k_hat_y_calc[i,1:((2)*(1+numofki67int))]);
                c_yfplo_calc[i] = intplus(1, k_hat_y_calc[i,:], numofki67int,  0, 1);
                
                // total counts
                c_total_y_calc[i] = sum(k_hat_y_calc[i, 1:switch1]);
                

                c_kilo_y_calc[i] = c_total_y_calc[i] - c_kihi_y_calc[i];
                c_yfphi_kihi_calc[i] = c_total_y_calc[i] - c_yfplo_calc[i] - sum(k_hat_y_calc[i,((2+0)*(1+numofki67int)+1):(switch1-1)]);
                
                // fractions of yfphi cells in ki67Hi
                f_yfphi_kihi_calc[i] = (c_yfphi_kihi_calc[i])/c_kihi_y_calc[i];
              
                // fractions of yfphi cells in KI67lo
                c_yfphi_kilo_calc[i] = (sum(k_hat_y_calc[i,(1+(2)*(1+numofki67int)):(switch1-1)]));
                f_yfphi_kilo_calc[i] = c_yfphi_kilo_calc[i]/c_kilo_y_calc[i];
        
                c_yfplo_kihi_calc[i] = (c_total_y_calc[i]-c_yfphi_kilo_calc[i]-c_yfphi_kihi_calc[i]-c_yfplo_kilo_calc[i]);
                
                // fractions of ki67hi cells in yfpHi
                f_kihi_yfphi_calc[i] = (c_yfphi_kihi_calc[i])/(c_yfphi_kihi_calc[i]+c_yfphi_kilo_calc[i]);
              
                // fractions of ki67hi cells in yfplo
                f_kihi_yfplo_calc[i] = c_yfplo_kihi_calc[i]/(c_yfplo_kihi_calc[i]+c_yfplo_kilo_calc[i]);
                     
                // fractions of yfpUHi
                f_yfphi_calc[i] = (c_total_y_calc[i]-c_yfplo_calc[i])/c_total_y_calc[i];
                
                // fractions of ki67Hi
                f_kihi_y_calc[i] = c_kihi_y_calc[i]/c_total_y_calc[i];
                
                norm_c_total_y_calc[i,1] = c_yfphi_kilo_calc[i]/c_total_y_calc[i];
                norm_c_total_y_calc[i,2] = c_yfphi_kihi_calc[i]/c_total_y_calc[i];
                norm_c_total_y_calc[i,3] = c_yfplo_kilo_calc[i]/c_total_y_calc[i];
                norm_c_total_y_calc[i,4] = (c_total_y_calc[i]-c_yfphi_kilo_calc[i]-c_yfphi_kihi_calc[i]-c_yfplo_kilo_calc[i])/c_total_y_calc[i];  
                    
                //clonal[i] = (lamdas[i,1]+ lamdas[i,2])/1000; 
                
               if (a_y_data==0.372){
                   f_pre_yfphi_calc[i]=expincsol(fulltime_y[i], a_y,b_y);}
               else{
                   f_pre_yfphi_calc[i]=expdecaysol(fulltime_y[i], a_y,b_y);}
   
            }

    vector[numObs_t] log_lik_t;
    vector[numObs_t] data_counts_t_hat;
    array[numObs_t,2] int  norm_f_total_t_hat;
    array[numObs_t] real  f_tomhi_hat;
    
    for (i in 1:numObs_t) {
        log_lik_t[i] = normal_lpdf(logit(tomhi[i]) | logit(f_tomhi[i]), sigma_tomhi);
        f_tomhi_hat[i] = inv_logit(normal_rng(logit(f_tomhi[i]),sigma_tomhi));
        norm_f_total_t_hat[i,:] = dirichlet_multinomial_rng(norm_c_total_t[i,:]/phi_inv_t,counts_t0_t); 
        data_counts_t_hat[i] = lognormal_rng(log(c_total_t[i]),scale_t);

        log_lik[i+numObs_y]=log_lik_t[i];
        }
   
     array[100] vector[2] k_hat_t_calc = solveode_t(init_conditions_t[1], fulltime_t, time_index_equilibrium, 
      alpha_A,
         delta_A,
         Source,
         100,
         a_t,
         b_t,
         timestart); 
         
            vector<lower = 0>[100]   c_total_t_calc;        // Total cell counts
            vector<lower = 0,upper=1>[100]  f_tomhi_calc;        // fraction of BRDUHi
            vector<lower = 0>[100]  c_tomlo_calc;       // fraction of ki67Hi 
            vector<lower = 0>[100]  f_pre_tomhi_calc;       // fraction of ki67Hi 
            
            array[100] vector[2] norm_c_total_t_calc;    // Simplex to store the normalized quadrant counts
            //vector[100] clonal;
            //array[100] vector[2] lamdas = ode_rk45(lamda,init_lamdas, 0.0, linspaced_array(100, 1, 1000), alpha_A, alpha_B, delta_A, delta_B, gamma);

            for (i in 1:100){
                // calculating
                c_tomlo_calc[i] = k_hat_t_calc[i,2];
                
                // total counts
                c_total_t_calc[i] = sum(k_hat_t_calc[i, 1:2]);

                // fractions of tomUHi
                f_tomhi_calc[i] = (c_total_t_calc[i]-c_tomlo_calc[i])/c_total_t_calc[i];
                
                norm_c_total_t_calc[i,1] = c_tomlo_calc[i]/c_total_t_calc[i];
                norm_c_total_t_calc[i,2] = f_tomhi_calc[i]; 
                
                f_pre_tomhi_calc[i]=  expdecaysol(fulltime_t[i], a_t,b_t);
                //clonal[i] = (lamdas[i,1]+ lamdas[i,2])/1000; 
   
            } 
}

