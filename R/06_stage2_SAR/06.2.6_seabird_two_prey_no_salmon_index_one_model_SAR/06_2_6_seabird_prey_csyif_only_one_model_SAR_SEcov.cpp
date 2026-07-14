
#include <TMB.hpp>
#include <algorithm>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;
  // SAR:
  // Data
  DATA_VECTOR(z_i);  // survival of fish i
  DATA_IVECTOR(t_i); // index for the year of fish i
  DATA_INTEGER(n_t); // number of years in the dataset
  
  // individual-level covariates
  DATA_IVECTOR(SRF_transported_i); // 0 if fish i is not SRF transported, 1 if it is
  DATA_IVECTOR(SRF_non_transported_i); // 0 if fish i is not SRF non-transported, 1 if it is
  DATA_IVECTOR(UCSF_i); // 0 if fish is not UCSF, 1 if it is
  DATA_VECTOR(outmigration_i); // outmigration date of fish i
  
  // indices of abundance from stage 1 model
  DATA_VECTOR(prey_field_index_t);
  DATA_VECTOR(rf_index_t);
  DATA_VECTOR(sosh_index_t);
  DATA_VECTOR(comu_index_t);
  
  // uncertainty in indices of abundance from stage 1 model
  DATA_MATRIX(Sigma_prey_field_index_t);  // covariance matrix for index of abundance for prey field
  DATA_MATRIX(Sigma_rf_index_t);  // covariance matrix for index of abundance for rf
  DATA_MATRIX(Sigma_sosh_index_t);  // covariance matrix for index of abundance for prey field
  DATA_MATRIX(Sigma_comu_index_t);  // covariance matrix for index of abundance for prey field
    
  // overlap metrics from stage 1 model
  DATA_VECTOR(pianka_o_csyif_prey_field_t);
  DATA_VECTOR(pianka_o_csyif_rf_t);
  DATA_VECTOR(pianka_o_csyif_sosh_t);
  DATA_VECTOR(pianka_o_csyif_comu_t);
  
  // uncertainty in overlap metrics from stage 1 model
  // DATA_MATRIX(Sigma_pianka_o_csyif_prey_field_t);  // covariance matrix for index of abundance for csyif
  // because we can't ADReport() this index to get the covariance, here we will just treat this as a fixed measurement error.
  
  DATA_MATRIX(Sigma_pianka_o_csyif_prey_field_t); // estimated standard error for overlap metric with prey field
  DATA_MATRIX(Sigma_pianka_o_csyif_rf_t); // estimated standard error for overlap metric with rf
  DATA_MATRIX(Sigma_pianka_o_csyif_sosh_t); // estimated standard error for overlap metric with sosh
  DATA_MATRIX(Sigma_pianka_o_csyif_comu_t); // estimated standard error for overlap metric with comu

  
  // Parameters
  PARAMETER( beta_UCSF );
  PARAMETER( beta_SRF_transported );
  PARAMETER( beta_SRF_non_transported );
  
  // effect of individual-level covariates
  PARAMETER(beta_outmigration_SRF); // effect of outmigration date for SRF
  PARAMETER(beta_outmigration_SRF2); // effect of outmigration date^2 for SRF
  PARAMETER(beta_outmigration_UCSF); // effect of outmigration date for UCSF
  PARAMETER(beta_outmigration_UCSF2); // effect of outmigration date^2 for UCSF
  
  // parameters for the influence of abundance of prey/predators
  PARAMETER(beta_prey_field_index);
  PARAMETER(beta_rf_index);
  PARAMETER(beta_sosh_index);
  PARAMETER(beta_comu_index);
    
  // parameters for the influence of overlap with prey/predators
  PARAMETER(beta_ov_prey_field);
  PARAMETER(beta_ov_rf);
  PARAMETER(beta_ov_sosh);
  PARAMETER(beta_ov_comu);

  // Random effects
  // latent variables for EiV approach
  PARAMETER_VECTOR(prey_field_index_t_latent);
  PARAMETER_VECTOR(pianka_o_csyif_prey_field_t_latent);
  PARAMETER_VECTOR(rf_index_t_latent);
  PARAMETER_VECTOR(pianka_o_csyif_rf_t_latent);
  PARAMETER_VECTOR(sosh_index_t_latent);
  PARAMETER_VECTOR(pianka_o_csyif_sosh_t_latent);
  PARAMETER_VECTOR(comu_index_t_latent);
  PARAMETER_VECTOR(pianka_o_csyif_comu_t_latent);
  
  
  // Objective function
  Type jnll = 0;
  
  // Estimate actual indices of abundance as latent random effects
  
  // MVN for errors in variables of indices of abundance
  
  // prey_field index of abundance
  vector<Type> residual_prey_field_index_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_prey_field_index(Sigma_prey_field_index_t);
  residual_prey_field_index_t = vector<Type>(prey_field_index_t - prey_field_index_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_prey_field_index(residual_prey_field_index_t); // then we use those residuals and evaluate them
  
  // rf index of abundance
  vector<Type> residual_rf_index_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_rf_index(Sigma_rf_index_t);
  residual_rf_index_t = vector<Type>(rf_index_t - rf_index_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_rf_index(residual_rf_index_t); // then we use those residuals and evaluate them
  
  // sosh index of abundance
  vector<Type> residual_sosh_index_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_sosh_index(Sigma_sosh_index_t);
  residual_sosh_index_t = vector<Type>(sosh_index_t - sosh_index_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_sosh_index(residual_sosh_index_t); // then we use those residuals and evaluate them
  
  // comu index of abundance
  vector<Type> residual_comu_index_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_comu_index(Sigma_comu_index_t);
  residual_comu_index_t = vector<Type>(comu_index_t - comu_index_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_comu_index(residual_comu_index_t); // then we use those residuals and evaluate them
  
  // csyif_prey_field overlap metric
  vector<Type> residual_pianka_o_csyif_prey_field_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_pianka_o_csyif_prey_field(Sigma_pianka_o_csyif_prey_field_t);
  residual_pianka_o_csyif_prey_field_t = vector<Type>(pianka_o_csyif_prey_field_t - pianka_o_csyif_prey_field_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_pianka_o_csyif_prey_field(residual_pianka_o_csyif_prey_field_t); // then we use those residuals and evaluate them
  
  // csyif_rf overlap metric
  vector<Type> residual_pianka_o_csyif_rf_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_pianka_o_csyif_rf(Sigma_pianka_o_csyif_rf_t);
  residual_pianka_o_csyif_rf_t = vector<Type>(pianka_o_csyif_rf_t - pianka_o_csyif_rf_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_pianka_o_csyif_rf(residual_pianka_o_csyif_rf_t); // then we use those residuals and evaluate them
  
  // csyif_sosh overlap metric
  vector<Type> residual_pianka_o_csyif_sosh_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_pianka_o_csyif_sosh(Sigma_pianka_o_csyif_sosh_t);
  residual_pianka_o_csyif_sosh_t = vector<Type>(pianka_o_csyif_sosh_t - pianka_o_csyif_sosh_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_pianka_o_csyif_sosh(residual_pianka_o_csyif_sosh_t); // then we use those residuals and evaluate them
  
  // csyif_comu overlap metric
  vector<Type> residual_pianka_o_csyif_comu_t(n_t);
  MVNORM_t<Type> neg_log_dmvnorm_pianka_o_csyif_comu(Sigma_pianka_o_csyif_comu_t);
  residual_pianka_o_csyif_comu_t = vector<Type>(pianka_o_csyif_comu_t - pianka_o_csyif_comu_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  jnll += neg_log_dmvnorm_pianka_o_csyif_comu(residual_pianka_o_csyif_comu_t); // then we use those residuals and evaluate the
  
  // SAR model, using overlap and indices as covariates
  for( int i=0; i<z_i.size(); i++){
    // jnll -= dbinom_robust(z_i(i), Type(1), beta_ov_ab*ov_ab_t_latent(t_i(i)) + epsilon_t(t_i(i)), true ); // this is a Bernoulli (by using Type(1) as the size argument)
    jnll -= dbinom_robust(z_i(i), Type(1),
      // individual level covariates
      beta_UCSF * UCSF_i(i) +
        beta_SRF_transported * SRF_transported_i(i) +
        beta_SRF_non_transported * SRF_non_transported_i(i) +
        beta_outmigration_SRF * (SRF_transported_i(i) + SRF_non_transported_i(i)) * outmigration_i(i) +
        beta_outmigration_SRF2 * pow((SRF_transported_i(i) + SRF_non_transported_i(i)) * outmigration_i(i),2) +
        beta_outmigration_UCSF * (UCSF_i(i)) * outmigration_i(i) +
        beta_outmigration_UCSF2 * pow((UCSF_i(i)) * outmigration_i(i),2) +
      // SDM derived outputs
      beta_prey_field_index*prey_field_index_t_latent(t_i(i)) +
      beta_ov_prey_field*pianka_o_csyif_prey_field_t_latent(t_i(i)) +
      beta_rf_index*rf_index_t_latent(t_i(i)) +
      beta_ov_rf*pianka_o_csyif_rf_t_latent(t_i(i)) +
      beta_sosh_index*sosh_index_t_latent(t_i(i)) +
      beta_ov_sosh*pianka_o_csyif_sosh_t_latent(t_i(i)) +
      beta_comu_index*comu_index_t_latent(t_i(i)) +
      beta_ov_comu*pianka_o_csyif_comu_t_latent(t_i(i)),
      
                          true ); // this is a Bernoulli (by using Type(1) as the size argument)
  }
  
  
  // Reporting: SAR
  REPORT( beta_UCSF );
  REPORT( beta_SRF_transported );
  REPORT( beta_SRF_non_transported );
  
  // effect of individual-level covariates
  REPORT(beta_outmigration_SRF); // effect of outmigration date for SRF
  REPORT(beta_outmigration_SRF2); // effect of outmigration date^2 for SRF
  REPORT(beta_outmigration_UCSF); // effect of outmigration date for UCSF
  REPORT(beta_outmigration_UCSF2); // effect of outmigration date^2 for UCSF
  
  // effect of derived outputs
  REPORT( beta_prey_field_index );
  REPORT( beta_rf_index );
  REPORT( beta_sosh_index );
  REPORT( beta_comu_index );
  
  REPORT( beta_ov_prey_field );
  REPORT( beta_ov_rf );
  REPORT( beta_ov_sosh );
  REPORT( beta_ov_comu );
  
  // report latent variables
  REPORT(prey_field_index_t_latent);
  REPORT(pianka_o_csyif_prey_field_t_latent);
  REPORT(rf_index_t_latent);
  REPORT(pianka_o_csyif_rf_t_latent);
  REPORT(sosh_index_t_latent);
  REPORT(pianka_o_csyif_sosh_t_latent);
  REPORT(comu_index_t_latent);
  REPORT(pianka_o_csyif_comu_t_latent);
  
  return jnll;
}