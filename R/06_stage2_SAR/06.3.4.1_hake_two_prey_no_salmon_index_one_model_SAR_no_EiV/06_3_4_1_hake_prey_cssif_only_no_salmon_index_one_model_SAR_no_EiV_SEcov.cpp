
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
  DATA_VECTOR(hake_index_t);
    
  // overlap metrics from stage 1 model
  DATA_VECTOR(pianka_o_cssif_prey_field_t);
  DATA_VECTOR(pianka_o_cssif_rf_t);
  DATA_VECTOR(pianka_o_cssif_hake_t);

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
  PARAMETER(beta_hake_index);
    
  // parameters for the influence of overlap with prey/predators
  PARAMETER(beta_ov_prey_field);
  PARAMETER(beta_ov_rf);
  PARAMETER(beta_ov_hake);

  // Random effects
  // no random effects
  
  
  // Objective function
  Type jnll = 0;
  
  
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
      beta_prey_field_index*prey_field_index_t(t_i(i)) +
      beta_ov_prey_field*pianka_o_cssif_prey_field_t(t_i(i)) +
      beta_rf_index*rf_index_t(t_i(i)) +
      beta_ov_rf*pianka_o_cssif_rf_t(t_i(i)) +
      beta_hake_index*hake_index_t(t_i(i)) +
      beta_ov_hake*pianka_o_cssif_hake_t(t_i(i)),
      
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
  REPORT( beta_hake_index );
  
  REPORT( beta_ov_prey_field );
  REPORT( beta_ov_rf );
  REPORT( beta_ov_hake );
  
  return jnll;
}