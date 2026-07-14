// Using the best-fit SDMs for both CSSIF, CSYIF, and hake, estimate the overlap
// between juvenile salmonids and hake


#include <TMB.hpp>
#include <algorithm>
#include "utils.h"

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;

// This model uses a lot of subscripts to keep our taxa straight


// Objective function
Type jnll = 0;

// #### SHARED DATA ####

DATA_INTEGER(n_t_jsoes); // number of years in the jsoes bongo/trawl datasets
DATA_INTEGER(n_t_shared); // number of shared years between all datasets
DATA_INTEGER(n_t_hake); // number of years in the hake dataset
DATA_IVECTOR(t_jsoes_indices); // indices of the jsoes years in the shared years
DATA_IVECTOR(t_hake_indices); // indices of the hake years in the shared years

DATA_MATRIX(temp_gt); // temperature at each location in each year
DATA_VECTOR(dist_g); // distance from shore at each location


// #### CSYIF SDM ####

// Data for csyif
DATA_VECTOR(D_i_csyif);  // density of csyif in measurement i
DATA_IVECTOR(t_i_csyif); // index for the year of measurement i (csyif)
DATA_VECTOR(weights_i_csyif); // optional weights for csyif - used for cAIC
DATA_VECTOR(temp_i_csyif);  // temperature for measurement i


// Projection matrices for csyif
DATA_SPARSE_MATRIX(A_is_csyif); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_csyif); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for csyif
PARAMETER( beta_temp_csyif );
PARAMETER( ln_tau_omega_csyif ); // Tau parameter for spatial effects (omega) for csyif
PARAMETER( ln_tau_epsilon_csyif ); // Tau parameter for spatiotemporal effects (epsilon) for csyif
PARAMETER( ln_kappa_csyif ); // Kappa term for SPDE for csyif
PARAMETER( ln_phi_csyif ); // phi term in tweedie for csyif
PARAMETER( finv_power_csyif ); // power parameter in tweedie for csyif

PARAMETER( logit_rhoE_csyif ); // AR1 parameter

PARAMETER_VECTOR( ln_H_input_csyif ); // anisotropy input for csyif

Type rhoE_csyif = invlogit( logit_rhoE_csyif ); // tranform AR1 parameter

// Random effects for csyif
PARAMETER_VECTOR( omega_s_csyif ); // vector of spatial random effects for csyif
PARAMETER_MATRIX( epsilon_st_csyif ); // vector of spatiotemporal random effects for csyif

int n_i_csyif = A_is_csyif.rows(); // number of csyif samples
int n_g_csyif = A_gs_csyif.rows(); // number of units in the projection grid for csyif

// priors for anisotropy parameters for csyif
Type ln_H_0_csyif_mean = 0.0;
Type ln_H_0_csyif_sd = 0.5;
Type ln_H_1_csyif_mean = 0.0;
Type ln_H_1_csyif_sd = 0.5;

jnll -= dnorm(ln_H_input_csyif(0), ln_H_0_csyif_mean, ln_H_0_csyif_sd, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_csyif(1), ln_H_1_csyif_mean, ln_H_1_csyif_sd, true); // Anisotropic correlation

// Anisotropy elements for csyif
matrix<Type> H_csyif( 2, 2 );
H_csyif(0,0) = exp(ln_H_input_csyif(0));
H_csyif(1,0) = ln_H_input_csyif(1);
H_csyif(0,1) = ln_H_input_csyif(1);
H_csyif(1,1) = (1+ln_H_input_csyif(1)*ln_H_input_csyif(1)) / exp(ln_H_input_csyif(0));

// implement anisotropy for csyif
Eigen::SparseMatrix<Type> Q_csyif;
// Using INLA
DATA_STRUCT( spatial_list_csyif, R_inla::spde_aniso_t );
// Build precision matrix for csyif
Q_csyif = R_inla::Q_spde( spatial_list_csyif, exp(ln_kappa_csyif), H_csyif );

// Derived quantities
// SDM derived quantities for csyif
Type Range_csyif = sqrt(8) / exp( ln_kappa_csyif );
Type SigmaO_csyif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_csyif) * exp(2*ln_kappa_csyif));
Type SigmaE_csyif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_csyif) * exp(2*ln_kappa_csyif));

// Probability of random effects for csyif
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_csyif), 1/exp(ln_tau_omega_csyif) )( omega_s_csyif );
// spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  if( t==0 ){
    jnll += SCALE( GMRF(Q_csyif), 1 / exp(ln_tau_epsilon_csyif) / pow( 1.0-pow(rhoE_csyif,2), 0.5 ) )( epsilon_st_csyif.col(t) );
  }else{
    jnll += SCALE( GMRF(Q_csyif), 1 / exp(ln_tau_epsilon_csyif) )( epsilon_st_csyif.col(t) - rhoE_csyif*epsilon_st_csyif.col(t-1) );
  }
}

// SDM projections for csyif
vector<Type> omega_g_csyif( n_g_csyif );
omega_g_csyif = A_gs_csyif * omega_s_csyif;
matrix<Type> epsilon_gt_csyif( n_g_csyif, n_t_jsoes );
epsilon_gt_csyif = A_gs_csyif * epsilon_st_csyif;

// Probability of data conditional on random effects for csyif
vector<Type> omega_i_csyif( n_i_csyif );
omega_i_csyif = A_is_csyif * omega_s_csyif;
matrix<Type> epsilon_it_csyif( n_i_csyif, n_t_jsoes );
epsilon_it_csyif = A_is_csyif * epsilon_st_csyif;

vector<Type> dhat_i_csyif( n_i_csyif );
dhat_i_csyif.setZero();


for( int i=0; i<D_i_csyif.size(); i++){
  dhat_i_csyif(i) = exp(beta_temp_csyif*temp_i_csyif(i) + omega_i_csyif(i) + epsilon_it_csyif(i, t_i_csyif(i)));
  jnll -= dtweedie( D_i_csyif(i), dhat_i_csyif(i), exp(ln_phi_csyif), Type(1.0)+invlogit(finv_power_csyif), true ) * weights_i_csyif(i);
}



// SDM: modeled density in the projection grid for csyif
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_csyif( n_g_csyif, n_t_jsoes );
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    ln_d_gt_csyif(g,t) = beta_temp_csyif*temp_gt(g,t) + omega_g_csyif(g) + epsilon_gt_csyif(g,t);
  }
}


// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_csyif( n_g_csyif, n_t_jsoes );

// loop through the array to convert to normal space
for (int i = 0; i < d_gt_csyif.rows(); i++) {
  for (int j = 0; j < d_gt_csyif.cols(); j++) {
    d_gt_csyif(i,j) = exp(ln_d_gt_csyif(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> csyif_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  csyif_index_of_abundance(t) = d_gt_csyif.col(t).sum();
}

REPORT(csyif_index_of_abundance);
ADREPORT(csyif_index_of_abundance);


// Reporting: csyif
REPORT( dhat_i_csyif );
REPORT( beta_temp_csyif );
REPORT( ln_tau_omega_csyif );
REPORT( ln_tau_epsilon_csyif );
REPORT( ln_kappa_csyif );
REPORT( ln_phi_csyif );
REPORT( finv_power_csyif );
REPORT(SigmaO_csyif);
REPORT(SigmaE_csyif);
REPORT(Range_csyif);
REPORT(ln_d_gt_csyif);
REPORT( H_csyif );
REPORT( logit_rhoE_csyif );

// #### CSSIF SDM ####

// Data for cssif
DATA_VECTOR(D_i_cssif);  // density of measurement i
DATA_IVECTOR(t_i_cssif); // index for the year of measurement i
DATA_VECTOR(weights_i_cssif); // optional weights - used for cAIC
DATA_VECTOR(temp_i_cssif);  // temperature for measurement i
DATA_VECTOR(dist_i_cssif);  // distance from shore for measurement i

// Data for smoothers for cssif
DATA_IVECTOR(b_smooth_start_cssif);
DATA_STRUCT(Zs_cssif, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
DATA_STRUCT(proj_Zs_cssif, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
DATA_MATRIX(Xs_cssif); // smoother linear effect matrix
DATA_MATRIX(proj_Xs_cssif); // smoother linear effect matrix

// Projection matrices for cssif
// A_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
DATA_SPARSE_MATRIX(A_is_cssif);
// A_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
DATA_SPARSE_MATRIX(A_gs_cssif);


// Parameters for cssif

// Parameters for SDM for cssif
PARAMETER( beta_temp_cssif );
PARAMETER( ln_tau_omega_cssif );
PARAMETER( ln_tau_epsilon_cssif );
PARAMETER( ln_kappa_cssif );
PARAMETER( ln_phi_cssif ); // phi term in tweedie
PARAMETER( finv_power_cssif ); // power parameter in tweedie

// AR1 parameters
// PARAMETER( logit_rhoE );

PARAMETER_VECTOR( ln_H_input_cssif ); // anisotropy input.

// smoother parameters
PARAMETER_VECTOR(bs_cssif); // smoother linear effects
PARAMETER_VECTOR(b_smooth_cssif);  // P-spline smooth parameters
PARAMETER_VECTOR(ln_smooth_sigma_cssif);  // variances of spline REs if included

// Random effects for cssif

PARAMETER_VECTOR( omega_s_cssif );
PARAMETER_MATRIX( epsilon_st_cssif );

// Objective function
int n_i_cssif = A_is_cssif.rows();
int n_g_cssif = A_gs_cssif.rows();

// p-splines/smoothers
vector<Type> eta_smooth_i_cssif(n_i_cssif);
eta_smooth_i_cssif.setZero();
for (int s = 0; s < b_smooth_start_cssif.size(); s++) { // iterate over # of smooth elements
  vector<Type> beta_s_cssif(Zs_cssif(s).cols());
  beta_s_cssif.setZero();
  for (int j = 0; j < beta_s_cssif.size(); j++) {
    beta_s_cssif(j) = b_smooth_cssif(b_smooth_start_cssif(s) + j);
    // PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
    jnll -= dnorm(beta_s_cssif(j), Type(0), exp(ln_smooth_sigma_cssif(s)), true);
  }
  eta_smooth_i_cssif += Zs_cssif(s) * beta_s_cssif;
}
// eta_smooth_i += Xs * vector<Type>(bs.col(1));
eta_smooth_i_cssif += Xs_cssif * bs_cssif;
REPORT(b_smooth_cssif);     // smooth coefficients for penalized splines
REPORT(ln_smooth_sigma_cssif); // standard deviations of smooth random effects, in log-space

// AR1 parameters
// Type rhoE = invlogit( logit_rhoE );

// priors for anisotropy parameters
// turned on in this model
Type ln_H_0_mean_cssif = 0.0;
Type ln_H_0_sd_cssif = 0.5;
Type ln_H_1_mean_cssif = 0.0;
Type ln_H_1_sd_cssif = 0.5;

jnll -= dnorm(ln_H_input_cssif(0), ln_H_0_mean_cssif, ln_H_0_sd_cssif, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_cssif(1), ln_H_1_mean_cssif, ln_H_1_sd_cssif, true); // Anisotropic correlation

// Anisotropy elements
matrix<Type> H_cssif( 2, 2 );
H_cssif(0,0) = exp(ln_H_input_cssif(0));
H_cssif(1,0) = ln_H_input_cssif(1);
H_cssif(0,1) = ln_H_input_cssif(1);
H_cssif(1,1) = (1+ln_H_input_cssif(1)*ln_H_input_cssif(1)) / exp(ln_H_input_cssif(0));
// H(0,0) = 1;
// H(0,1) = 0;
// H(1,0) = 0;
// H(1,1) = 1;

// implement anisotropy
Eigen::SparseMatrix<Type> Q_cssif;
// Using INLA
DATA_STRUCT( spatial_list_cssif, R_inla::spde_aniso_t );
// Build precision
Q_cssif = R_inla::Q_spde( spatial_list_cssif, exp(ln_kappa_cssif), H_cssif );

// Derived quantities
// SDM derived quantities for species A
Type Range_cssif = sqrt(8) / exp( ln_kappa_cssif );
Type SigmaO_cssif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_cssif) * exp(2*ln_kappa_cssif));
Type SigmaE_cssif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_cssif) * exp(2*ln_kappa_cssif));


// Probability of random effects

// SDM:
// Eigen::SparseMatrix<Type> Q_a = exp(4*ln_kappa_a)*M0_a + Type(2.0)*exp(2*ln_kappa_a)*M1_a + M2_a;
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_cssif), 1/exp(ln_tau_omega_cssif) )( omega_s_cssif );
// spatio-temporal random effect - scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  jnll += SCALE( GMRF(Q_cssif), 1/exp(ln_tau_epsilon_cssif) )( epsilon_st_cssif.col(t) );
}


// SDM projections
vector<Type> omega_g_cssif( n_g_cssif );
omega_g_cssif = A_gs_cssif * omega_s_cssif;
matrix<Type> epsilon_gt_cssif( n_g_cssif, n_t_jsoes );
epsilon_gt_cssif = A_gs_cssif * epsilon_st_cssif;

// Probability of data conditional on random effects

// SDM species A
vector<Type> omega_i_cssif( n_i_cssif );
omega_i_cssif = A_is_cssif * omega_s_cssif;
matrix<Type> epsilon_it_cssif( n_i_cssif, n_t_jsoes );
epsilon_it_cssif = A_is_cssif * epsilon_st_cssif;

vector<Type> dhat_i_cssif( n_i_cssif );
dhat_i_cssif.setZero();


for( int i=0; i<D_i_cssif.size(); i++){
  dhat_i_cssif(i) = exp( epsilon_it_cssif(i, t_i_cssif(i)) + omega_i_cssif(i) +
    beta_temp_cssif*temp_i_cssif(i) + // beta_dist*dist_i(i));
    eta_smooth_i_cssif(i));
  jnll -= dtweedie( D_i_cssif(i), dhat_i_cssif(i), exp(ln_phi_cssif), Type(1.0)+invlogit(finv_power_cssif), true ) * weights_i_cssif(i);
}



// SDM: modeled density in the projection grid for species A

// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_cssif( n_g_cssif, n_t_jsoes );

// estimate the smooth effects in the projection grid
// smoothers for projection grid
vector<Type> proj_smooth_i_cssif(n_g_cssif*n_t_jsoes);
proj_smooth_i_cssif.setZero();
for (int s = 0; s < b_smooth_start_cssif.size(); s++) { // iterate over # of smooth elements
  vector<Type> beta_s_cssif(proj_Zs_cssif(s).cols());
  beta_s_cssif.setZero();
  for (int j = 0; j < beta_s_cssif.size(); j++) {
    beta_s_cssif(j) = b_smooth_cssif(b_smooth_start_cssif(s) + j);
  }
  proj_smooth_i_cssif += proj_Zs_cssif(s) * beta_s_cssif;
}
proj_smooth_i_cssif += proj_Xs_cssif * bs_cssif;

// add smoothed effects, fixed effects, and random effects
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_cssif; g++){
    // ln_d_gt(g,t) = beta_t(t) + beta_temp*temp_gt(g,t) + beta_dist*dist_g(g) + omega_g(g) + epsilon_gt(g,t) ;
    ln_d_gt_cssif(g,t) = omega_g_cssif(g) + epsilon_gt_cssif(g,t) + beta_temp_cssif*temp_gt(g,t) + proj_smooth_i_cssif(g + n_g_cssif * t) ;
  }
}

// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_cssif( n_g_cssif, n_t_jsoes );

// loop through the array to convert to normal space
for (int i = 0; i < d_gt_cssif.rows(); i++) {
  for (int j = 0; j < d_gt_cssif.cols(); j++) {
    d_gt_cssif(i,j) = exp(ln_d_gt_cssif(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> cssif_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  cssif_index_of_abundance(t) = d_gt_cssif.col(t).sum();
}

REPORT(cssif_index_of_abundance);
ADREPORT(cssif_index_of_abundance);



// Reporting: SDM
REPORT( dhat_i_cssif );
// REPORT( beta_t );
REPORT( beta_temp_cssif );
// REPORT( beta_dist );
REPORT( ln_tau_omega_cssif );
REPORT( ln_tau_epsilon_cssif );
REPORT( ln_kappa_cssif );
REPORT( ln_phi_cssif );
REPORT( finv_power_cssif );
REPORT(SigmaO_cssif);
REPORT(SigmaE_cssif);
REPORT(Range_cssif);
REPORT(ln_d_gt_cssif);
REPORT( H_cssif );
// REPORT( logit_rhoE );


// #### HAKE ####

  // Data
  DATA_VECTOR(D_i_hake);  // density of measurement i
  DATA_IVECTOR(t_i_hake); // index for the year of measurement i
  DATA_VECTOR(weights_i_hake); // optional weights - used for cAIC
  DATA_VECTOR(temp_i_hake);  // temperature for measurement i
  DATA_VECTOR(dist_i_hake);  // distance from shore for measurement i
  
  // Data for smoothers
  DATA_IVECTOR(b_smooth_start_hake);
  DATA_STRUCT(Zs_hake, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(proj_Zs_hake, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs_hake); // smoother linear effect matrix
  DATA_MATRIX(proj_Xs_hake); // smoother linear effect matrix
  
  // Projection matrices for species A
  // A_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
  DATA_SPARSE_MATRIX(A_is_hake);
  // A_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
  DATA_SPARSE_MATRIX(A_gs_hake);
  
  
  // Parameters
  
  // Parameters for SDM for species A
  PARAMETER( ln_tau_epsilon_hake );
  PARAMETER( ln_kappa_hake );
  PARAMETER( ln_phi_hake ); // phi term in tweedie
  PARAMETER( finv_power_hake ); // power parameter in tweedie
  
  // AR1 parameters
  PARAMETER( logit_rhoE_hake );
  
  // PARAMETER_VECTOR( ln_H_input ); // anisotropy input.
  
  // smoother parameters
  PARAMETER_VECTOR(bs_hake); // smoother linear effects
  PARAMETER_VECTOR(b_smooth_hake);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma_hake);  // variances of spline REs if included
  
  // Random effects
  
  // PARAMETER_VECTOR( omega_s );
  PARAMETER_MATRIX( epsilon_st_hake );
  
  int n_i_hake = A_is_hake.rows();
  int n_g_hake = A_gs_hake.rows();
  
  // p-splines/smoothers
  vector<Type> eta_smooth_i_hake(n_i_hake);
  eta_smooth_i_hake.setZero();
  for (int s = 0; s < b_smooth_start_hake.size(); s++) { // iterate over # of smooth elements
    vector<Type> beta_s_hake(Zs_hake(s).cols());
    beta_s_hake.setZero();
    for (int j = 0; j < beta_s_hake.size(); j++) {
      beta_s_hake(j) = b_smooth_hake(b_smooth_start_hake(s) + j);
      // PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
      jnll -= dnorm(beta_s_hake(j), Type(0), exp(ln_smooth_sigma_hake(s)), true);
    }
    eta_smooth_i_hake += Zs_hake(s) * beta_s_hake;
  }
  // eta_smooth_i += Xs * vector<Type>(bs.col(1));
  eta_smooth_i_hake += Xs_hake * bs_hake;
  REPORT(b_smooth_hake);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma_hake); // standard deviations of smooth random effects, in log-space
  
  // AR1 parameters
  Type rhoE_hake = invlogit( logit_rhoE_hake );
  
  // priors for anisotropy parameters
  // turned OFF in this model
  // Type ln_H_0_mean = 0.0;
  // Type ln_H_0_sd = 0.5;
  // Type ln_H_1_mean = 0.0;
  // Type ln_H_1_sd = 0.5;
  
  // jnll -= dnorm(ln_H_input(0), ln_H_0_mean, ln_H_0_sd, true); // Northings anisotropy
  // jnll -= dnorm(ln_H_input(1), ln_H_1_mean, ln_H_1_sd, true); // Anisotropic correlation
  
  // Anisotropy elements
  matrix<Type> H_hake( 2, 2 );
  // H(0,0) = exp(ln_H_input(0));
  // H(1,0) = ln_H_input(1);
  // H(0,1) = ln_H_input(1);
  // H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  H_hake(0,0) = 1;
  H_hake(0,1) = 0;
  H_hake(1,0) = 0;
  H_hake(1,1) = 1;
  
  // implement anisotropy
  Eigen::SparseMatrix<Type> Q_hake;
  // Using INLA
  DATA_STRUCT( spatial_list_hake, R_inla::spde_aniso_t );
  // Build precision
  Q_hake = R_inla::Q_spde( spatial_list_hake, exp(ln_kappa_hake), H_hake );

  // Derived quantities
  // SDM derived quantities for species A
  Type Range_hake = sqrt(8) / exp( ln_kappa_hake );
  // Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega) * exp(2*ln_kappa));
  Type SigmaE_hake = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_hake) * exp(2*ln_kappa_hake));
  
  
  // Probability of random effects

  // SDM:
  // Eigen::SparseMatrix<Type> Q_a = exp(4*ln_kappa_a)*M0_a + Type(2.0)*exp(2*ln_kappa_a)*M1_a + M2_a;
  // spatial random effect - scaled by ln_tau_omega
  // jnll += SCALE( GMRF(Q), 1/exp(ln_tau_omega) )( omega_s );
  // spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
  for( int t=0; t<n_t_hake; t++){
    if( t==0 ){
      jnll += SCALE( GMRF(Q_hake), 1 / exp(ln_tau_epsilon_hake) / pow( 1.0-pow(rhoE_hake,2), 0.5 ) )( epsilon_st_hake.col(t) );
    }else{
      jnll += SCALE( GMRF(Q_hake), 1 / exp(ln_tau_epsilon_hake) )( epsilon_st_hake.col(t) - rhoE_hake*epsilon_st_hake.col(t-1) );
    }
  }


  
  // SDM projections
  // vector<Type> omega_g( n_g );
  // omega_g = A_gs * omega_s;
  matrix<Type> epsilon_gt_hake( n_g_hake, n_t_hake );
  epsilon_gt_hake = A_gs_hake * epsilon_st_hake;
  
  // Probability of data conditional on random effects
  
  // SDM species A
  // vector<Type> omega_i( n_i );
  // omega_i = A_is * omega_s;
  matrix<Type> epsilon_it_hake( n_i_hake, n_t_hake );
  epsilon_it_hake = A_is_hake * epsilon_st_hake;

  vector<Type> dhat_i_hake( n_i_hake );
  dhat_i_hake.setZero();


  for( int i=0; i<D_i_hake.size(); i++){
    dhat_i_hake(i) = exp( epsilon_it_hake(i, t_i_hake(i))  +
      // beta_temp*temp_i(i) + beta_dist*dist_i(i));
      eta_smooth_i_hake(i));
    jnll -= dtweedie( D_i_hake(i), dhat_i_hake(i), exp(ln_phi_hake), Type(1.0)+invlogit(finv_power_hake), true ) * weights_i_hake(i);
  }
  
  
  
  // SDM: modeled density in the projection grid for species A
  
  // Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
  array<Type> ln_d_gt_hake( n_g_hake, n_t_hake );
  
  // estimate the smooth effects in the projection grid
  // smoothers for projection grid
  vector<Type> proj_smooth_i_hake(n_g_hake*n_t_hake);
  proj_smooth_i_hake.setZero();
      for (int s = 0; s < b_smooth_start_hake.size(); s++) { // iterate over # of smooth elements
        vector<Type> beta_s_hake(proj_Zs_hake(s).cols());
        beta_s_hake.setZero();
        for (int j = 0; j < beta_s_hake.size(); j++) {
          beta_s_hake(j) = b_smooth_hake(b_smooth_start_hake(s) + j);
        }
        proj_smooth_i_hake += proj_Zs_hake(s) * beta_s_hake;
      }
      proj_smooth_i_hake += proj_Xs_hake * bs_hake;
  
  // add smoothed effects, fixed effects, and random effects
  for( int t=0; t<n_t_hake; t++){
    for( int g=0; g<n_g_hake; g++){
      // ln_d_gt(g,t) = beta_t(t) + beta_temp*temp_gt(g,t) + beta_dist*dist_g(g) + epsilon_gt(g,t) ;
      ln_d_gt_hake(g,t) = epsilon_gt_hake(g,t) + proj_smooth_i_hake(g + n_g_hake * t) ;
    }
  }
  
  // SDM: Estimate an index of abundance
  // sum across all cells, multiply by the area of each cell
  // convert to normal space
  array<Type> d_gt_hake( n_g_hake, n_t_hake );
  
  // loop through the array to convert to normal space
  for (int i = 0; i < d_gt_hake.rows(); i++) {
    for (int j = 0; j < d_gt_hake.cols(); j++) {
      d_gt_hake(i,j) = exp(ln_d_gt_hake(i,j));
    }
  }
  
  // create a vector to store the index of abundance and calculate it for each year
  vector<Type> hake_index_of_abundance(n_t_hake);
  for( int t=0; t<n_t_hake; t++){
    hake_index_of_abundance(t) = d_gt_hake.col(t).sum();
  }
  
  REPORT(hake_index_of_abundance);
  ADREPORT(hake_index_of_abundance);
  
  
  // Reporting: SDM
  REPORT( dhat_i_hake );
  REPORT( ln_tau_epsilon_hake );
  REPORT( ln_kappa_hake );
  REPORT( ln_phi_hake );
  REPORT( finv_power_hake );
  // REPORT(SigmaO);
  REPORT(SigmaE_hake);
  REPORT(Range_hake);
  REPORT(ln_d_gt_hake);
  // REPORT( H );
  REPORT( logit_rhoE_hake );



// PREDATORS

// #### CALCULATE OVERLAP WITH hake ####
// Use Pianka's O: 0-1 index that estimates correlation betweeen prey and predator taxa
// Calculate this separately for subyearlings and yearlings

// Start by defining the denominators for CSYIF and CSSIF - only need to do this once, and it gets re-used
array<Type> d_gt_csyif_sq( n_g_csyif, n_t_jsoes ); // array to store the squared proportional densities of csyif in each cell (for denominators)
array<Type> d_gt_cssif_sq( n_g_cssif, n_t_jsoes ); // array to store the squared proportional densities of csyif in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_csyif_sq(g,t) = pow(d_gt_csyif(g,t)/d_gt_csyif.col((t)).sum(),2); // calculate value for each cell
    d_gt_cssif_sq(g,t) = pow(d_gt_cssif(g,t)/d_gt_cssif.col((t)).sum(),2); // calculate value for each cell
  }
}

// Define denominator for hake - only need to do this once, and it gets re-used
array<Type> d_gt_hake_sq( n_g_csyif, n_t_hake); // array to store the squared proportional densities of hake in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_hake; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_hake_sq(g,t) = pow(d_gt_hake(g,t)/d_gt_hake.col(t).sum(),2); // calculate value for each cell
  }
}

// #### csyif and hake ####
vector<Type> pianka_o_csyif_hake_t( n_t_shared ); // vector to store the value per year
array<Type> ov_gt_num_csyif_hake( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// calculate the numerator values
for( int t=0; t<n_t_shared; t++){
  for( int g=0; g<n_g_csyif; g++){
    ov_gt_num_csyif_hake(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_hake(g,t_hake_indices(t))/d_gt_hake.col(t_hake_indices(t)).sum(); // calculate value for each cell
  }
}


// calculate the index for each year
for( int t=0; t<n_t_shared; t++){
  pianka_o_csyif_hake_t(t) = ov_gt_num_csyif_hake.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_hake_sq.col(t_hake_indices(t)).sum());
}

// Reporting: overlap for csyif and hake
REPORT(pianka_o_csyif_hake_t); // report overlap for full JSOES domain
// ADREPORT(pianka_o_csyif_hake_t); // report overlap for full JSOES domain


// #### cssif and hake ####
vector<Type> pianka_o_cssif_hake_t( n_t_shared ); // vector to store the value per year
array<Type> ov_gt_num_cssif_hake( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// calculate the numerator values
for( int t=0; t<n_t_shared; t++){
  for( int g=0; g<n_g_cssif; g++){
    ov_gt_num_cssif_hake(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_hake(g,t_hake_indices(t))/d_gt_hake.col(t_hake_indices(t)).sum(); // calculate value for each cell
  }
}


// calculate the index for each year
for( int t=0; t<n_t_shared; t++){
  pianka_o_cssif_hake_t(t) = ov_gt_num_cssif_hake.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_hake_sq.col(t_hake_indices(t)).sum());
}

// Reporting: overlap for cssif and hake
REPORT(pianka_o_cssif_hake_t); // report overlap for full JSOES domain
// ADREPORT(pianka_o_cssif_hake_t); // report overlap for full JSOES domain

return jnll;
}