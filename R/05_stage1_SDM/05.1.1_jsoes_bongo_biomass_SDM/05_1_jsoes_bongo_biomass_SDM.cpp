// Using the best-fit SDMs for both CSSIF, CSYIF and five prey items, estimate the overlap
// between juvenile salmonids and the aggregate prey field


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
DATA_IVECTOR(t_jsoes_indices); // indices of the jsoes years in the shared yaers

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


// #### Cancer Crab Larvae SDM ####

// Data for cancer_crab_larvae
DATA_VECTOR(D_i_cancer_crab_larvae);  // density of cancer_crab_larvae in measurement i
DATA_IVECTOR(t_i_cancer_crab_larvae); // index for the year of measurement i (cancer_crab_larvae)
DATA_VECTOR(weights_i_cancer_crab_larvae); // optional weights for cancer_crab_larvae - used for cAIC


// Projection matrices for cancer_crab_larvae
DATA_SPARSE_MATRIX(A_is_cancer_crab_larvae); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_cancer_crab_larvae); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for cancer_crab_larvae
PARAMETER( ln_tau_omega_cancer_crab_larvae ); // Tau parameter for spatial effects (omega) for cancer_crab_larvae
PARAMETER( ln_tau_epsilon_cancer_crab_larvae ); // Tau parameter for spatiotemporal effects (epsilon) for cancer_crab_larvae
PARAMETER( ln_kappa_cancer_crab_larvae ); // Kappa term for SPDE for cancer_crab_larvae
PARAMETER( ln_phi_cancer_crab_larvae ); // phi term in tweedie for cancer_crab_larvae
PARAMETER( finv_power_cancer_crab_larvae ); // power parameter in tweedie for cancer_crab_larvae
PARAMETER( logit_rhoE_cancer_crab_larvae ); // AR1 parameter

PARAMETER_VECTOR( ln_H_input_cancer_crab_larvae ); // anisotropy input for cancer_crab_larvae

Type rhoE_cancer_crab_larvae = invlogit( logit_rhoE_cancer_crab_larvae ); // tranform AR1 parameter

// Random effects for cancer_crab_larvae
PARAMETER_VECTOR( omega_s_cancer_crab_larvae ); // vector of spatial random effects for cancer_crab_larvae
PARAMETER_MATRIX( epsilon_st_cancer_crab_larvae ); // vector of spatiotemporal random effects for cancer_crab_larvae

int n_i_cancer_crab_larvae = A_is_cancer_crab_larvae.rows(); // number of cancer_crab_larvae samples
int n_g_cancer_crab_larvae = A_gs_cancer_crab_larvae.rows(); // number of units in the projection grid for cancer_crab_larvae

// priors for anisotropy parameters for cancer_crab_larvae
Type ln_H_0_cancer_crab_larvae_mean = 0.0;
Type ln_H_0_cancer_crab_larvae_sd = 0.5;
Type ln_H_1_cancer_crab_larvae_mean = 0.0;
Type ln_H_1_cancer_crab_larvae_sd = 0.5;

jnll -= dnorm(ln_H_input_cancer_crab_larvae(0), ln_H_0_cancer_crab_larvae_mean, ln_H_0_cancer_crab_larvae_sd, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_cancer_crab_larvae(1), ln_H_1_cancer_crab_larvae_mean, ln_H_1_cancer_crab_larvae_sd, true); // Anisotropic correlation

// Anisotropy elements for cancer_crab_larvae
matrix<Type> H_cancer_crab_larvae( 2, 2 );
H_cancer_crab_larvae(0,0) = exp(ln_H_input_cancer_crab_larvae(0));
H_cancer_crab_larvae(1,0) = ln_H_input_cancer_crab_larvae(1);
H_cancer_crab_larvae(0,1) = ln_H_input_cancer_crab_larvae(1);
H_cancer_crab_larvae(1,1) = (1+ln_H_input_cancer_crab_larvae(1)*ln_H_input_cancer_crab_larvae(1)) / exp(ln_H_input_cancer_crab_larvae(0));

// implement anisotropy for cancer_crab_larvae
Eigen::SparseMatrix<Type> Q_cancer_crab_larvae;
// Using INLA
DATA_STRUCT( spatial_list_cancer_crab_larvae, R_inla::spde_aniso_t );
// Build precision matrix for cancer_crab_larvae
Q_cancer_crab_larvae = R_inla::Q_spde( spatial_list_cancer_crab_larvae, exp(ln_kappa_cancer_crab_larvae), H_cancer_crab_larvae );

// Derived quantities
// SDM derived quantities for cancer_crab_larvae
Type Range_cancer_crab_larvae = sqrt(8) / exp( ln_kappa_cancer_crab_larvae );
Type SigmaO_cancer_crab_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_cancer_crab_larvae) * exp(2*ln_kappa_cancer_crab_larvae));
Type SigmaE_cancer_crab_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_cancer_crab_larvae) * exp(2*ln_kappa_cancer_crab_larvae));

// Probability of random effects for cancer_crab_larvae
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_cancer_crab_larvae), 1/exp(ln_tau_omega_cancer_crab_larvae) )( omega_s_cancer_crab_larvae );
// spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  if( t==0 ){
    jnll += SCALE( GMRF(Q_cancer_crab_larvae), 1 / exp(ln_tau_epsilon_cancer_crab_larvae) / pow( 1.0-pow(rhoE_cancer_crab_larvae,2), 0.5 ) )( epsilon_st_cancer_crab_larvae.col(t) );
  }else{
    jnll += SCALE( GMRF(Q_cancer_crab_larvae), 1 / exp(ln_tau_epsilon_cancer_crab_larvae) )( epsilon_st_cancer_crab_larvae.col(t) - rhoE_cancer_crab_larvae*epsilon_st_cancer_crab_larvae.col(t-1) );
  }
}

// SDM projections for cancer_crab_larvae
vector<Type> omega_g_cancer_crab_larvae( n_g_cancer_crab_larvae );
omega_g_cancer_crab_larvae = A_gs_cancer_crab_larvae * omega_s_cancer_crab_larvae;
matrix<Type> epsilon_gt_cancer_crab_larvae( n_g_cancer_crab_larvae, n_t_jsoes );
epsilon_gt_cancer_crab_larvae = A_gs_cancer_crab_larvae * epsilon_st_cancer_crab_larvae;

// Probability of data conditional on random effects for cancer_crab_larvae
vector<Type> omega_i_cancer_crab_larvae( n_i_cancer_crab_larvae );
omega_i_cancer_crab_larvae = A_is_cancer_crab_larvae * omega_s_cancer_crab_larvae;
matrix<Type> epsilon_it_cancer_crab_larvae( n_i_cancer_crab_larvae, n_t_jsoes );
epsilon_it_cancer_crab_larvae = A_is_cancer_crab_larvae * epsilon_st_cancer_crab_larvae;

vector<Type> dhat_i_cancer_crab_larvae( n_i_cancer_crab_larvae );
dhat_i_cancer_crab_larvae.setZero();


for( int i=0; i<D_i_cancer_crab_larvae.size(); i++){
  dhat_i_cancer_crab_larvae(i) = exp(omega_i_cancer_crab_larvae(i) + epsilon_it_cancer_crab_larvae(i, t_i_cancer_crab_larvae(i)));
  jnll -= dtweedie( D_i_cancer_crab_larvae(i), dhat_i_cancer_crab_larvae(i), exp(ln_phi_cancer_crab_larvae), Type(1.0)+invlogit(finv_power_cancer_crab_larvae), true ) * weights_i_cancer_crab_larvae(i);
}



// SDM: modeled density in the projection grid for cancer_crab_larvae
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_cancer_crab_larvae( n_g_cancer_crab_larvae, n_t_jsoes );
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_cancer_crab_larvae; g++){
    ln_d_gt_cancer_crab_larvae(g,t) = omega_g_cancer_crab_larvae(g) + epsilon_gt_cancer_crab_larvae(g,t);
  }
}

// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_cancer_crab_larvae( n_g_cancer_crab_larvae, n_t_jsoes );

// loop through the array to convert to normal space
for (int i = 0; i < d_gt_cancer_crab_larvae.rows(); i++) {
  for (int j = 0; j < d_gt_cancer_crab_larvae.cols(); j++) {
    d_gt_cancer_crab_larvae(i,j) = exp(ln_d_gt_cancer_crab_larvae(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> cancer_crab_larvae_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  cancer_crab_larvae_index_of_abundance(t) = d_gt_cancer_crab_larvae.col(t).sum();
}

REPORT(cancer_crab_larvae_index_of_abundance);
// ADREPORT(cancer_crab_larvae_index_of_abundance);

// Reporting: cancer_crab_larvae
REPORT( dhat_i_cancer_crab_larvae );
REPORT( ln_tau_omega_cancer_crab_larvae );
REPORT( ln_tau_epsilon_cancer_crab_larvae );
REPORT( ln_kappa_cancer_crab_larvae );
REPORT( ln_phi_cancer_crab_larvae );
REPORT( finv_power_cancer_crab_larvae );
REPORT(SigmaO_cancer_crab_larvae);
REPORT(SigmaE_cancer_crab_larvae);
REPORT(Range_cancer_crab_larvae);
REPORT(ln_d_gt_cancer_crab_larvae);
REPORT( H_cancer_crab_larvae );
REPORT( logit_rhoE_cancer_crab_larvae );

// #### Non-Cancer Crab Larvae SDM ####

// Data for non_cancer_crab_larvae
DATA_VECTOR(D_i_non_cancer_crab_larvae);  // density of non_cancer_crab_larvae in measurement i
DATA_IVECTOR(t_i_non_cancer_crab_larvae); // index for the year of measurement i (non_cancer_crab_larvae)
DATA_VECTOR(weights_i_non_cancer_crab_larvae); // optional weights for non_cancer_crab_larvae - used for cAIC
DATA_VECTOR(dist_i_non_cancer_crab_larvae);  // distance from shore for measurement i


// Projection matrices for non_cancer_crab_larvae
DATA_SPARSE_MATRIX(A_is_non_cancer_crab_larvae); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_non_cancer_crab_larvae); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for non_cancer_crab_larvae
PARAMETER_VECTOR( beta_t_non_cancer_crab_larvae ); // fixed effects for year on density of non_cancer_crab_larvae
PARAMETER( beta_dist_non_cancer_crab_larvae ); // fixed effect for distance from shore on non_cancer_crab_larvae
PARAMETER( ln_tau_omega_non_cancer_crab_larvae ); // Tau parameter for spatial effects (omega) for non_cancer_crab_larvae
PARAMETER( ln_tau_epsilon_non_cancer_crab_larvae ); // Tau parameter for spatiotemporal effects (epsilon) for non_cancer_crab_larvae
PARAMETER( ln_kappa_non_cancer_crab_larvae ); // Kappa term for SPDE for non_cancer_crab_larvae
PARAMETER( ln_phi_non_cancer_crab_larvae ); // phi term in tweedie for non_cancer_crab_larvae
PARAMETER( finv_power_non_cancer_crab_larvae ); // power parameter in tweedie for non_cancer_crab_larvae

PARAMETER( logit_rhoE_non_cancer_crab_larvae ); // AR1 parameter for non_cancer_crab_larvae

// Random effects for non_cancer_crab_larvae
PARAMETER_VECTOR( omega_s_non_cancer_crab_larvae ); // vector of spatial random effects for non_cancer_crab_larvae
PARAMETER_MATRIX( epsilon_st_non_cancer_crab_larvae ); // vector of spatiotemporal random effects for non_cancer_crab_larvae

Type rhoE_non_cancer_crab_larvae = invlogit( logit_rhoE_non_cancer_crab_larvae ); // transform AR1 parameter

int n_i_non_cancer_crab_larvae = A_is_non_cancer_crab_larvae.rows(); // number of non_cancer_crab_larvae samples
int n_g_non_cancer_crab_larvae = A_gs_non_cancer_crab_larvae.rows(); // number of units in the projection grid for non_cancer_crab_larvae

// Anisotropy turned off in this model
// Anisotropy elements
matrix<Type> H_non_cancer_crab_larvae( 2, 2 );
// H(0,0) = exp(ln_H_input(0));
// H(1,0) = ln_H_input(1);
// H(0,1) = ln_H_input(1);
// H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
H_non_cancer_crab_larvae(0,0) = 1;
H_non_cancer_crab_larvae(0,1) = 0;
H_non_cancer_crab_larvae(1,0) = 0;
H_non_cancer_crab_larvae(1,1) = 1;

// implement anisotropy
Eigen::SparseMatrix<Type> Q_non_cancer_crab_larvae;
// Using INLA
DATA_STRUCT( spatial_list_non_cancer_crab_larvae, R_inla::spde_aniso_t );
// Build precision
Q_non_cancer_crab_larvae = R_inla::Q_spde( spatial_list_non_cancer_crab_larvae, exp(ln_kappa_non_cancer_crab_larvae), H_non_cancer_crab_larvae );

// Derived quantities
// SDM derived quantities for non_cancer_crab_larvae
Type Range_non_cancer_crab_larvae = sqrt(8) / exp( ln_kappa_non_cancer_crab_larvae );
Type SigmaO_non_cancer_crab_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_non_cancer_crab_larvae) * exp(2*ln_kappa_non_cancer_crab_larvae));
Type SigmaE_non_cancer_crab_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_non_cancer_crab_larvae) * exp(2*ln_kappa_non_cancer_crab_larvae));

// Probability of random effects for non_cancer_crab_larvae
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_non_cancer_crab_larvae), 1/exp(ln_tau_omega_non_cancer_crab_larvae) )( omega_s_non_cancer_crab_larvae );
// spatio-temporal random effect - scaled by ln_tau_epsilon
// spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  if( t==0 ){
    jnll += SCALE( GMRF(Q_non_cancer_crab_larvae), 1 / exp(ln_tau_epsilon_non_cancer_crab_larvae) / pow( 1.0-pow(rhoE_non_cancer_crab_larvae,2), 0.5 ) )( epsilon_st_non_cancer_crab_larvae.col(t) );
  }else{
    jnll += SCALE( GMRF(Q_non_cancer_crab_larvae), 1 / exp(ln_tau_epsilon_non_cancer_crab_larvae) )( epsilon_st_non_cancer_crab_larvae.col(t) - rhoE_non_cancer_crab_larvae*epsilon_st_non_cancer_crab_larvae.col(t-1) );
  }
}

// SDM projections for non_cancer_crab_larvae
vector<Type> omega_g_non_cancer_crab_larvae( n_g_non_cancer_crab_larvae );
omega_g_non_cancer_crab_larvae = A_gs_non_cancer_crab_larvae * omega_s_non_cancer_crab_larvae;
matrix<Type> epsilon_gt_non_cancer_crab_larvae( n_g_non_cancer_crab_larvae, n_t_jsoes );
epsilon_gt_non_cancer_crab_larvae = A_gs_non_cancer_crab_larvae * epsilon_st_non_cancer_crab_larvae;

// Probability of data conditional on random effects for non_cancer_crab_larvae
vector<Type> omega_i_non_cancer_crab_larvae( n_i_non_cancer_crab_larvae );
omega_i_non_cancer_crab_larvae = A_is_non_cancer_crab_larvae * omega_s_non_cancer_crab_larvae;
matrix<Type> epsilon_it_non_cancer_crab_larvae( n_i_non_cancer_crab_larvae, n_t_jsoes );
epsilon_it_non_cancer_crab_larvae = A_is_non_cancer_crab_larvae * epsilon_st_non_cancer_crab_larvae;

vector<Type> dhat_i_non_cancer_crab_larvae( n_i_non_cancer_crab_larvae );
dhat_i_non_cancer_crab_larvae.setZero();


for( int i=0; i<D_i_non_cancer_crab_larvae.size(); i++){
  dhat_i_non_cancer_crab_larvae(i) = exp( beta_t_non_cancer_crab_larvae(t_i_non_cancer_crab_larvae(i)) + beta_dist_non_cancer_crab_larvae*dist_i_non_cancer_crab_larvae(i) + 
    omega_i_non_cancer_crab_larvae(i) + epsilon_it_non_cancer_crab_larvae(i, t_i_non_cancer_crab_larvae(i)));
  jnll -= dtweedie( D_i_non_cancer_crab_larvae(i), dhat_i_non_cancer_crab_larvae(i), exp(ln_phi_non_cancer_crab_larvae), Type(1.0)+invlogit(finv_power_non_cancer_crab_larvae), true ) * weights_i_non_cancer_crab_larvae(i);
}



// SDM: modeled density in the projection grid for non_cancer_crab_larvae
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_non_cancer_crab_larvae( n_g_non_cancer_crab_larvae, n_t_jsoes );
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_non_cancer_crab_larvae; g++){
    ln_d_gt_non_cancer_crab_larvae(g,t) = beta_t_non_cancer_crab_larvae(t) + beta_dist_non_cancer_crab_larvae*dist_g(g) + omega_g_non_cancer_crab_larvae(g) + epsilon_gt_non_cancer_crab_larvae(g,t);
  }
}

// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_non_cancer_crab_larvae( n_g_non_cancer_crab_larvae, n_t_jsoes );

// loop through the array to convert to normal space
for (int i = 0; i < d_gt_non_cancer_crab_larvae.rows(); i++) {
  for (int j = 0; j < d_gt_non_cancer_crab_larvae.cols(); j++) {
    d_gt_non_cancer_crab_larvae(i,j) = exp(ln_d_gt_non_cancer_crab_larvae(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> non_cancer_crab_larvae_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  non_cancer_crab_larvae_index_of_abundance(t) = d_gt_non_cancer_crab_larvae.col(t).sum();
}

REPORT(non_cancer_crab_larvae_index_of_abundance);
// ADREPORT(non_cancer_crab_larvae_index_of_abundance);

// Reporting: non_cancer_crab_larvae
REPORT( dhat_i_non_cancer_crab_larvae );
REPORT( beta_t_non_cancer_crab_larvae );
REPORT( beta_dist_non_cancer_crab_larvae );
REPORT( ln_tau_omega_non_cancer_crab_larvae );
REPORT( ln_tau_epsilon_non_cancer_crab_larvae );
REPORT( ln_kappa_non_cancer_crab_larvae );
REPORT( ln_phi_non_cancer_crab_larvae );
REPORT( finv_power_non_cancer_crab_larvae );
REPORT(SigmaO_non_cancer_crab_larvae);
REPORT(SigmaE_non_cancer_crab_larvae);
REPORT(Range_non_cancer_crab_larvae);
REPORT(ln_d_gt_non_cancer_crab_larvae);
REPORT( H_non_cancer_crab_larvae );

// #### Shrimp Larvae SDM ####

// Data for shrimp_larvae
DATA_VECTOR(D_i_shrimp_larvae);  // density of shrimp_larvae in measurement i
DATA_IVECTOR(t_i_shrimp_larvae); // index for the year of measurement i (shrimp_larvae)
DATA_VECTOR(weights_i_shrimp_larvae); // optional weights for shrimp_larvae - used for cAIC
DATA_VECTOR(dist_i_shrimp_larvae); // distance from shore for shrimp_larvae


// Projection matrices for shrimp_larvae
DATA_SPARSE_MATRIX(A_is_shrimp_larvae); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_shrimp_larvae); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain

// Data for shrimp larvae smoothers
DATA_IVECTOR(b_smooth_start_shrimp_larvae);
DATA_STRUCT(Zs_shrimp_larvae, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
DATA_STRUCT(proj_Zs_shrimp_larvae, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
DATA_MATRIX(Xs_shrimp_larvae); // smoother linear effect matrix
DATA_MATRIX(proj_Xs_shrimp_larvae); // smoother linear effect matrix

// Parameters for shrimp_larvae
PARAMETER( ln_tau_omega_shrimp_larvae ); // Tau parameter for spatial effects (omega) for shrimp_larvae
PARAMETER( ln_tau_epsilon_shrimp_larvae ); // Tau parameter for spatiotemporal effects (epsilon) for shrimp_larvae
PARAMETER( ln_kappa_shrimp_larvae ); // Kappa term for SPDE for shrimp_larvae
PARAMETER( ln_phi_shrimp_larvae ); // phi term in tweedie for shrimp_larvae
PARAMETER( finv_power_shrimp_larvae ); // power parameter in tweedie for shrimp_larvae

PARAMETER( logit_rhoE_shrimp_larvae ); // AR1 parameter for shrimp_larvae

// smoother parameters for shrimp_larvae
PARAMETER_VECTOR(bs_shrimp_larvae); // smoother linear effects
PARAMETER_VECTOR(b_smooth_shrimp_larvae);  // P-spline smooth parameters
PARAMETER_VECTOR(ln_smooth_sigma_shrimp_larvae);  // variances of spline REs if included

// Random effects for shrimp_larvae
PARAMETER_VECTOR( omega_s_shrimp_larvae ); // vector of spatial random effects for shrimp_larvae
PARAMETER_MATRIX( epsilon_st_shrimp_larvae ); // vector of spatiotemporal random effects for shrimp_larvae

Type rhoE_shrimp_larvae = invlogit( logit_rhoE_shrimp_larvae ); // transform AR1 parameter

int n_i_shrimp_larvae = A_is_shrimp_larvae.rows(); // number of shrimp_larvae samples
int n_g_shrimp_larvae = A_gs_shrimp_larvae.rows(); // number of units in the projection grid for shrimp_larvae

// p-splines/smoothers
vector<Type> eta_smooth_i_shrimp_larvae(n_i_shrimp_larvae);
eta_smooth_i_shrimp_larvae.setZero();
for (int s = 0; s < b_smooth_start_shrimp_larvae.size(); s++) { // iterate over # of smooth elements
  vector<Type> beta_s_shrimp_larvae(Zs_shrimp_larvae(s).cols());
  beta_s_shrimp_larvae.setZero();
  for (int j = 0; j < beta_s_shrimp_larvae.size(); j++) {
    beta_s_shrimp_larvae(j) = b_smooth_shrimp_larvae(b_smooth_start_shrimp_larvae(s) + j);
    // PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
    jnll -= dnorm(beta_s_shrimp_larvae(j), Type(0), exp(ln_smooth_sigma_shrimp_larvae(s)), true);
  }
  eta_smooth_i_shrimp_larvae += Zs_shrimp_larvae(s) * beta_s_shrimp_larvae;
}
// eta_smooth_i += Xs * vector<Type>(bs.col(1));
eta_smooth_i_shrimp_larvae += Xs_shrimp_larvae * bs_shrimp_larvae;
REPORT(b_smooth_shrimp_larvae);     // smooth coefficients for penalized splines
REPORT(ln_smooth_sigma_shrimp_larvae); // standard deviations of smooth random effects, in log-space

// Anisotropy elements for shrimp_larvae
// Anisotropy turned off in this model
matrix<Type> H_shrimp_larvae( 2, 2 );
// H(0,0) = exp(ln_H_input(0));
// H(1,0) = ln_H_input(1);
// H(0,1) = ln_H_input(1);
// H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
H_shrimp_larvae(0,0) = 1;
H_shrimp_larvae(0,1) = 0;
H_shrimp_larvae(1,0) = 0;
H_shrimp_larvae(1,1) = 1;

// implement anisotropy
Eigen::SparseMatrix<Type> Q_shrimp_larvae;
// Using INLA
DATA_STRUCT( spatial_list_shrimp_larvae, R_inla::spde_aniso_t );
// Build precision
Q_shrimp_larvae = R_inla::Q_spde( spatial_list_shrimp_larvae, exp(ln_kappa_shrimp_larvae), H_shrimp_larvae );

// Derived quantities
// SDM derived quantities for shrimp_larvae
Type Range_shrimp_larvae = sqrt(8) / exp( ln_kappa_shrimp_larvae );
Type SigmaO_shrimp_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_shrimp_larvae) * exp(2*ln_kappa_shrimp_larvae));
Type SigmaE_shrimp_larvae = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_shrimp_larvae) * exp(2*ln_kappa_shrimp_larvae));

// Probability of random effects for shrimp_larvae
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_shrimp_larvae), 1/exp(ln_tau_omega_shrimp_larvae) )( omega_s_shrimp_larvae );
// spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  if( t==0 ){
    jnll += SCALE( GMRF(Q_shrimp_larvae), 1 / exp(ln_tau_epsilon_shrimp_larvae) / pow( 1.0-pow(rhoE_shrimp_larvae,2), 0.5 ) )( epsilon_st_shrimp_larvae.col(t) );
  }else{
    jnll += SCALE( GMRF(Q_shrimp_larvae), 1 / exp(ln_tau_epsilon_shrimp_larvae) )( epsilon_st_shrimp_larvae.col(t) - rhoE_shrimp_larvae*epsilon_st_shrimp_larvae.col(t-1) );
  }
}

// SDM projections for shrimp_larvae
vector<Type> omega_g_shrimp_larvae( n_g_shrimp_larvae );
omega_g_shrimp_larvae = A_gs_shrimp_larvae * omega_s_shrimp_larvae;
matrix<Type> epsilon_gt_shrimp_larvae( n_g_shrimp_larvae, n_t_jsoes );
epsilon_gt_shrimp_larvae = A_gs_shrimp_larvae * epsilon_st_shrimp_larvae;

// Probability of data conditional on random effects for shrimp_larvae
vector<Type> omega_i_shrimp_larvae( n_i_shrimp_larvae );
omega_i_shrimp_larvae = A_is_shrimp_larvae * omega_s_shrimp_larvae;
matrix<Type> epsilon_it_shrimp_larvae( n_i_shrimp_larvae, n_t_jsoes );
epsilon_it_shrimp_larvae = A_is_shrimp_larvae * epsilon_st_shrimp_larvae;

vector<Type> dhat_i_shrimp_larvae( n_i_shrimp_larvae );
dhat_i_shrimp_larvae.setZero();


for( int i=0; i<D_i_shrimp_larvae.size(); i++){
  dhat_i_shrimp_larvae(i) = exp( omega_i_shrimp_larvae(i) + epsilon_it_shrimp_larvae(i, t_i_shrimp_larvae(i)) + eta_smooth_i_shrimp_larvae(i));
  jnll -= dtweedie( D_i_shrimp_larvae(i), dhat_i_shrimp_larvae(i), exp(ln_phi_shrimp_larvae), Type(1.0)+invlogit(finv_power_shrimp_larvae), true ) * weights_i_shrimp_larvae(i);
}



// SDM: modeled density in the projection grid for shrimp_larvae
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_shrimp_larvae( n_g_shrimp_larvae, n_t_jsoes );

// estimate the smooth effects in the projection grid
// smoothers for projection grid
vector<Type> proj_smooth_i_shrimp_larvae(n_g_shrimp_larvae*n_t_jsoes);
proj_smooth_i_shrimp_larvae.setZero();
for (int s = 0; s < b_smooth_start_shrimp_larvae.size(); s++) { // iterate over # of smooth elements
  vector<Type> beta_s_shrimp_larvae(proj_Zs_shrimp_larvae(s).cols());
  beta_s_shrimp_larvae.setZero();
  for (int j = 0; j < beta_s_shrimp_larvae.size(); j++) {
    beta_s_shrimp_larvae(j) = b_smooth_shrimp_larvae(b_smooth_start_shrimp_larvae(s) + j);
  }
  proj_smooth_i_shrimp_larvae += proj_Zs_shrimp_larvae(s) * beta_s_shrimp_larvae;
}
proj_smooth_i_shrimp_larvae += proj_Xs_shrimp_larvae * bs_shrimp_larvae;

// add smoothed effects, fixed effects, and random effects
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_shrimp_larvae; g++){
    ln_d_gt_shrimp_larvae(g,t) = omega_g_shrimp_larvae(g) + epsilon_gt_shrimp_larvae(g,t) + proj_smooth_i_shrimp_larvae(g + n_g_shrimp_larvae * t);
  }
}

// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_shrimp_larvae( n_g_shrimp_larvae, n_t_jsoes );

// loop through the array to convert to normal space
for (int i = 0; i < d_gt_shrimp_larvae.rows(); i++) {
  for (int j = 0; j < d_gt_shrimp_larvae.cols(); j++) {
    d_gt_shrimp_larvae(i,j) = exp(ln_d_gt_shrimp_larvae(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> shrimp_larvae_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  shrimp_larvae_index_of_abundance(t) = d_gt_shrimp_larvae.col(t).sum();
}

REPORT(shrimp_larvae_index_of_abundance);
// ADREPORT(shrimp_larvae_index_of_abundance);

// Reporting: shrimp_larvae
REPORT( dhat_i_shrimp_larvae );
REPORT( ln_tau_omega_shrimp_larvae );
REPORT( ln_tau_epsilon_shrimp_larvae );
REPORT( ln_kappa_shrimp_larvae );
REPORT( ln_phi_shrimp_larvae );
REPORT( finv_power_shrimp_larvae );
REPORT(SigmaO_shrimp_larvae);
REPORT(SigmaE_shrimp_larvae);
REPORT(Range_shrimp_larvae);
REPORT(ln_d_gt_shrimp_larvae);
REPORT( logit_rhoE_shrimp_larvae );

// #### Hyperiid Amphipods SDM ####

// Data for hyperiid_amphipods
DATA_VECTOR(D_i_hyperiid_amphipods);  // density of hyperiid_amphipods in measurement i
DATA_IVECTOR(t_i_hyperiid_amphipods); // index for the year of measurement i (hyperiid_amphipods)
DATA_VECTOR(weights_i_hyperiid_amphipods); // optional weights for hyperiid_amphipods - used for cAIC


// Projection matrices for hyperiid_amphipods
DATA_SPARSE_MATRIX(A_is_hyperiid_amphipods); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_hyperiid_amphipods); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for hyperiid_amphipods
PARAMETER( ln_tau_omega_hyperiid_amphipods ); // Tau parameter for spatial effects (omega) for hyperiid_amphipods
PARAMETER( ln_tau_epsilon_hyperiid_amphipods ); // Tau parameter for spatiotemporal effects (epsilon) for hyperiid_amphipods
PARAMETER( ln_kappa_hyperiid_amphipods ); // Kappa term for SPDE for hyperiid_amphipods
PARAMETER( ln_phi_hyperiid_amphipods ); // phi term in tweedie for hyperiid_amphipods
PARAMETER( finv_power_hyperiid_amphipods ); // power parameter in tweedie for hyperiid_amphipods


PARAMETER_VECTOR( ln_H_input_hyperiid_amphipods ); // anisotropy input for hyperiid_amphipods

// Random effects for hyperiid_amphipods
PARAMETER_VECTOR( omega_s_hyperiid_amphipods ); // vector of spatial random effects for hyperiid_amphipods
PARAMETER_MATRIX( epsilon_st_hyperiid_amphipods ); // vector of spatiotemporal random effects for hyperiid_amphipods


int n_i_hyperiid_amphipods = A_is_hyperiid_amphipods.rows(); // number of hyperiid_amphipods samples
int n_g_hyperiid_amphipods = A_gs_hyperiid_amphipods.rows(); // number of units in the projection grid for hyperiid_amphipods

// priors for anisotropy parameters for hyperiid_amphipods
Type ln_H_0_hyperiid_amphipods_mean = 0.0;
Type ln_H_0_hyperiid_amphipods_sd = 0.5;
Type ln_H_1_hyperiid_amphipods_mean = 0.0;
Type ln_H_1_hyperiid_amphipods_sd = 0.5;

jnll -= dnorm(ln_H_input_hyperiid_amphipods(0), ln_H_0_hyperiid_amphipods_mean, ln_H_0_hyperiid_amphipods_sd, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_hyperiid_amphipods(1), ln_H_1_hyperiid_amphipods_mean, ln_H_1_hyperiid_amphipods_sd, true); // Anisotropic correlation

// Anisotropy elements for hyperiid_amphipods
matrix<Type> H_hyperiid_amphipods( 2, 2 );
H_hyperiid_amphipods(0,0) = exp(ln_H_input_hyperiid_amphipods(0));
H_hyperiid_amphipods(1,0) = ln_H_input_hyperiid_amphipods(1);
H_hyperiid_amphipods(0,1) = ln_H_input_hyperiid_amphipods(1);
H_hyperiid_amphipods(1,1) = (1+ln_H_input_hyperiid_amphipods(1)*ln_H_input_hyperiid_amphipods(1)) / exp(ln_H_input_hyperiid_amphipods(0));

// implement anisotropy for hyperiid_amphipods
Eigen::SparseMatrix<Type> Q_hyperiid_amphipods;
// Using INLA
DATA_STRUCT( spatial_list_hyperiid_amphipods, R_inla::spde_aniso_t );
// Build precision matrix for hyperiid_amphipods
Q_hyperiid_amphipods = R_inla::Q_spde( spatial_list_hyperiid_amphipods, exp(ln_kappa_hyperiid_amphipods), H_hyperiid_amphipods );

// Derived quantities
// SDM derived quantities for hyperiid_amphipods
Type Range_hyperiid_amphipods = sqrt(8) / exp( ln_kappa_hyperiid_amphipods );
Type SigmaO_hyperiid_amphipods = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_hyperiid_amphipods) * exp(2*ln_kappa_hyperiid_amphipods));
Type SigmaE_hyperiid_amphipods = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_hyperiid_amphipods) * exp(2*ln_kappa_hyperiid_amphipods));

// Probability of random effects for hyperiid_amphipods
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_hyperiid_amphipods), 1/exp(ln_tau_omega_hyperiid_amphipods) )( omega_s_hyperiid_amphipods );
// spatio-temporal random effect - scaled by ln_tau_epsilon
for( int t=0; t<n_t_jsoes; t++){
  jnll += SCALE( GMRF(Q_hyperiid_amphipods), 1/exp(ln_tau_epsilon_hyperiid_amphipods) )( epsilon_st_hyperiid_amphipods.col(t) );
}

// SDM projections for hyperiid_amphipods
vector<Type> omega_g_hyperiid_amphipods( n_g_hyperiid_amphipods );
omega_g_hyperiid_amphipods = A_gs_hyperiid_amphipods * omega_s_hyperiid_amphipods;
matrix<Type> epsilon_gt_hyperiid_amphipods( n_g_hyperiid_amphipods, n_t_jsoes );
epsilon_gt_hyperiid_amphipods = A_gs_hyperiid_amphipods * epsilon_st_hyperiid_amphipods;

// Probability of data conditional on random effects for hyperiid_amphipods
vector<Type> omega_i_hyperiid_amphipods( n_i_hyperiid_amphipods );
omega_i_hyperiid_amphipods = A_is_hyperiid_amphipods * omega_s_hyperiid_amphipods;
matrix<Type> epsilon_it_hyperiid_amphipods( n_i_hyperiid_amphipods, n_t_jsoes );
epsilon_it_hyperiid_amphipods = A_is_hyperiid_amphipods * epsilon_st_hyperiid_amphipods;

vector<Type> dhat_i_hyperiid_amphipods( n_i_hyperiid_amphipods );
dhat_i_hyperiid_amphipods.setZero();


for( int i=0; i<D_i_hyperiid_amphipods.size(); i++){
  dhat_i_hyperiid_amphipods(i) = exp(omega_i_hyperiid_amphipods(i) + epsilon_it_hyperiid_amphipods(i, t_i_hyperiid_amphipods(i)));
  jnll -= dtweedie( D_i_hyperiid_amphipods(i), dhat_i_hyperiid_amphipods(i), exp(ln_phi_hyperiid_amphipods), Type(1.0)+invlogit(finv_power_hyperiid_amphipods), true ) * weights_i_hyperiid_amphipods(i);
}



// SDM: modeled density in the projection grid for hyperiid_amphipods
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_hyperiid_amphipods( n_g_hyperiid_amphipods, n_t_jsoes );
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_hyperiid_amphipods; g++){
    ln_d_gt_hyperiid_amphipods(g,t) = omega_g_hyperiid_amphipods(g) + epsilon_gt_hyperiid_amphipods(g,t);
  }
}

// SDM: Estimate an index of abundance
// sum across all cells, multiply by the area of each cell
// convert to normal space
array<Type> d_gt_hyperiid_amphipods( n_g_hyperiid_amphipods, n_t_jsoes );
// 
// loop through the array to convert to normal space
for (int i = 0; i < d_gt_hyperiid_amphipods.rows(); i++) {
  for (int j = 0; j < d_gt_hyperiid_amphipods.cols(); j++) {
    d_gt_hyperiid_amphipods(i,j) = exp(ln_d_gt_hyperiid_amphipods(i,j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> hyperiid_amphipods_index_of_abundance(n_t_jsoes);
for( int t=0; t<n_t_jsoes; t++){
  hyperiid_amphipods_index_of_abundance(t) = d_gt_hyperiid_amphipods.col(t).sum();
}

REPORT(hyperiid_amphipods_index_of_abundance);
// ADREPORT(hyperiid_amphipods_index_of_abundance);

// Reporting: hyperiid_amphipods
REPORT( dhat_i_hyperiid_amphipods );
REPORT( ln_tau_omega_hyperiid_amphipods );
REPORT( ln_tau_epsilon_hyperiid_amphipods );
REPORT( ln_kappa_hyperiid_amphipods );
REPORT( ln_phi_hyperiid_amphipods );
REPORT( finv_power_hyperiid_amphipods );
REPORT(SigmaO_hyperiid_amphipods);
REPORT(SigmaE_hyperiid_amphipods);
REPORT(Range_hyperiid_amphipods);
REPORT(ln_d_gt_hyperiid_amphipods);
REPORT( H_hyperiid_amphipods );

// #### summarize into an aggregate prey field
// loop through the projection grid and populate with the sum of different prey items

// create an empty array to store the prey field
// any n_g object is fine because we're projecting to the same area (the JSOES domain) for all taxa
array<Type> d_gt_prey_field( n_g_hyperiid_amphipods, n_t_shared );

// loop through the array to summarize the prey field
for (int i = 0; i < d_gt_prey_field.rows(); i++) {
  for (int j = 0; j < d_gt_prey_field.cols(); j++) {
    d_gt_prey_field(i,j) = d_gt_cancer_crab_larvae(i,t_jsoes_indices(j)) + d_gt_non_cancer_crab_larvae(i,t_jsoes_indices(j)) +
      d_gt_shrimp_larvae(i,t_jsoes_indices(j)) + d_gt_hyperiid_amphipods(i,t_jsoes_indices(j));
  }
}

// create a vector to store the index of abundance and calculate it for each year
vector<Type> prey_field_index_of_abundance(n_t_shared);
for( int t=0; t<n_t_shared; t++){
  prey_field_index_of_abundance(t) = d_gt_prey_field.col(t).sum();
}

REPORT(prey_field_index_of_abundance);
ADREPORT(prey_field_index_of_abundance);


// PREDATORS


// #### CALCULATE OVERLAP WITH INDIVIDUAL PREY TAXA ####
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

// Define denominator for each of the prey items - only need to do this once, and it gets re-used
array<Type> d_gt_prey_field_sq( n_g_csyif, n_t_shared ); // array to store the squared proportional densities of prey_field in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_shared; t++){
  for( int g=0; g<n_g_cssif; g++){
    d_gt_prey_field_sq(g,t) = pow(d_gt_prey_field(g,t)/d_gt_prey_field.col(t).sum(),2); // calculate value for each cell
  }
}
array<Type> d_gt_cancer_crab_larvae_sq( n_g_csyif, n_t_jsoes ); // array to store the squared proportional densities of cancer_crab_larvae in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_cancer_crab_larvae_sq(g,t) = pow(d_gt_cancer_crab_larvae(g,t)/d_gt_cancer_crab_larvae.col(t).sum(),2); // calculate value for each cell
  }
}
array<Type> d_gt_hyperiid_amphipods_sq( n_g_csyif, n_t_jsoes ); // array to store the squared proportional densities of hyperiid_amphipods in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_hyperiid_amphipods_sq(g,t) = pow(d_gt_hyperiid_amphipods(g,t)/d_gt_hyperiid_amphipods.col(t).sum(),2); // calculate value for each cell
  }
}
array<Type> d_gt_non_cancer_crab_larvae_sq( n_g_csyif, n_t_jsoes ); // array to store the squared proportional densities of non_cancer_crab_larvae in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_non_cancer_crab_larvae_sq(g,t) = pow(d_gt_non_cancer_crab_larvae(g,t)/d_gt_non_cancer_crab_larvae.col(t).sum(),2); // calculate value for each cell
  }
}
array<Type> d_gt_shrimp_larvae_sq( n_g_csyif, n_t_jsoes ); // array to store the squared proportional densities of shrimp_larvae in each cell (for denominators)
// calculate the squared densities for the denominator
for( int t=0; t<n_t_jsoes; t++){
  for( int g=0; g<n_g_csyif; g++){
    d_gt_shrimp_larvae_sq(g,t) = pow(d_gt_shrimp_larvae(g,t)/d_gt_shrimp_larvae.col(t).sum(),2); // calculate value for each cell
  }
}


// // #### csyif and cancer_crab_larvae ####
// vector<Type> pianka_o_csyif_cancer_crab_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_csyif_cancer_crab_larvae( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// 
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_csyif; g++){
//     ov_gt_num_csyif_cancer_crab_larvae(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_cancer_crab_larvae(g,t_jsoes_indices(t))/d_gt_cancer_crab_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_csyif_cancer_crab_larvae_t(t) = ov_gt_num_csyif_cancer_crab_larvae.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_cancer_crab_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for csyif and cancer_crab_larvae
// REPORT(pianka_o_csyif_cancer_crab_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_csyif_cancer_crab_larvae_t); // report overlap for full JSOES domain
// 
// // #### csyif and hyperiid_amphipods ####
// vector<Type> pianka_o_csyif_hyperiid_amphipods_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_csyif_hyperiid_amphipods( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_csyif; g++){
//     ov_gt_num_csyif_hyperiid_amphipods(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_hyperiid_amphipods(g,t_jsoes_indices(t))/d_gt_hyperiid_amphipods.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_csyif_hyperiid_amphipods_t(t) = ov_gt_num_csyif_hyperiid_amphipods.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_hyperiid_amphipods_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for csyif and hyperiid_amphipods
// REPORT(pianka_o_csyif_hyperiid_amphipods_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_csyif_hyperiid_amphipods_t); // report overlap for full JSOES domain
// 
// 
// // #### csyif and non_cancer_crab_larvae ####
// vector<Type> pianka_o_csyif_non_cancer_crab_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_csyif_non_cancer_crab_larvae( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_csyif; g++){
//     ov_gt_num_csyif_non_cancer_crab_larvae(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_non_cancer_crab_larvae(g,t_jsoes_indices(t))/d_gt_non_cancer_crab_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_csyif_non_cancer_crab_larvae_t(t) = ov_gt_num_csyif_non_cancer_crab_larvae.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_non_cancer_crab_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for csyif and non_cancer_crab_larvae
// REPORT(pianka_o_csyif_non_cancer_crab_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_csyif_non_cancer_crab_larvae_t); // report overlap for full JSOES domain
// 
// // #### csyif and shrimp_larvae ####
// vector<Type> pianka_o_csyif_shrimp_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_csyif_shrimp_larvae( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_csyif; g++){
//     ov_gt_num_csyif_shrimp_larvae(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_shrimp_larvae(g,t_jsoes_indices(t))/d_gt_shrimp_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_csyif_shrimp_larvae_t(t) = ov_gt_num_csyif_shrimp_larvae.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_shrimp_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for csyif and shrimp_larvae
// REPORT(pianka_o_csyif_shrimp_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_csyif_shrimp_larvae_t); // report overlap for full JSOES domain
// 
// #### csyif and aggregate prey field ####
vector<Type> pianka_o_csyif_prey_field_t( n_t_shared ); // vector to store the value per year
array<Type> ov_gt_num_csyif_prey_field( n_g_csyif, n_t_shared ); // array to store the numerator values in each cell
// calculate the numerator values
for( int t=0; t<n_t_shared; t++){
  for( int g=0; g<n_g_csyif; g++){
    ov_gt_num_csyif_prey_field(g,t) = d_gt_csyif(g,t_jsoes_indices(t))/d_gt_csyif.col(t_jsoes_indices(t)).sum() * d_gt_prey_field(g,t)/d_gt_prey_field.col(t).sum(); // calculate value for each cell
  }
}


// calculate the index for each year
for( int t=0; t<n_t_shared; t++){
  pianka_o_csyif_prey_field_t(t) = ov_gt_num_csyif_prey_field.col(t).sum()/sqrt(d_gt_csyif_sq.col(t_jsoes_indices(t)).sum()*d_gt_prey_field_sq.col(t).sum());
}

// Reporting: overlap for csyif and prey_field
REPORT(pianka_o_csyif_prey_field_t); // report overlap for full JSOES domain
// ADREPORT(pianka_o_csyif_prey_field_t); // report overlap for full JSOES domain
// 
// // #### cssif and cancer_crab_larvae ####
// vector<Type> pianka_o_cssif_cancer_crab_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_cssif_cancer_crab_larvae( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// 
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_cssif; g++){
//     ov_gt_num_cssif_cancer_crab_larvae(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_cancer_crab_larvae(g,t_jsoes_indices(t))/d_gt_cancer_crab_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_cssif_cancer_crab_larvae_t(t) = ov_gt_num_cssif_cancer_crab_larvae.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_cancer_crab_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for cssif and cancer_crab_larvae
// REPORT(pianka_o_cssif_cancer_crab_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_cssif_cancer_crab_larvae_t); // report overlap for full JSOES domain
// 
// // #### cssif and hyperiid_amphipods ####
// vector<Type> pianka_o_cssif_hyperiid_amphipods_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_cssif_hyperiid_amphipods( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_cssif; g++){
//     ov_gt_num_cssif_hyperiid_amphipods(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_hyperiid_amphipods(g,t_jsoes_indices(t))/d_gt_hyperiid_amphipods.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_cssif_hyperiid_amphipods_t(t) = ov_gt_num_cssif_hyperiid_amphipods.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_hyperiid_amphipods_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for cssif and hyperiid_amphipods
// REPORT(pianka_o_cssif_hyperiid_amphipods_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_cssif_hyperiid_amphipods_t); // report overlap for full JSOES domain
// 
// 
// // #### cssif and non_cancer_crab_larvae ####
// vector<Type> pianka_o_cssif_non_cancer_crab_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_cssif_non_cancer_crab_larvae( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_cssif; g++){
//     ov_gt_num_cssif_non_cancer_crab_larvae(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_non_cancer_crab_larvae(g,t_jsoes_indices(t))/d_gt_non_cancer_crab_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_cssif_non_cancer_crab_larvae_t(t) = ov_gt_num_cssif_non_cancer_crab_larvae.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_non_cancer_crab_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for cssif and non_cancer_crab_larvae
// REPORT(pianka_o_cssif_non_cancer_crab_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_cssif_non_cancer_crab_larvae_t); // report overlap for full JSOES domain
// 
// // #### cssif and shrimp_larvae ####
// vector<Type> pianka_o_cssif_shrimp_larvae_t( n_t_shared ); // vector to store the value per year
// array<Type> ov_gt_num_cssif_shrimp_larvae( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// // calculate the numerator values
// for( int t=0; t<n_t_shared; t++){
//   for( int g=0; g<n_g_cssif; g++){
//     ov_gt_num_cssif_shrimp_larvae(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_shrimp_larvae(g,t_jsoes_indices(t))/d_gt_shrimp_larvae.col(t_jsoes_indices(t)).sum(); // calculate value for each cell
//   }
// }
// 
// 
// // calculate the index for each year
// for( int t=0; t<n_t_shared; t++){
//   pianka_o_cssif_shrimp_larvae_t(t) = ov_gt_num_cssif_shrimp_larvae.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_shrimp_larvae_sq.col(t_jsoes_indices(t)).sum());
// }
// 
// // Reporting: overlap for cssif and shrimp_larvae
// REPORT(pianka_o_cssif_shrimp_larvae_t); // report overlap for full JSOES domain
// // ADREPORT(pianka_o_cssif_shrimp_larvae_t); // report overlap for full JSOES domain
// 
// #### cssif and aggregate prey field ####
vector<Type> pianka_o_cssif_prey_field_t( n_t_shared ); // vector to store the value per year
array<Type> ov_gt_num_cssif_prey_field( n_g_cssif, n_t_shared ); // array to store the numerator values in each cell
// calculate the numerator values
for( int t=0; t<n_t_shared; t++){
  for( int g=0; g<n_g_cssif; g++){
    ov_gt_num_cssif_prey_field(g,t) = d_gt_cssif(g,t_jsoes_indices(t))/d_gt_cssif.col(t_jsoes_indices(t)).sum() * d_gt_prey_field(g,t)/d_gt_prey_field.col(t).sum(); // calculate value for each cell
  }
}


// calculate the index for each year
for( int t=0; t<n_t_shared; t++){
  pianka_o_cssif_prey_field_t(t) = ov_gt_num_cssif_prey_field.col(t).sum()/sqrt(d_gt_cssif_sq.col(t_jsoes_indices(t)).sum()*d_gt_prey_field_sq.col(t).sum());
}

// Reporting: overlap for cssif and prey_field
REPORT(pianka_o_cssif_prey_field_t); // report overlap for full JSOES domain
// ADREPORT(pianka_o_cssif_prey_field_t); // report overlap for full JSOES domain
return jnll;
}