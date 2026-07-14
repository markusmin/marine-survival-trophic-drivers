// Best fit RE model, with density ~ factor(year) + s(dist)


#include <TMB.hpp>
#include <algorithm>
#include "utils.h"

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;
//using namespace R_inla;  // Not loaded globally, but used below for anisotropy

  // Data
  DATA_VECTOR(D_i);  // density of measurement i
  DATA_IVECTOR(t_i); // index for the year of measurement i
  DATA_IVECTOR(year_factor_t_i); // index for the as.factor(year) parameter
  DATA_INTEGER(n_t); // number of years in the dataset
  DATA_VECTOR(weights_i); // optional weights - used for cAIC
  // DATA_VECTOR(temp_i);  // temperature for measurement i
  DATA_VECTOR(dist_i);  // distance from shore for measurement i
  
  // Data for smoothers
  DATA_IVECTOR(b_smooth_start);
  DATA_STRUCT(Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(proj_Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_MATRIX(proj_Xs); // smoother linear effect matrix
  
  // Projection matrices for species A
  // A_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
  DATA_SPARSE_MATRIX(A_is);
  // A_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
  DATA_SPARSE_MATRIX(A_gs);

  // DATA_MATRIX(temp_gt); // temperature at each location in each year
  DATA_VECTOR(dist_g); // distance from shore at each location
  
  
  // Parameters
  
  // Parameters for SDM for species A
  PARAMETER_VECTOR( beta_t );
  // PARAMETER( beta_temp );
  // PARAMETER( beta_dist );
  // PARAMETER( ln_tau_omega );
  PARAMETER( ln_tau_epsilon );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi ); // phi term in tweedie
  PARAMETER( finv_power ); // power parameter in tweedie
  
  // AR1 parameters
  PARAMETER( logit_rhoE );
  
  // PARAMETER_VECTOR( ln_H_input ); // anisotropy input.
  
  // smoother parameters
  PARAMETER_VECTOR(bs); // smoother linear effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  
  // Random effects
  
  // PARAMETER_VECTOR( omega_s );
  PARAMETER_MATRIX( epsilon_st );
  
  // Objective function
  Type jnll = 0;
  int n_i = A_is.rows();
  int n_g = A_gs.rows();
  
  // p-splines/smoothers
  vector<Type> eta_smooth_i(n_i);
  eta_smooth_i.setZero();
  for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
    vector<Type> beta_s(Zs(s).cols());
    beta_s.setZero();
    for (int j = 0; j < beta_s.size(); j++) {
      beta_s(j) = b_smooth(b_smooth_start(s) + j);
      // PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
      jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
    }
    eta_smooth_i += Zs(s) * beta_s;
  }
  // eta_smooth_i += Xs * vector<Type>(bs.col(1));
  eta_smooth_i += Xs * bs;
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  
  // AR1 parameters
  Type rhoE = invlogit( logit_rhoE );
  
  // priors for anisotropy parameters
  // turned OFF in this model
  // Type ln_H_0_mean = 0.0;
  // Type ln_H_0_sd = 0.5;
  // Type ln_H_1_mean = 0.0;
  // Type ln_H_1_sd = 0.5;
  
  // jnll -= dnorm(ln_H_input(0), ln_H_0_mean, ln_H_0_sd, true); // Northings anisotropy
  // jnll -= dnorm(ln_H_input(1), ln_H_1_mean, ln_H_1_sd, true); // Anisotropic correlation
  
  // Anisotropy elements
  matrix<Type> H( 2, 2 );
  // H(0,0) = exp(ln_H_input(0));
  // H(1,0) = ln_H_input(1);
  // H(0,1) = ln_H_input(1);
  // H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  H(0,0) = 1;
  H(0,1) = 0;
  H(1,0) = 0;
  H(1,1) = 1;
  
  // implement anisotropy
  Eigen::SparseMatrix<Type> Q;
  // Using INLA
  DATA_STRUCT( spatial_list, R_inla::spde_aniso_t );
  // Build precision
  Q = R_inla::Q_spde( spatial_list, exp(ln_kappa), H );

  // Derived quantities
  // SDM derived quantities for species A
  Type Range = sqrt(8) / exp( ln_kappa );
  // Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon) * exp(2*ln_kappa));
  
  
  // Probability of random effects

  // SDM:
  // Eigen::SparseMatrix<Type> Q_a = exp(4*ln_kappa_a)*M0_a + Type(2.0)*exp(2*ln_kappa_a)*M1_a + M2_a;
  // spatial random effect - scaled by ln_tau_omega
  // jnll += SCALE( GMRF(Q), 1/exp(ln_tau_omega) )( omega_s );
  // spatio-temporal random effect - autocorrelated and scaled by ln_tau_epsilon
  for( int t=0; t<n_t; t++){
    if( t==0 ){
      jnll += SCALE( GMRF(Q), 1 / exp(ln_tau_epsilon) / pow( 1.0-pow(rhoE,2), 0.5 ) )( epsilon_st.col(t) );
    }else{
      jnll += SCALE( GMRF(Q), 1 / exp(ln_tau_epsilon) )( epsilon_st.col(t) - rhoE*epsilon_st.col(t-1) );
    }
  }


  
  // SDM projections
  // vector<Type> omega_g( n_g );
  // omega_g = A_gs * omega_s;
  matrix<Type> epsilon_gt( n_g, n_t );
  epsilon_gt = A_gs * epsilon_st;
  
  // Probability of data conditional on random effects
  
  // SDM species A
  // vector<Type> omega_i( n_i );
  // omega_i = A_is * omega_s;
  matrix<Type> epsilon_it( n_i, n_t );
  epsilon_it = A_is * epsilon_st;

  vector<Type> dhat_i( n_i );
  dhat_i.setZero();


  for( int i=0; i<D_i.size(); i++){
    dhat_i(i) = exp( beta_t(year_factor_t_i(t_i(i))) + epsilon_it(i, t_i(i))  +
      // beta_temp*temp_i(i) + beta_dist*dist_i(i));
      eta_smooth_i(i));
    jnll -= dtweedie( D_i(i), dhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true ) * weights_i(i);
  }
  
  
  
  // SDM: modeled density in the projection grid for species A
  
  // Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
  array<Type> ln_d_gt( n_g, n_t );
  
  // estimate the smooth effects in the projection grid
  // smoothers for projection grid
  vector<Type> proj_smooth_i(n_g*n_t);
  proj_smooth_i.setZero();
      for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
        vector<Type> beta_s(proj_Zs(s).cols());
        beta_s.setZero();
        for (int j = 0; j < beta_s.size(); j++) {
          beta_s(j) = b_smooth(b_smooth_start(s) + j);
        }
        proj_smooth_i += proj_Zs(s) * beta_s;
      }
      proj_smooth_i += proj_Xs * bs;
  
  // add smoothed effects, fixed effects, and random effects
  for( int t=0; t<n_t; t++){
    if( year_factor_t_i(t)==-999 ){
      for( int g=0; g<n_g; g++){
        ln_d_gt(g,t) = epsilon_gt(g,t) + proj_smooth_i(g + n_g * t) ;
      }  
    }else{
      for( int g=0; g<n_g; g++){
        ln_d_gt(g,t) = beta_t(year_factor_t_i(t)) + epsilon_gt(g,t) + proj_smooth_i(g + n_g * t) ;
      }
    }
  }
  
  // Reporting: SDM
  REPORT( dhat_i );
  REPORT( beta_t );
  // REPORT( beta_temp );
  // REPORT( beta_dist );
  // REPORT( ln_tau_omega );
  REPORT( ln_tau_epsilon );
  REPORT( ln_kappa );
  REPORT( ln_phi );
  REPORT( finv_power );
  // REPORT(SigmaO);
  REPORT(SigmaE);
  REPORT(Range);
  REPORT(ln_d_gt);
  // REPORT( H );
  REPORT( logit_rhoE );
  
  return jnll;
}