// Model 3: Main effects for year and space, with the independent interaction of space and year
// An intercept that is independent by year that represents the median density in the year
// Spatial random fields
// Spatiotemporal random fields, that are independent by year


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
  DATA_INTEGER(n_t); // number of years in the dataset
  DATA_VECTOR(weights_i); // optional weights - used for cAIC
  DATA_VECTOR(temp_i);  // temperature for measurement i
  DATA_VECTOR(dist_i);  // distance from shore for measurement i
  
  DATA_STRUCT(Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(proj_Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_MATRIX(proj_Xs); // smoother linear effect matrix
  
  // Projection matrices for species A
  // A_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
  DATA_SPARSE_MATRIX(A_is);
  // A_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
  DATA_SPARSE_MATRIX(A_gs);

  DATA_MATRIX(temp_gt); // temperature at each location in each year
  DATA_VECTOR(dist_g); // distance from shore at each location
  
  // Data for smoothers
  DATA_IVECTOR(b_smooth_start);
  
  
  // Parameters
  
  // Parameters for SDM for species A
  // PARAMETER_VECTOR( beta_t );
  // PARAMETER( beta_temp );
  // PARAMETER( beta_dist );
  //PARAMETER( ln_tau_omega );
  //PARAMETER( ln_tau_epsilon );
  //PARAMETER( ln_kappa );
  PARAMETER( ln_phi ); // phi term in tweedie
  PARAMETER( finv_power ); // power parameter in tweedie
  
  // Parameters for smooth
  // PARAMETER_ARRAY(bs); // smoother linear effects
  // PARAMETER_ARRAY(b_smooth);  // P-spline smooth parameters
  // PARAMETER_ARRAY(ln_smooth_sigma);  // variances of spline REs if included
  PARAMETER_VECTOR(bs); // smoother linear effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  
  // Objective function
  Type jnll = 0;
  int n_i = A_is.rows();
  int n_g = A_gs.rows();

  vector<Type> dhat_i( n_i );
  dhat_i.setZero();
  
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


  for( int i=0; i<D_i.size(); i++){
    // dhat_i(i) = exp( beta_t(t_i(i)) + // + omega_i(i) + epsilon_it(i, t_i(i)) +
    //   beta_temp*temp_i(i) + beta_dist*dist_i(i));
    dhat_i(i) = exp( eta_smooth_i(i));
    jnll -= dtweedie( D_i(i), dhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true ) * weights_i(i);
  }
  
  
  
  // SDM: modeled density in the projection grid for species A
  // Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
  // array<Type> ln_d_gt( n_g, n_t );
  // 
  // // smoothers for projection grid
  // vector<Type> proj_smooth_i(n_g);
  // proj_smooth_i.setZero();
  //     for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
  //       vector<Type> beta_s(proj_Zs(s).cols());
  //       beta_s.setZero();
  //       for (int j = 0; j < beta_s.size(); j++) {
  //         beta_s(j) = b_smooth(b_smooth_start(s) + j);
  //       }
  //       proj_smooth_i += proj_Zs(s) * vector<Type>(beta_s);
  //     }
  //     proj_smooth_i += proj_Xs * vector<Type>(bs.col(1));
  //     proj_fe += proj_smooth_i;
  
  
  
  // for( int t=0; t<n_t; t++){
  //   for( int g=0; g<n_g; g++){
  //     ln_d_gt(g,t) = beta_t(t) + beta_temp*temp_gt(g,t) + beta_dist*dist_g(g) ;
  //   }
  // }
  // 
  // Reporting: SDM
  REPORT( dhat_i );
  // REPORT( beta_t );
  // REPORT( beta_temp );
  // REPORT( beta_dist );
  //REPORT( ln_tau_omega );
  //REPORT( ln_tau_epsilon );
  //REPORT( ln_kappa );
  //REPORT( ln_phi );
  REPORT( finv_power );
  //REPORT(SigmaO);
  //REPORT(SigmaE);
  // REPORT(Range);
  // REPORT(ln_d_gt);
  //REPORT( H );
  
  return jnll;
}