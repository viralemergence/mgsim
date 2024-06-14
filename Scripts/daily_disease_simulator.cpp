#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Helper Function for Breeding Season SIRI Simulation
//'
//' This function is an internal one for my Frankenmodel (epizootic + metaRange)
//' which performs daily simulations of disease dynamics and demography for all
//' populations at one breeding season timestep.
//'
//' @name daily_siri_summer
//'
//' @param Sj_abundance A matrix of susceptible juvenile abundances.
//' @param Sa_abundance A matrix of susceptible adult abundances.
//' @param I1j_abundance A matrix of juveniles infected for the first time.
//' @param I1a_abundance A matrix of adults infected for the first time.
//' @param Rj_abundance A matrix of recovered juveniles.
//' @param Ra_abundance A matrix of recovered adults.
//' @param I2j_abundance A matrix of juveniles infected for the second time.
//' @param I2a_abundance A matrix of adults infected for the second time.
//' @param fecundity Seasonal fecundity of adults.
//' @param transmission_Sj_summer Transmission rate from susceptible juveniles to first-time infected juveniles in summer.
//' A single numeric value.
//' @param transmission_Sa_summer Transmission rate from susceptible adults to first-time infected adults in summer.
//' A single numeric value.
//' @param transmission_Rj_summer Transmission rate from recovered juveniles to second-time infected juveniles in summer.
//' A single numeric value.
//' @param transmission_Ra_summer Transmission rate from recovered infected adults to second-time infected adults in summer.
//' A single numeric value.
//' @param recovery_I1j_summer Recovery rate from first-time infected juveniles to recovered juveniles in summer.
//' A single numeric value.
//' @param recovery_I1a_summer Recovery rate from first-time infected adults to recovered adults in summer.
//' A single numeric value.
//' @param recovery_I2j_summer Recovery rate from second-time infected juveniles to recovered juveniles in summer.
//' A single numeric value.
//' @param recovery_I2a_summer Recovery rate from second-time infected adults to recovered adults in summer.
//' A single numeric value.
//' @param mortality_Sj_summer Mortality rate of susceptible juveniles in summer.
//' A single numeric value.
//' @param mortality_Sa_summer Mortality rate of susceptible adults in summer.
//' A single numeric value.
//' @param mortality_I1j_summer Mortality rate of first-time infected juveniles in summer.
//' A single numeric value.
//' @param mortality_I1a_summer Mortality rate of first-time infected adults in summer.
//' A single numeric value.
//' @param mortality_Rj_summer Mortality rate of recovered juveniles in summer.
//' A single numeric value.
//' @param mortality_Ra_summer Mortality rate of recovered adults in summer.
//' A single numeric value.
//' @param mortality_I2j_summer Mortality rate of second-time infected juveniles in summer.
//' A single numeric value.
//' @param mortality_I2a_summer Mortality rate of second-time infected adults in summer.
//' A single numeric value.
//' @param season_length A vector of season lengths in days.
//' @param abundance_threshold A vector of quasi-extinction thresholds below which a
//' population becomes extinct.
//' @param density_max The maximum population density for a population.
//' @param habitat_suitability A numeric vector that indicates the habitat suitabilities
//' for the populations.
//' @return A matrix of 8 rows by N populations, where N is the length of the
//' input population matrices.
// [[Rcpp::export]]
arma::mat daily_siri_summer(
  arma::mat Sj_abundance,
  arma::mat Sa_abundance,
  arma::mat I1j_abundance,
  arma::mat I1a_abundance,
  arma::mat Rj_abundance,
  arma::mat Ra_abundance,
  arma::mat I2j_abundance,
  arma::mat I2a_abundance,
  double fecundity,
  double transmission_Sj_summer,
  double transmission_Sa_summer,
  double transmission_Rj_summer,
  double transmission_Ra_summer,
  double recovery_I1j_summer,
  double recovery_I1a_summer,
  double recovery_I2j_summer,
  double recovery_I2a_summer,
  double mortality_Sj_summer,
  double mortality_Sa_summer,
  double mortality_I1j_summer,
  double mortality_I1a_summer,
  double mortality_Rj_summer,
  double mortality_Ra_summer,
  double mortality_I2j_summer,
  double mortality_I2a_summer,
  arma::vec season_length,
  double abundance_threshold,
  double density_max,
  arma::vec habitat_suitability
) {
  // Prep the data
  int n_stages = 8;
  int n_pops = Sj_abundance.n_elem;
  arma::mat pop_matrix(n_stages, n_pops);
  arma::vec carrying_capacity = arma::round(density_max * habitat_suitability);
  // Find and handle non-finite values in season_length
  season_length.elem( arma::find_nonfinite(season_length) ).zeros();
  carrying_capacity.elem( arma::find_nonfinite(carrying_capacity) ).zeros();

  for (int p = 0; p < n_pops; p++) {
    if (season_length(p) == 0) {
      continue;
    }
    if (carrying_capacity(p) == 0) {
      continue;
    }
    // Bypass if the population size is 0
    int population_size = Sj_abundance(p) + Sa_abundance(p) + I1j_abundance(p) + I1a_abundance(p) + Rj_abundance(p) + Ra_abundance(p) + I2j_abundance(p) + I2a_abundance(p);
    if (population_size == 0) {
      continue;
    }
    // Prep the parameters
    double birth_rate = fecundity / season_length(p);
    arma::vec mortality(n_stages);
    arma::vec dd_mortality(n_stages);
    mortality(0) = mortality_Sj_summer / season_length(p);
    mortality(1) = mortality_Sa_summer / season_length(p);
    mortality(2) = mortality_I1j_summer;
    mortality(3) = mortality_I1a_summer;
    mortality(4) = mortality_Rj_summer / season_length(p);
    mortality(5) = mortality_Ra_summer / season_length(p);
    mortality(6) = mortality_I2j_summer;
    mortality(7) = mortality_I2a_summer;

    arma::mat state(n_stages, season_length(p) + 1);
    state(0, 0) = Sj_abundance(p);
    state(1, 0) = Sa_abundance(p);
    state(2, 0) = I1j_abundance(p);
    state(3, 0) = I1a_abundance(p);
    state(4, 0) = Rj_abundance(p);
    state(5, 0) = Ra_abundance(p);
    state(6, 0) = I2j_abundance(p);
    state(7, 0) = I2a_abundance(p);

    for (int t = 0; t < season_length(p); t++) {
      // Unpack states
      double Sj = state(0, t);
      double Sa = state(1, t);
      double I1j = state(2, t);
      double I1a = state(3, t);
      double Rj = state(4, t);
      double Ra = state(5, t);
      double I2j = state(6, t);
      double I2a = state(7, t);

      double N = std::min(arma::accu(state.col(t)), carrying_capacity(p));

      if (N < abundance_threshold) {
        state.cols(t + 1, season_length(p)).zeros();
        break;
      }

      for (int i = 0; i < n_stages; i++) {
        dd_mortality(i) = std::min((1.0 + N / carrying_capacity(p)) * mortality(i), 1.0);
      }

      double adults = Sa + I1a + Ra + I2a;
      double birth_rate_adj = birth_rate * (1.0 - N / carrying_capacity(p));

      double infected = I1j + I1a + I2j + I2a;

      double new_juv = Rcpp::rpois(1, adults * birth_rate_adj)[0];

      bool hasInfections = infected > 0;

      if (!hasInfections) {
        // Only update Sj and Sa, skipping all infection and recovery calculations
        state(0, t + 1) = Sj + new_juv - Rcpp::rbinom(1, Sj + new_juv, dd_mortality(0))[0];
        state(1, t + 1) = Sa - Rcpp::rbinom(1, Sa, dd_mortality(1))[0];
        state(2, t + 1) = I1j;
        state(3, t + 1) = I1a;
        state(4, t + 1) = Rj - Rcpp::rbinom(1, Rj, dd_mortality(4))[0];
        state(5, t + 1) = Ra - Rcpp::rbinom(1, Ra, dd_mortality(5))[0];
        state(6, t + 1) = I2j;
        state(7, t + 1) = I2a;
      } else {
        double infection1_juv = std::min(Rcpp::rbinom(1, Sj * infected, transmission_Sj_summer)[0], Sj);
        double infection1_adult = std::min(Rcpp::rbinom(1, Sa * infected, transmission_Sa_summer)[0], Sa);
        double susceptible_adult_death = Rcpp::rbinom(1, Sa - infection1_adult, dd_mortality(1))[0];
        double susceptible_juvenile_death = Rcpp::rbinom(1, Sj + new_juv - infection1_juv, dd_mortality(0))[0];
        double infected1_juvenile_death = Rcpp::rbinom(1, I1j + infection1_juv, dd_mortality(2))[0];
        double infected1_adult_death = Rcpp::rbinom(1, I1a + infection1_adult, dd_mortality(3))[0];
        double recovery1_juv = std::min(Rcpp::rbinom(1, I1j + infection1_juv - infected1_juvenile_death, recovery_I1j_summer)[0], I1j + infection1_juv - infected1_juvenile_death);
        double recovery1_adult = std::min(Rcpp::rbinom(1, I1a + infection1_adult - infected1_adult_death, recovery_I1a_summer)[0], I1a + infection1_adult - infected1_adult_death);
        double infection2_juv = std::min(Rcpp::rbinom(1, (Rj + recovery1_juv) * infected, transmission_Rj_summer)[0], Rj + recovery1_juv);
        double infection2_adult = std::min(Rcpp::rbinom(1, (Ra + recovery1_adult) * infected, transmission_Ra_summer)[0], Ra + recovery1_adult);
        double infected2_juvenile_death = Rcpp::rbinom(1, I2j + infection2_juv, dd_mortality(6))[0];
        double infected2_adult_death = Rcpp::rbinom(1, I2a + infection2_adult, dd_mortality(7))[0];
        double recovery2_juv = std::min(Rcpp::rbinom(1, I2j + infection2_juv - infected2_juvenile_death, recovery_I2j_summer)[0], I2j + infection2_juv - infected2_juvenile_death);
        double recovery2_adult = std::min(Rcpp::rbinom(1, I2a + infection2_adult - infected2_adult_death, recovery_I2a_summer)[0], I2a + infection2_adult - infected2_adult_death);
        double recovered_juvenile_death = Rcpp::rbinom(1, Rj + recovery1_juv + recovery2_juv - infection2_juv, dd_mortality(4))[0];
        double recovered_adult_death = Rcpp::rbinom(1, Ra + recovery1_adult + recovery2_adult - infection2_adult, dd_mortality(5))[0];

        // Update state for the next time step
        state(0, t + 1) = Sj + new_juv - infection1_juv - susceptible_juvenile_death;
        state(1, t + 1) = Sa - infection1_adult - susceptible_adult_death;
        state(2, t + 1) = I1j + infection1_juv - recovery1_juv - infected1_juvenile_death;
        state(3, t + 1) = I1a + infection1_adult - recovery1_adult - infected1_adult_death;
        state(4, t + 1) = Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death;
        state(5, t + 1) = Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death;
        state(6, t + 1) = I2j + infection2_juv - recovery2_juv - infected2_juvenile_death;
        state(7, t + 1) = I2a + infection2_adult - recovery2_adult - infected2_adult_death;
      }

    }
    pop_matrix.col(p) = state.col(season_length(p));
  }
  return pop_matrix;
}

//' Helper Function for Non-Breeding Season SIRI Simulation
 //'
 //' This function is an internal one for my Frankenmodel (epizootic + metaRange)
 //' which performs daily simulations of disease dynamics and demography for all
 //' populations at one non-breeding season timestep.
 //'
 //' @name daily_siri_winter
 //'
 //' @param Sj_abundance A matrix of susceptible juvenile abundances.
 //' @param Sa_abundance A matrix of susceptible adult abundances.
 //' @param I1j_abundance A matrix of juveniles infected for the first time.
 //' @param I1a_abundance A matrix of adults infected for the first time.
 //' @param Rj_abundance A matrix of recovered juveniles.
 //' @param Ra_abundance A matrix of recovered adults.
 //' @param I2j_abundance A matrix of juveniles infected for the second time.
 //' @param I2a_abundance A matrix of adults infected for the second time.
 //' @param transmission_Sj_winter Transmission rate from susceptible juveniles to first-time infected juveniles in winter.
 //' A single numeric value.
 //' @param transmission_Sa_winter Transmission rate from susceptible adults to first-time infected adults in winter.
 //' A single numeric value.
 //' @param transmission_Rj_winter Transmission rate from recovered juveniles to second-time infected juveniles in winter.
 //' A single numeric value.
 //' @param transmission_Ra_winter Transmission rate from recovered infected adults to second-time infected adults in winter.
 //' A single numeric value.
 //' @param recovery_I1j_winter Recovery rate from first-time infected juveniles to recovered juveniles in winter.
 //' A single numeric value.
 //' @param recovery_I1a_winter Recovery rate from first-time infected adults to recovered adults in winter.
 //' A single numeric value.
 //' @param recovery_I2j_winter Recovery rate from second-time infected juveniles to recovered juveniles in winter.
 //' A single numeric value.
 //' @param recovery_I2a_winter Recovery rate from second-time infected adults to recovered adults in winter.
 //' A single numeric value.
 //' @param mortality_Sj_winter Mortality rate of susceptible juveniles in winter.
 //' A single numeric value.
 //' @param mortality_Sa_winter Mortality rate of susceptible adults in winter.
 //' A single numeric value.
 //' @param mortality_I1j_winter Mortality rate of first-time infected juveniles in winter.
 //' A single numeric value.
 //' @param mortality_I1a_winter Mortality rate of first-time infected adults in winter.
 //' A single numeric value.
 //' @param mortality_Rj_winter Mortality rate of recovered juveniles in winter.
 //' A single numeric value.
 //' @param mortality_Ra_winter Mortality rate of recovered adults in winter.
 //' A single numeric value.
 //' @param mortality_I2j_winter Mortality rate of second-time infected juveniles in winter.
 //' A single numeric value.
 //' @param mortality_I2a_winter Mortality rate of second-time infected adults in winter.
 //' A single numeric value.
 //' @param season_length A vector of season lengths in days.
 //' @param abundance_threshold A vector of quasi-extinction thresholds below which a
 //' population becomes extinct.
 //' @param density_max The maximum population density for a population.
 //' @param habitat_suitability A numeric vector that indicates the habitat suitabilities
 //' for the populations.
 //' @return A matrix of 8 rows by N populations, where N is the length of the
 //' input population matrices.
 // [[Rcpp::export]]
 arma::mat daily_siri_winter(
     arma::mat Sj_abundance,
     arma::mat Sa_abundance,
     arma::mat I1j_abundance,
     arma::mat I1a_abundance,
     arma::mat Rj_abundance,
     arma::mat Ra_abundance,
     arma::mat I2j_abundance,
     arma::mat I2a_abundance,
     double transmission_Sj_winter,
     double transmission_Sa_winter,
     double transmission_Rj_winter,
     double transmission_Ra_winter,
     double recovery_I1j_winter,
     double recovery_I1a_winter,
     double recovery_I2j_winter,
     double recovery_I2a_winter,
     double mortality_Sj_winter,
     double mortality_Sa_winter,
     double mortality_I1j_winter,
     double mortality_I1a_winter,
     double mortality_Rj_winter,
     double mortality_Ra_winter,
     double mortality_I2j_winter,
     double mortality_I2a_winter,
     arma::vec season_length,
     double abundance_threshold,
     double density_max,
     arma::vec habitat_suitability
 ) {
   // Prep the data
   int n_stages = 8;
   int n_pops = Sj_abundance.n_elem;
   arma::mat pop_matrix(n_stages, n_pops);
   arma::vec carrying_capacity = density_max * habitat_suitability;
   arma::vec winter_length = 365 - season_length;  // Calculate the winter length

   // Find and handle non-finite values in season_length and carrying_capacity
   winter_length.elem(arma::find_nonfinite(winter_length)).zeros();
   carrying_capacity.elem(arma::find_nonfinite(carrying_capacity)).zeros();

   for (int p = 0; p < n_pops; p++) {
     if (winter_length(p) == 0) {
       continue;
     }
     if (carrying_capacity(p) == 0) {
       continue;
     }
     // Bypass if the population size is 0
     double population_size = Sj_abundance(p) + Sa_abundance(p) + I1j_abundance(p) + I1a_abundance(p) + Rj_abundance(p) + Ra_abundance(p) + I2j_abundance(p) + I2a_abundance(p);
     if (population_size == 0) {
       continue;
     }
     // Prep the parameters
     arma::vec mortality(n_stages);
     mortality(0) = mortality_Sj_winter / winter_length(p);
     mortality(1) = mortality_Sa_winter / winter_length(p);
     mortality(2) = mortality_I1j_winter;
     mortality(3) = mortality_I1a_winter;
     mortality(4) = mortality_Rj_winter / winter_length(p);
     mortality(5) = mortality_Ra_winter / winter_length(p);
     mortality(6) = mortality_I2j_winter;
     mortality(7) = mortality_I2a_winter;
     arma::vec dd_mortality = mortality;

     arma::mat state(n_stages, winter_length(p) + 1);
     state(0, 0) = Sj_abundance(p);
     state(1, 0) = Sa_abundance(p);
     state(2, 0) = I1j_abundance(p);
     state(3, 0) = I1a_abundance(p);
     state(4, 0) = Rj_abundance(p);
     state(5, 0) = Ra_abundance(p);
     state(6, 0) = I2j_abundance(p);
     state(7, 0) = I2a_abundance(p);

     for (int t = 0; t < winter_length(p); t++) {
       // Unpack states
       double Sj = state(0, t);
       double Sa = state(1, t);
       double I1j = state(2, t);
       double I1a = state(3, t);
       double Rj = state(4, t);
       double Ra = state(5, t);
       double I2j = state(6, t);
       double I2a = state(7, t);

       double N = std::min(arma::accu(state.col(t)), carrying_capacity(p));

       if (N < abundance_threshold) {
         state.cols(t + 1, winter_length(p)).zeros();
         break;
       }

       for (int i = 0; i < n_stages; i++) {
         dd_mortality(i) = std::min((1.0 + N / carrying_capacity(p)) * mortality(i), 1.0);
       }

       double infected = I1j + I1a + I2j + I2a;

       bool hasInfections = infected > 0;

       if (!hasInfections) {
         // Only update Sj and Sa, skipping all infection and recovery calculations
         state(0, t + 1) = Sj - Rcpp::rbinom(1, Sj, dd_mortality(0))[0];
         state(1, t + 1) = Sa - Rcpp::rbinom(1, Sa, dd_mortality(1))[0];
         state(2, t + 1) = I1j;
         state(3, t + 1) = I1a;
         state(4, t + 1) = Rj - Rcpp::rbinom(1, Rj, dd_mortality(4))[0];
         state(5, t + 1) = Ra - Rcpp::rbinom(1, Ra, dd_mortality(5))[0];
         state(6, t + 1) = I2j;
         state(7, t + 1) = I2a;
       } else {
         double infection1_juv = std::min(Rcpp::rbinom(1, Sj * infected, transmission_Sj_winter)[0], Sj);
         double infection1_adult = std::min(Rcpp::rbinom(1, Sa * infected, transmission_Sa_winter)[0], Sa);
         double susceptible_adult_death = Rcpp::rbinom(1, Sa - infection1_adult, dd_mortality(1))[0];
         double susceptible_juvenile_death = Rcpp::rbinom(1, Sj - infection1_juv, dd_mortality(0))[0];
         double infected1_juvenile_death = Rcpp::rbinom(1, I1j + infection1_juv, dd_mortality(2))[0];
         double infected1_adult_death = Rcpp::rbinom(1, I1a + infection1_adult, dd_mortality(3))[0];
         double recovery1_juv = std::min(Rcpp::rbinom(1, I1j + infection1_juv - infected1_juvenile_death, recovery_I1j_winter)[0], I1j + infection1_juv - infected1_juvenile_death);
         double recovery1_adult = std::min(Rcpp::rbinom(1, I1a + infection1_adult - infected1_adult_death, recovery_I1a_winter)[0], I1a + infection1_adult - infected1_adult_death);
         double infection2_juv = std::min(Rcpp::rbinom(1, (Rj + recovery1_juv) * infected, transmission_Rj_winter)[0], Rj + recovery1_juv);
         double infection2_adult = std::min(Rcpp::rbinom(1, (Ra + recovery1_adult) * infected, transmission_Ra_winter)[0], Ra + recovery1_adult);
         double infected2_juvenile_death = Rcpp::rbinom(1, I2j + infection2_juv, dd_mortality(6))[0];
         double infected2_adult_death = Rcpp::rbinom(1, I2a + infection2_adult, dd_mortality(7))[0];
         double recovery2_juv = std::min(Rcpp::rbinom(1, I2j + infection2_juv - infected2_juvenile_death, recovery_I2j_winter)[0], I2j + infection2_juv - infected2_juvenile_death);
         double recovery2_adult = std::min(Rcpp::rbinom(1, I2a + infection2_adult - infected2_adult_death, recovery_I2a_winter)[0], I2a + infection2_adult - infected2_adult_death);
         double recovered_juvenile_death = Rcpp::rbinom(1, Rj + recovery1_juv + recovery2_juv - infection2_juv, dd_mortality(4))[0];
         double recovered_adult_death = Rcpp::rbinom(1, Ra + recovery1_adult + recovery2_adult - infection2_adult, dd_mortality(5))[0];

         // Update state for the next time step
         state(0, t + 1) = Sj - infection1_juv - susceptible_juvenile_death;
         state(1, t + 1) = Sa - infection1_adult - susceptible_adult_death;
         state(2, t + 1) = I1j + infection1_juv - recovery1_juv - infected1_juvenile_death;
         state(3, t + 1) = I1a + infection1_adult - recovery1_adult - infected1_adult_death;
         state(4, t + 1) = Rj + recovery1_juv + recovery2_juv - infection2_juv - recovered_juvenile_death;
         state(5, t + 1) = Ra + recovery1_adult + recovery2_adult - infection2_adult - recovered_adult_death;
         state(6, t + 1) = I2j + infection2_juv - recovery2_juv - infected2_juvenile_death;
         state(7, t + 1) = I2a + infection2_adult - recovery2_adult - infected2_adult_death;
       }
     }

     pop_matrix.col(p) = state.col(winter_length(p));
   }
   return pop_matrix;
 }


