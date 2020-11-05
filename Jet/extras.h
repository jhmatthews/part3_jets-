/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief extra definitions used to modelling jets with pluto.

  \author J. Matthews (james.matthews@physics.ox.ac.uk)
  \date   May 15, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

/* routines for jet and ambient medium */
double GetJetParams(double *var_array, double eta);
double GetAmbientDensity(double x1, double x2, double x3, double beta, double r_c, double rho0);

#define LINELEN 132
void Init_Variable_Jet(char filename[LINELEN]);

#define NMAX_TIME 100000
struct timeseries
{
  int npts;
  double time[NMAX_TIME];
  double flux[NMAX_TIME];
  double t0; // offset in time
  int locator;
}
TimeSeries;