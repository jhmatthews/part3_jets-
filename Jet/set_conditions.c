/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Variable jet routines.

  Functions for creating a variable jet.

  \author J. Matthews (james.matthews@physics.ox.ac.uk)
  \date   May 18, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "extras.h"
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
#define TIME_SERIES 1

static int init_jet;

void Init_Variable_Jet(char filename[LINELEN])
{
  char dummy[LINELEN];
  FILE *fptr;
  int n;
  double t, g;

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    print ("Could not open filename %s\n", filename);
    exit (0);
  }

  n = 0;
  while (fgets (dummy, LINELEN, fptr) != NULL)
  {
    if ((dummy[0] != '#'))
    {
      sscanf (dummy, "%le %le", &t, &g);
      /* convert from Myrs to seconds then to simulation unit time */
      // print ("%le %le\n", t, f);
      TimeSeries.time[n] = t * 1e6 * 86400.0 * 365.25 / UNIT_TIME;
      TimeSeries.gmm[n] = g;

      n++;
    }
  }
  TimeSeries.npts = n - 1;
  TimeSeries.locator = 0;
  init_jet = 1;
}

/* ///////////////////////////////////////////////////////////////////// */
double GetJetParams(double *var_array, double eta) 
/*! 
  Specify the properties inside the jet - can be used to set a variable jet.
*/
/* ///////////////////////////////////////////////////////////////////// */
{
  double v_j, rho, v_j_cubed, interpolated_gmm, time_physical;
  double power, width_physical, rho_physical; 
  double mach_number = g_inputParam[JET_MACH_NUMBER];
  char filename[LINELEN];
  /* improve this */
  sprintf (filename, "gmm_Q%.0f_sig%.1f_eta%.0f_seed%.0f.dat", log10(g_inputParam[Q0]), g_inputParam[SIGMA], log10(g_inputParam[ETA]), g_inputParam[RNG_SEED]);

  if (TIME_SERIES && init_jet == 0) {
    print ("Reading from %s\n", filename);
    Init_Variable_Jet(filename);
    print ("Initialised Variable Jet");
  }

  /* get density in simulation units */
  rho = 1.0 / eta;
  width_physical = g_inputParam[JET_WIDTH] * UNIT_LENGTH;
  rho_physical = rho * UNIT_DENSITY;

  if (TIME_SERIES)
  {
    /* advance the locator */
    time_physical = g_time * UNIT_TIME;

    while (g_time > TimeSeries.time[TimeSeries.locator] && TimeSeries.locator < TimeSeries.npts)
    {
      TimeSeries.locator++;
    }
    
    /* check if we've run out of points */
    if (TimeSeries.locator == TimeSeries.npts)
    {
      print ("Ran out of Time series!!");
      print ("Time %8.4e %8.4e\n", time_physical, TimeSeries.time[TimeSeries.locator]);
      print ("%8.4e %8.4e %d\n", g_time, TimeSeries.time[TimeSeries.locator], TimeSeries.locator);
      exit (0);
    }

    interpolated_gmm = TimeSeries.gmm[TimeSeries.locator] * (TimeSeries.time[TimeSeries.locator + 1]-g_time);
    interpolated_gmm += TimeSeries.gmm[TimeSeries.locator + 1] * (g_time - TimeSeries.time[TimeSeries.locator]);
    interpolated_gmm /= TimeSeries.time[TimeSeries.locator + 1] - TimeSeries.time[TimeSeries.locator];
    interpolated_gmm = TimeSeries.gmm[TimeSeries.locator];

    /* this is in physical units */
    v_j = sqrt(1.0 - (1.0 / interpolated_gmm / interpolated_gmm));

    /* copy to array */
    var_array[RHO] = rho;
    var_array[VX2] = v_j;

    /* set pressure by jet mach number */
    var_array[PRS] = v_j * v_j / mach_number / mach_number / (5.0/3.0) * rho;

    // print ("JM: %8.4e %8.4e %8.4e %8.4e\n", power, v_j, g_time, time_physical);
    //print ("Velocity is %8.4e power %8.4e %8.4e %d\n", v_j, power, TimeSeries.flux[TimeSeries.locator], TimeSeries.locator);
    //print ("Time %8.4e %8.4e\n", time_physical, TimeSeries.time[TimeSeries.locator]);
  }  
  else 
  {
    power = 1e45;
    //v_j = g_inputParam[V_OVER_C] * CONST_c / UNIT_VELOCITY;  
    v_j_cubed = power / rho_physical / width_physical / width_physical / CONST_PI;
    v_j = pow(v_j_cubed, 1./3.);

    /* copy to array */
    var_array[RHO] = rho;
    var_array[VX2] = v_j;

    /* set pressure by jet mach number */
    var_array[PRS] = v_j * v_j / mach_number / mach_number / (5.0/3.0) * rho;
  }

  return v_j;
}

/* ///////////////////////////////////////////////////////////////////// */
double GetAmbientDensity(double x1, double x2, double x3, double beta, double r_c, double rho0) 
/*! 
  Get the ambient density at a location x1, x2 according to an isothermal
  King profile of the form
  \f[
      \rho(r) = \rho_0 \left[1 + \left(\frac{r}{r_c}\right)^2 \right]^{-3 \beta/2}
  \f]

  where \f$r_c\f$ is the core radius, \f$\rho_0\f$ is the base density 
  and \f$\beta\f$

  \author J. Matthews (james.matthews@physics.ox.ac.uk)
  \date   September 27, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
{
  double r, rho, exponent;
  if ((beta < 0) || (beta > 2)) {
      print ("!! Error: beta must be >0,<2.\n");
      exit(0);
  }
  #if GEOMETRY == CARTESIAN
    r = sqrt(x1*x1 + x2*x2 + x3*x3); /* spherical radius in cart. coords */
  #elif GEOMETRY == CYLINDRICAL
    r = sqrt(x1*x1 + x2*x2); /* spherical radius in cyl. coords */
  #elif GEOMETRY == SPHERICAL
    r = x1; /* spherical radius in sph. coords */
  #endif

  exponent = -3.0 * beta / 2.0;
  rho = 1.0 + pow(r/r_c, 2.0);
  rho = pow(rho, exponent);
  rho *= rho0;

  return rho;
}