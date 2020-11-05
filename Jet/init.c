/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Hydrodynamic jet propagation in 3D

  This problem considers the propagation of a hydrodynamic jet into a
  an ambient medium with constant pressure. 
  The ambient density in sim units is \f$ \rho_a(y) = \eta (y/a_0)^{-\beta} \f$ 
  where \f$\eta\f$ the maximum ambient/jet density ratio in the active phase.
  The jet inflow is set through a user-defined boundary condition at the 
  lower z-boundary. A simple top-hat injection nozzle is used.
  
  The basic configuration is defined in terms of the following parameters:

  - <tt>ETA</tt>:   density ratio between ambient and jet;
  - <tt>V_OVER_C</tt>:  jet velocity, \f$v_j/c\f$ in active phase;
  - <tt>BETA</tt>:  Density power law index, \f$\beta \f$;
  - <tt>SCALE_LENGTH</tt>:  density scale length, \f$a_0 \f$;

  defined in \c pluto.ini. The jet variability is set by

  - <tt>RISE_FRACTION</tt>:  rise fraction, \f$f_{R} \f$;
  - <tt>DUTY_FRACTION</tt>:  duty fraction, \f$f_{DC}\f$;
  - <tt>PERIOD</tt>:  Variability period \f$T \f$;
  - <tt>FQ_DENSITY</tt>:  quiescent density (fraction of $\rho_j$);
  - <tt>FQ_VELOCITY</tt>: quiescent velocity (fraction of $v_j$).

  As shown in the following figure. The jet variability is defined using what is fundamentally a square wave
  with a sinusoidal rise and fall period. This parameterisation has sufficient
  flexibility to create a completely square wave (by setting RISE_FRACTION to zero)
  or a completely sinusoidal wave (by setting RISE_FRACTION to 0.5*DUTY_FRACTION).

 
  \image html pulses.png  "An example of a jet velocity pulse" width=50%

  The jet variability is defined using what is fundamentally a square wave
  with a sinusoidal rise and fall period. This parameterisation has sufficient
  flexibility to create a completely square wave (by setting RISE_FRACTION to zero)
  or a completely sinusoidal wave (by setting RISE_FRACTION to 0.5*DUTY_FRACTION).
  

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "extras.h"

static double Profile(double r, double w_jet, int nv);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  /* get some input parameters */
  double rho0 = 1.0;
  double beta = g_inputParam[BETA];
  double r_c = g_inputParam[CORE_RADIUS];
  double sound_speed = g_inputParam[CS_A];
  double rnd, drho;
  
  double scrh;

  /* normally use Taub-Matthews EOS which has cell by cell gamma */
  #if EOS == IDEAL
   g_gamma = 5.0/3.0;
  #endif
  
  v[RHO] = GetAmbientDensity(x1, x2, x3, beta, r_c, rho0);

  /* apply small density perturbation */
  rnd = (double)(rand()) / ((double)RAND_MAX + 1.0);
  drho = (-1 + 2.0 * rnd) * 1e-10 * v[RHO]; 
  v[RHO] += drho;

  /* initially at rest */
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  /* set ambient pressure such that sound speed is unit velocity */
  v[PRS] = sound_speed * sound_speed / (5./3.) * v[RHO];

  /* set a pressure floor - not clear if this is needed */
  g_smallPressure = v[PRS] / 500.0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}



void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  /* here you can save information on e.g. tracer particles */
}


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions in the lower boundary ghost
 *  zones.  The profile is top-hat: 
 *  \f[
 *     V_{ij} = \left\{\begin{array}{ll}
 *     V_{\rm jet} & \quad\mathrm{for}\quad r_i < 1 \\ \noalign{\medskip}
 *     \mathrm{Reflect}(V)  & \quad\mathrm{otherwise}
 *    \end{array}\right.
 *  \f]
 * where \f$ V_{\rm jet} = (\rho,v,p)_{\rm jet} = (1,M,1/\Gamma)\f$ and
 * \c M is the flow Mach number (the unit velocity is the jet sound speed, 
 * so \f$ v = M\f$).
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double vjet[NVAR], vout[NVAR], w_jet, r, eta;

  /* basic jet parameters from input file */
  w_jet = g_inputParam[JET_WIDTH];
  eta = g_inputParam[ETA];

  //CheckGrid(grid);
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  if (side == X2_BEG){
    GetJetParams(vjet, eta);
    X2_BEG_LOOP(k,j,i){
      VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG-j-1][i];
      vout[VX2] *= -1.0;

      r = x1[i];
      for (nv = 0; nv < NVAR; nv++)
      {
        d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,w_jet,nv);
        //if (r <= w_jet) d->Vc[nv][k][j][i] = vjet[nv];
        //else d->Vc[nv][k][j][i] = vout[nv];
      }

      /* assign a passive scalar */
      if (r <= w_jet)
        d->Vc[TRC][k][j][i] = 1.0;
    }
  }
}



#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * This is currently set up to operate as a stabilising force for the 
 * King profile with constant sound speed.
 *
 *********************************************************************** */
{
  double gs, rs;
  #if GEOMETRY == CARTESIAN
    rs = sqrt(x1*x1 + x2*x2 + x3*x3); /* spherical radius in cart. coords */
  #elif GEOMETRY == CYLINDRICAL
    rs = sqrt(x1*x1 + x2*x2); /* spherical radius in cyl. coords */
  #elif GEOMETRY == SPHERICAL
    rs = x1; /* spherical radius in sph. coords */
  #endif

  double beta = g_inputParam[BETA];
  double rho0 = 1.0;
  double rcore = g_inputParam[CORE_RADIUS];
  double cs = g_inputParam[CS_A];

  double exponent = (-3.0 * beta/2.0) - 1.0;
  double drhodr = 1.0 + pow(rs/rcore, 2.0);
  drhodr = pow(drhodr, exponent);
  drhodr *= rs * rho0 * (-3.0 * beta / rcore / rcore);

  double dPdr = drhodr * cs * cs / 1.666666667;
  double rho = GetAmbientDensity(x1, x2, x3, beta, rcore, rho0);
  gs = dPdr / rho; /* spherical gravity */

  #if GEOMETRY == CARTESIAN
    g[IDIR] = gs*x1/rs;
    g[JDIR] = gs*x2/rs;
    g[KDIR] = gs*x3/rs;
  #elif GEOMETRY == CYLINDRICAL
    g[IDIR] = gs*x1/rs;
    g[JDIR] = gs*x2/rs;
    g[KDIR] = 0.0;
  #elif GEOMETRY == SPHERICAL
    g[IDIR] = gs;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
  #endif
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 * __TODO__ include a Dark Matter potential here.
 * at the moment this deals with the CR pressure
 *
 *********************************************************************** */
{

  return 0.0;
}
#endif

double Profile(double r, double w_jet, int nv)
/* 
 * cosh smoothing profile 
 *
 *********************************************************************** */
{
  int xn = 14;
  double r0 = w_jet;

  if (nv == RHO) r0 = 1.1 * w_jet;

  #if GEOMETRY == SPHERICAL
   r0 = 5.0/180.0*CONST_PI;
  #endif
  return 1.0/cosh(pow(r/r0,xn));
}
