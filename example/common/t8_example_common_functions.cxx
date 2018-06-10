/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_example_common_functions.cxx Provide real valued functions
  * that are used in more than one example. */

#include <example/common/t8_example_common.h>
#include <t8_vec.h>
#include <complex.h>

T8_EXTERN_C_BEGIN ();

double
t8_levelset_sphere (const double x[3], double t, void *data)
{
  t8_levelset_sphere_data_t *ls_data = (t8_levelset_sphere_data_t *) data;

  T8_ASSERT (ls_data->radius > 0);
  return t8_vec_dist (x, ls_data->M) - ls_data->radius;
}

double
t8_scalar3d_constant_one (const double x[3], double t)
{
  return 1;
}

double
t8_scalar3d_constant_zero (const double x[3], double t)
{
  return 0;
}

double
t8_scalar3d_project_x (const double x[3], double t)
{
  return x[0];
}

double
t8_scalar3d_exp_distribution (const double x[3], double t)
{
  double              dummy, X;

  /* Get fractional part of t. t is thus periodically
   * mapped to the unit interval */
  t = modf (t, &dummy);
  X = x[0] - .5;
  return exp (-4 * X * X);
}

/* This function is =1 if 0.25 <= x <= 0.75 and 0 else */
double
t8_scalar3d_step_function (const double x[3], double t)
{
  return 0.25 <= x[0] && x[0] <= 0.75;
}

/* This function is =1 if 0.25 <= x <= 0.75,
 * it is 0 outside of 0.25-eps and 0.75+eps,
 * it interpolates linearly in between. */
double
t8_scalar3d_almost_step_function (const double x[3], double t)
{
  double              eps = 0.1;

  /* interpolate in [0.25-eps,0.25+eps] */
  if (0.25 - eps < x[0] && x[0] < 0.25) {
    return (x[0] - 0.25 + eps) / (eps);
  }
  /* interpolate in [0.75-eps,0.75+eps] */
  else if (0.75 < x[0] && x[0] < 0.75 + eps) {
    return 1 - (x[0] - 0.75) / (eps);
  }
  /* 1 inside [0.25,0.75], 0 outside */
  return 0.25 <= x[0] && x[0] <= 0.75;
}

double
t8_scalar3d_sinx (const double x[3], double t)
{
  return sin (2 * M_PI * x[0]) + 1;
}

double
t8_scalar3d_sinx_cosy (const double x[3], double t)
{
  return sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]);
}

double
t8_scalar3d_sinx_cosy_z (const double x[3], double t)
{
  return 10 * sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]) * x[3];
}

double
t8_scalar3d_sint (const double x[3], double t)
{
  return sin (2 * M_PI * t);
}

/* general level set function for a sphere with given midpoint and radius. */
static double
t8_scalar3d_sphere (const double x[3], double M[3], double radius)
{

  /* Compute M - x */
  t8_vec_axpy (x, M, -1);

  /* return |M-x| - radius */

  return t8_vec_norm (M) - radius;
}

double
t8_scalar3d_sphere_75_radius (const double x[3], double t)
{
  double              M[3] = { 0, 0, 0 };
  return t8_scalar3d_sphere (x, M, 0.75);
}

double
t8_scalar3d_sphere_05_midpoint_375_radius (const double x[3], double t)
{
  double              M[3] = { 0.5, 0.5, 0.5 };

  return t8_scalar3d_sphere (x, M, 0.375);
}

double
t8_scalar3d_sphere_03_midpoint_25_radius (const double x[3], double t)
{
  double              M[3] = { 0.3, 0.3, 0.3 };

  return t8_scalar3d_sphere (x, M, 0.25);
}

double
t8_scalar3d_sphere_05_0z_midpoint_375_radius (const double x[3], double t)
{
  double              M[3] = { 0.5, 0.5, 0 };

  return t8_scalar3d_sphere (x, M, 0.375);
}

void
t8_flow_constant_one_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = x_out[1] = x_out[2] = 1;
}

void
t8_flow_constant_one_x_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = 1;
  x_out[1] = x_out[2] = 0;
}

void
t8_flow_constant_one_xy_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = 1;
  x_out[1] = 0.8;
  x_out[2] = 0;
}

void
t8_flow_constant_one_xyz_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = 1;
  x_out[1] = 0.8;
  x_out[2] = 0.9;
}

void
t8_flow_rotation_2d (const double x_in[3], double t, double x_out[3])
{
  double              x = x_in[0], y = x_in[1];

  x -= 0.5;
  y -= 0.5;

  x_out[0] = y;
  x_out[1] = -x;
  x_out[2] = 0;

  t8_vec_ax (x_out, 2 * M_PI);
}

void
t8_flow_compressible (const double x_in[3], double t, double x_out[3])
{
  x_out[0] = (1. / 2 - x_in[0]);
  x_out[1] = 0;
  x_out[2] = 0;
}

 /* The following function is a incompressible flow on the unit cube.
  * It is constructed from any function f with f(0) = f(1) = 0.
  */

/* Function with f(0) = f(1) = 0 */
static double
t8_incomp_cube_f (double x)
{
  return 2 * (1. - x) * x;
}

/* The derivative of f */
static double
t8_incomp_cube_df (double x)
{

  return 2. - 4. * x;
}

static double
t8_incomp_cube_f_sin (double x)
{
  return sin (M_PI * x);
}

static double
t8_incomp_cube_df_sin (double x)
{
  return M_PI * cos (M_PI * x);
}

void
t8_flow_incomp_cube_flow (const double x[3], double t, double x_out[3])
{
  double              (*f) (double) = t8_incomp_cube_f_sin;
  double              (*df) (double) = t8_incomp_cube_df_sin;

  x_out[0] = f (x[0]) * (df (x[1]) - df (x[2]));
  x_out[1] = -1. * f (x[1]) * df (x[0]);
  x_out[2] = f (x[2]) * df (x[0]);

  t8_vec_ax (x_out, 1. / 2);
  /* We reverse the flow at time 0.5 */
  if (t > 0.5) {
    t8_vec_ax (x_out, -1);
  }
}

/* Convert the first two entries of a vector into 2D polar
 * coordinates. x = r cos(phi)
 *              y = r sin(phi)
 *
 * On output: polar[0] = r, polar[1] = phi
 */
static void
t8_flow_2d_polar_coords (const double x[3], double polar[2])
{
  polar[0] = sqrt (SC_SQR (x[0]) + SC_SQR (x[1]));
  polar[1] = atan2 (x[1], x[0]);
}

/* Convert a 2D vector from polar coordinates to cartesian
 * coordinates.
 * On input: polar[0] = r, polar[1] = phi
 *
 * On output: cart[0] = r cos(phi)
 *            cart[1] = r sin(phi)
 *
 */
static void
t8_flow_2d_cart_coords (const double polar_values[2],
                        const double polar_coords[2], double cart[2])
{
  cart[0] = cos (polar_coords[1]) * polar_values[0]
    - sin (polar_coords[1]) * polar_values[1];
  cart[1] = sin (polar_coords[1]) * polar_values[0]
    + cos (polar_coords[1]) * polar_values[1];
}

/* 2d flow around a circle with radius R = 1 and
 * constant inflow with x-speed U = 1. */
void
t8_flow_around_circle (const double x[3], double t, double x_out[3])
{
  double              polar[2];
  double              polar_speed[2];
  //const double        R = 0.15;
  const double        R = 1;

  /* transform [0,1] x [0,1] coordinates to [-0.5,0.5] x [-0.5,0.5] */
  //t8_vec_axb (x, x_out, 3, -1.6);
  //////////////
    t8_vec_axb (x, x_out, 1, 0.2);
  //////////////


  /* Set the z-coordinate to zero */
  x_out[2] = 0;

  if (t8_vec_norm (x_out) < R) {
    /* Set the velocity inside the circle to 0 */
    x_out[0] = x_out[1] = x_out[2] = 0;
    return;
  }
  /* Convert x,y coordinates to polar r,phi coordinates */
  t8_flow_2d_polar_coords (x_out, polar);
  /* Compute v_r (r,phi) = U (1-R^2/r^2)cos(phi) */
  polar_speed[0] = (1 - SC_SQR (R) / SC_SQR (polar[0])) * cos (polar[1]);
  /* Compute v_phi(r,phi) = -U (1+ R^2/r^2) sin (phi) */
  polar_speed[1] = -(1 + SC_SQR (R) / SC_SQR (polar[0])) * sin (polar[1]);
  /* Convert the coordinates back to cartesian */
  t8_flow_2d_cart_coords ((const double *) polar_speed, polar, x_out);
  /* Set the z-coordinate to zero */
  x_out[2] = 0;
}

/* Compute out = x*y for two complex numbers x and y */
static void
t8_flow_cmplx_mult (const double x[2], const double y[2], double out[2])
{
    out[0] = x[0]*y[0] - x[1]*y[1];
    out[1] = x[1]*y[0] + x[0]*y[1];
}


/* Compute out = x/y for two complex numbers x and y.
 * If y = 0 then out is set to 0 */
static void
t8_flow_cmplx_div (const double x[2], const double y[2], double out[2])
{
    double norm_ysq = y[0] * y[0] + y[1]*y[1];

    if (norm_ysq == 0) {
        /* prevent div by 0 */
        out[0] = out[1] = 0;
    }
    out[0] = (x[0]*y[0] + x[1]*y[1])/norm_ysq;
    out[1] = (x[1]*y[0] - x[0]*y[1])/norm_ysq;
}

/* Apply the transformation J(z) = u(z)/(1-lambda/z^2) (Joukowski)
 * to the flow  u around the circle to obtain a flow
 * around a NACA airfoil. Here we interpret the (x[0],x[1]) coordinates
 * as z = x + iy.
 * The x[3] coordinate is not modified.
 * See Modeling the Fluid Flow around Airfoils Using
 * Conformal Mapping - Nitin R. Kapania, Katherine Terracciano, Shannon Taylor
 */
void
t8_flow_circle_to_naca (const double x[3], double t, double x_out[3])
{
    int i;
    double norm_sq;
    const double lambda_sq = 4;
    double temp[3];
    double z_sq[2];
    double z_sqml[2];
    double u[2];
    double complex x_cmplx;
    double complex inverse_joukowski;
    double z[3];


    /* transform [0,1] x [0,1] coordinates to [-0.5,0.5] x [-0.5,0.5] */
    t8_vec_axb (x, temp, 3, -1.5);
    /* Compute the inverse Joukowski transform to map from the
     * Joukowski airfoil to the circle.
     * z = (w - sqrt(w^2-4))/2
     */
    x_cmplx = temp[0] + temp[1] * _Complex_I;
    inverse_joukowski = (x_cmplx - csqrt (x_cmplx*x_cmplx -4)) / 2;
    /* Convert complex  to vector */
    z[0] = creal (inverse_joukowski);
    z[1] = cimag (inverse_joukowski);
    z[2] = 0;

    /* Compute flow around the circle */
    t8_flow_around_circle (z, t, u);

    /* Since
     *    u / (1-l/z^2) = u / ((z^2 -l)/z^2) = u * (z^2/(z^2-l))
     * we compute u * (z^2/(z^2-l))
     */

    /* Compute z^2 */
    t8_flow_cmplx_mult (z, z, z_sq);

    ///////////////////////////
    z_sqml[0] = z_sq[0] - lambda_sq;
    z_sqml[1] = z_sq[1];

    /* Compute temp = z_sq/(z_sq-1) */
    t8_flow_cmplx_div (z_sq, z_sqml, temp);
    /* Compute x_out * temp */
    t8_flow_cmplx_mult (u, temp, x_out);
    /* Set third coordinate to 0 */
    x_out[2] = 0;
    return;
    /////////////////////////////////////////


    /* norm^2 of z^2 */
    norm_sq = z_sq[0] * z_sq[0] + z_sq[1] * z_sq[1];
    if (norm_sq == 0) {
        /* prevent div by 0 */
        x_out[0] = x_out[1] = 0;
        return;
    }

    /* temp = 1 - lambda/z^2 */
    temp[0] = 1- lambda_sq * z_sq[0]/norm_sq;
    temp[1] = -lambda_sq * z_sq[1]/norm_sq;
    /* norm of temp */
    norm_sq = temp[0] * temp[0] + temp[1] *temp[1];
    if (norm_sq == 0) {
        /* prevent div by 0 */
        x_out[0] = x_out[1] = 0;
        return;
    }
    /* Compute z_sq = x_out/(1-l/z^2) */
    z_sq[0] = (x_out[0]*temp[0]+x_out[1]*temp[1])/norm_sq;
    z_sq[1] = (x_out[1]*temp[0]-x_out[0]*temp[1])/norm_sq;

    x_out[0] = z_sq[0];
    x_out[1] = z_sq[1];
}


/* The following functions model a solution to the stokes equation on
 * a spherical shell. See
 * Analytical solution for viscous incompressible Stokes flow in a
 * spherical shell
 * by Cedric Thieulot
 */

static void
t8_flow_stokes_sphere_alpha_beta (double R_1, double R_2, double gamma, int m,
                                  double *alpha, double *beta)
{
  /* We define two constants alpha and beta */
  *alpha =
    gamma * (m + 1) * (pow (R_1, -3) - pow (R_2, -3)) / (pow (R_1, -m - 4) -
                                                         pow (R_2, -m - 4));
  *beta =
    -3 * gamma * (pow (R_1, m + 1) - pow (R_2, m + 1)) / (pow (R_1, m + 4) -
                                                          pow (R_2, m + 4));

}

/* A component of the flow that depends on the inner radius R_2, the outer radius R_1,
 * a constant gamma, and a control parameter m with m != -1, m != -4 */
static double
t8_flow_stokes_sphere_g_component (double radius, double alpha, double beta,
                                   double gamma, int m)
{
  T8_ASSERT (m != -1 && m != -4);

  return -2 / (radius * radius) * (-alpha / (m + 1) * pow (radius, -m - 1) +
                                   beta / 3 * pow (radius, 3) + gamma);
}

static double
t8_flow_stokes_sphere_f_component (double radius, double alpha, double beta,
                                   int m)
{
  return alpha * pow (radius, -m - 3) + beta * radius;
}

void
t8_flow_stokes_flow_sphere_shell (const double x[3], double t, double x_out[])
{
  double              radius;
  double              theta, phi;
  double              alpha, beta;
  double              vel_r;
  double              vel_theta;
  double              vel_phi;
  const double        r_1 = .5, r_2 = 1, gamma = 1, m = 3;

#if 1
  /* translate unit cube to cube centered around origin */
  ((double *) x)[0] -= 0.5;
  ((double *) x)[1] -= 0.5;
  ((double *) x)[2] -= 0.5;
  ((double *) x)[0] *= 2;
  ((double *) x)[1] *= 2;
  ((double *) x)[2] *= 2;
#endif
  /* Compute spherical coordinates */
  radius = t8_vec_norm (x);
  theta = acos (x[2] / radius);
#if 1
  /* Phi component, not used */
  phi = atan2 (x[1], x[0]);
#endif

  if (radius < r_1) {
    /* If there are points in the geometry that lie in the inside radius,
     * set the flow to zero. */
    x_out[0] = x_out[1] = x_out[2] = 0;
    return;
  }

  /* Compute alpha and beta */
  t8_flow_stokes_sphere_alpha_beta (r_1, r_2, gamma, m, &alpha, &beta);
  /* Compute radial velocity and theta velocity */
  vel_r =
    t8_flow_stokes_sphere_g_component (radius, alpha, beta, gamma,
                                       m) * cos (theta);
  vel_theta =
    t8_flow_stokes_sphere_f_component (radius, alpha, beta, m) * sin (theta);
  /* Set phi velocity */
  vel_phi = 0;

  /* Compute euclidean coordinates */
  x_out[0] =
    vel_r * sin (theta) * cos (phi) + vel_theta * cos (theta) * cos (phi)
    - vel_phi * sin (phi);
  x_out[1] =
    vel_r * sin (theta) * sin (phi) + vel_theta * cos (theta) * sin (phi) +
    vel_phi * cos (phi);
  x_out[2] = vel_r * cos (theta) - vel_theta * cos (theta);
}

T8_EXTERN_C_END ();
