/* orbital.h - Orbital Mechanics Functions
 * C89 portable for DOS, Linux
 * 
 * Essential calculations for lunar/Mars landing missions:
 *   - Kepler's laws
 *   - Orbital velocity and period
 *   - Hohmann transfer
 *   - Delta-V calculations
 *   - Escape velocity
 *   - Sphere of influence
 *   - Gravitational parameters
 */

#ifndef ORBITAL_H
#define ORBITAL_H

#include "config.h"

#ifdef HAVE_ORBITAL

#include "apf.h"

/* ========== Physical Constants ========== */

/* Gravitational constant G (m^3 kg^-1 s^-2) */
#define G_CONST "6.67430e-11"

/* Standard gravitational parameters (mu = GM, km^3/s^2) */
#define MU_SUN     "1.32712440018e11"   /* Sun */
#define MU_EARTH   "3.986004418e5"      /* Earth */
#define MU_MOON    "4.9048695e3"        /* Moon */
#define MU_MARS    "4.282837e4"         /* Mars */
#define MU_JUPITER "1.26686534e8"       /* Jupiter */

/* Radii (km) */
#define R_EARTH    "6371.0"
#define R_MOON     "1737.4"
#define R_MARS     "3389.5"
#define R_SUN      "696340.0"

/* Mean orbital distances (km) */
#define DIST_EARTH_MOON   "384400.0"
#define DIST_EARTH_SUN    "1.496e8"
#define DIST_MARS_SUN     "2.279e8"

/* ========== Orbital Calculations ========== */

/* Orbital velocity at radius r around body with mu */
/* v = sqrt(mu/r) for circular orbit */
void orbital_velocity(apf *v, const apf *mu, const apf *r);

/* Escape velocity at radius r */
/* v_esc = sqrt(2*mu/r) */
void escape_velocity(apf *v, const apf *mu, const apf *r);

/* Orbital period (Kepler's 3rd law) */
/* T = 2*pi*sqrt(a^3/mu) */
void orbital_period(apf *T, const apf *mu, const apf *a);

/* Semi-major axis from period */
/* a = (mu * T^2 / (4*pi^2))^(1/3) */
void semimajor_from_period(apf *a, const apf *mu, const apf *T);

/* Vis-viva equation: velocity at r for orbit with semi-major axis a */
/* v = sqrt(mu * (2/r - 1/a)) */
void vis_viva(apf *v, const apf *mu, const apf *r, const apf *a);

/* Specific orbital energy */
/* E = -mu / (2*a) */
void specific_energy(apf *E, const apf *mu, const apf *a);

/* Hohmann transfer delta-V from circular orbit r1 to r2 */
/* Returns total delta-V (sum of both burns) */
void hohmann_deltav(apf *dv, const apf *mu, const apf *r1, const apf *r2);

/* Hohmann transfer time */
/* t = pi * sqrt((r1+r2)^3 / (8*mu)) */
void hohmann_time(apf *t, const apf *mu, const apf *r1, const apf *r2);

/* Sphere of influence radius */
/* r_SOI = a * (m/M)^(2/5) */
void sphere_of_influence(apf *r_soi, const apf *a, const apf *m_body, const apf *m_central);

/* Delta-V for plane change at velocity v, angle theta (radians) */
/* dv = 2 * v * sin(theta/2) */
void plane_change_deltav(apf *dv, const apf *v, const apf *theta);

/* Tsiolkovsky rocket equation */
/* dv = Isp * g0 * ln(m0/mf) */
/* or: mf/m0 = exp(-dv / (Isp * g0)) */
void tsiolkovsky_deltav(apf *dv, const apf *isp, const apf *m0, const apf *mf);
void tsiolkovsky_mass_ratio(apf *ratio, const apf *isp, const apf *dv);

/* Gravity assist delta-V (hyperbolic flyby) */
/* Maximum deflection delta-V */
void gravity_assist_deltav(apf *dv, const apf *v_inf, const apf *mu, const apf *r_periapsis);

/* ========== Kepler's Equation ========== */

/* Mean anomaly from eccentric anomaly */
/* M = E - e * sin(E) */
void mean_anomaly(apf *M, const apf *E, const apf *e);

/* Eccentric anomaly from mean anomaly (Newton-Raphson iteration) */
/* Solves M = E - e * sin(E) for E */
int eccentric_anomaly(apf *E, const apf *M, const apf *e);

/* True anomaly from eccentric anomaly */
/* tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2) */
void true_anomaly(apf *nu, const apf *E, const apf *e);

/* Radius from true anomaly (orbit equation) */
/* r = a * (1 - e^2) / (1 + e * cos(nu)) */
void radius_from_anomaly(apf *r, const apf *a, const apf *e, const apf *nu);

/* ========== Constants Access ========== */

/* Get gravitational parameter by body name */
int get_mu(apf *mu, const char *body);

/* Get radius by body name */
int get_radius(apf *r, const char *body);

/* Process orbital command */
void cmd_orbital(const char *args);

#endif /* HAVE_ORBITAL */
#endif /* ORBITAL_H */
