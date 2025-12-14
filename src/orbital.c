/* orbital.c - Orbital Mechanics Functions
 * C89 portable for DOS, Linux
 * 
 * Essential calculations for lunar/Mars landing missions
 */

#include "config.h"

#ifdef HAVE_ORBITAL

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "orbital.h"
#include "apf.h"
#include "apfx.h"

/* External display digits */
extern int display_digits;

/* Standard gravity (m/s^2) */
#define G0 "9.80665"

/* ========== Helper Functions ========== */

static void get_pi(apf *pi) {
    apf_from_str(pi, "3.14159265358979323846264338327950288419");
}

static void get_two_pi(apf *two_pi) {
    apf pi;
    get_pi(&pi);
    apf_add(two_pi, &pi, &pi);
}

/* ========== Orbital Calculations ========== */

/* Circular orbital velocity: v = sqrt(mu/r) */
void orbital_velocity(apf *v, const apf *mu, const apf *r) {
    apf ratio;
    apf_div(&ratio, mu, r);
    apf_sqrt(v, &ratio);
}

/* Escape velocity: v_esc = sqrt(2*mu/r) */
void escape_velocity(apf *v, const apf *mu, const apf *r) {
    apf two, ratio, tmp;
    apf_from_int(&two, 2);
    apf_mul(&tmp, &two, mu);
    apf_div(&ratio, &tmp, r);
    apf_sqrt(v, &ratio);
}

/* Orbital period: T = 2*pi*sqrt(a^3/mu) */
void orbital_period(apf *T, const apf *mu, const apf *a) {
    apf two_pi, a3, ratio, sqrt_ratio;
    apf tmp;
    
    get_two_pi(&two_pi);
    
    /* a^3 */
    apf_mul(&tmp, a, a);
    apf_mul(&a3, &tmp, a);
    
    /* a^3 / mu */
    apf_div(&ratio, &a3, mu);
    
    /* sqrt(a^3/mu) */
    apf_sqrt(&sqrt_ratio, &ratio);
    
    /* 2*pi*sqrt(...) */
    apf_mul(T, &two_pi, &sqrt_ratio);
}

/* Semi-major axis from period: a = (mu * T^2 / (4*pi^2))^(1/3) */
void semimajor_from_period(apf *a, const apf *mu, const apf *T) {
    apf pi, four_pi_sq, T_sq, numer, ratio, one_third;
    
    get_pi(&pi);
    
    /* 4*pi^2 */
    apf_mul(&four_pi_sq, &pi, &pi);
    {
        apf four;
        apf_from_int(&four, 4);
        apf_mul(&four_pi_sq, &four_pi_sq, &four);
    }
    
    /* T^2 */
    apf_mul(&T_sq, T, T);
    
    /* mu * T^2 */
    apf_mul(&numer, mu, &T_sq);
    
    /* mu * T^2 / (4*pi^2) */
    apf_div(&ratio, &numer, &four_pi_sq);
    
    /* ^(1/3) */
    apf_from_str(&one_third, "0.333333333333333333333333333333333333");
    apfx_pow(a, &ratio, &one_third);
}

/* Vis-viva equation: v = sqrt(mu * (2/r - 1/a)) */
void vis_viva(apf *v, const apf *mu, const apf *r, const apf *a) {
    apf two, two_over_r, one_over_a, diff, product;
    
    apf_from_int(&two, 2);
    
    /* 2/r */
    apf_div(&two_over_r, &two, r);
    
    /* 1/a */
    {
        apf one;
        apf_from_int(&one, 1);
        apf_div(&one_over_a, &one, a);
    }
    
    /* 2/r - 1/a */
    apf_sub(&diff, &two_over_r, &one_over_a);
    
    /* mu * (2/r - 1/a) */
    apf_mul(&product, mu, &diff);
    
    /* sqrt(...) */
    apf_sqrt(v, &product);
}

/* Specific orbital energy: E = -mu / (2*a) */
void specific_energy(apf *E, const apf *mu, const apf *a) {
    apf two, two_a;
    
    apf_from_int(&two, 2);
    apf_mul(&two_a, &two, a);
    apf_div(E, mu, &two_a);
    apf_neg(E, E);
}

/* Hohmann transfer delta-V */
void hohmann_deltav(apf *dv, const apf *mu, const apf *r1, const apf *r2) {
    apf v1_circ, v2_circ, v1_trans, v2_trans;
    apf a_trans, r_sum;
    apf dv1, dv2;
    apf two;
    
    apf_from_int(&two, 2);
    
    /* Transfer orbit semi-major axis: a = (r1 + r2) / 2 */
    apf_add(&r_sum, r1, r2);
    apf_div(&a_trans, &r_sum, &two);
    
    /* Circular velocities */
    orbital_velocity(&v1_circ, mu, r1);
    orbital_velocity(&v2_circ, mu, r2);
    
    /* Transfer velocities at periapsis and apoapsis */
    vis_viva(&v1_trans, mu, r1, &a_trans);
    vis_viva(&v2_trans, mu, r2, &a_trans);
    
    /* Delta-Vs */
    apf_sub(&dv1, &v1_trans, &v1_circ);
    apf_abs(&dv1, &dv1);
    
    apf_sub(&dv2, &v2_circ, &v2_trans);
    apf_abs(&dv2, &dv2);
    
    /* Total delta-V */
    apf_add(dv, &dv1, &dv2);
}

/* Hohmann transfer time: t = pi * sqrt((r1+r2)^3 / (8*mu)) */
void hohmann_time(apf *t, const apf *mu, const apf *r1, const apf *r2) {
    apf pi, eight, r_sum, r_sum3, denom, ratio;
    apf tmp;
    
    get_pi(&pi);
    apf_from_int(&eight, 8);
    
    /* r1 + r2 */
    apf_add(&r_sum, r1, r2);
    
    /* (r1+r2)^3 */
    apf_mul(&tmp, &r_sum, &r_sum);
    apf_mul(&r_sum3, &tmp, &r_sum);
    
    /* 8*mu */
    apf_mul(&denom, &eight, mu);
    
    /* (r1+r2)^3 / (8*mu) */
    apf_div(&ratio, &r_sum3, &denom);
    
    /* sqrt(...) */
    apf_sqrt(&tmp, &ratio);
    
    /* pi * sqrt(...) */
    apf_mul(t, &pi, &tmp);
}

/* Sphere of influence: r_SOI = a * (m/M)^(2/5) */
void sphere_of_influence(apf *r_soi, const apf *a, const apf *m_body, const apf *m_central) {
    apf ratio, exponent, powered;
    
    /* m/M */
    apf_div(&ratio, m_body, m_central);
    
    /* 2/5 = 0.4 */
    apf_from_str(&exponent, "0.4");
    
    /* (m/M)^0.4 */
    apfx_pow(&powered, &ratio, &exponent);
    
    /* a * (m/M)^0.4 */
    apf_mul(r_soi, a, &powered);
}

/* Plane change delta-V: dv = 2 * v * sin(theta/2) */
void plane_change_deltav(apf *dv, const apf *v, const apf *theta) {
    apf two, half_theta, sin_half;
    
    apf_from_int(&two, 2);
    
    /* theta/2 */
    apf_div(&half_theta, theta, &two);
    
    /* sin(theta/2) */
    apfx_sin(&sin_half, &half_theta);
    
    /* 2 * v * sin(theta/2) */
    apf_mul(dv, &two, v);
    apf_mul(dv, dv, &sin_half);
}

/* Tsiolkovsky: dv = Isp * g0 * ln(m0/mf) */
void tsiolkovsky_deltav(apf *dv, const apf *isp, const apf *m0, const apf *mf) {
    apf g0, ve, ratio, ln_ratio;
    
    apf_from_str(&g0, G0);
    
    /* Exhaust velocity: ve = Isp * g0 */
    apf_mul(&ve, isp, &g0);
    
    /* Mass ratio: m0/mf */
    apf_div(&ratio, m0, mf);
    
    /* ln(m0/mf) */
    apfx_log(&ln_ratio, &ratio);
    
    /* dv = ve * ln(m0/mf) */
    apf_mul(dv, &ve, &ln_ratio);
}

/* Tsiolkovsky: mass ratio = exp(-dv / (Isp * g0)) */
void tsiolkovsky_mass_ratio(apf *ratio, const apf *isp, const apf *dv) {
    apf g0, ve, exponent;
    
    apf_from_str(&g0, G0);
    
    /* Exhaust velocity */
    apf_mul(&ve, isp, &g0);
    
    /* -dv / ve */
    apf_div(&exponent, dv, &ve);
    apf_neg(&exponent, &exponent);
    
    /* exp(-dv/ve) */
    apfx_exp(ratio, &exponent);
}

/* Gravity assist maximum delta-V */
void gravity_assist_deltav(apf *dv, const apf *v_inf, const apf *mu, const apf *r_periapsis) {
    /* Maximum deflection is 180 degrees (theoretical limit) */
    /* Actual deflection: delta = 2 * arcsin(1 / (1 + r_p * v_inf^2 / mu)) */
    /* For maximum delta-V contribution: dv = 2 * v_inf */
    apf two;
    apf_from_int(&two, 2);
    apf_mul(dv, &two, v_inf);
    
    /* Suppress unused parameter warnings */
    (void)mu;
    (void)r_periapsis;
}

/* ========== Kepler's Equation ========== */

/* Mean anomaly from eccentric anomaly: M = E - e * sin(E) */
void mean_anomaly(apf *M, const apf *E, const apf *e) {
    apf sin_E, e_sin_E;
    
    apfx_sin(&sin_E, E);
    apf_mul(&e_sin_E, e, &sin_E);
    apf_sub(M, E, &e_sin_E);
}

/* Eccentric anomaly from mean anomaly (Newton-Raphson) */
int eccentric_anomaly(apf *E, const apf *M, const apf *e) {
    apf E_new, f, df, delta, tolerance;
    apf sin_E, cos_E, one;
    int iter;
    
    apf_from_int(&one, 1);
    apf_from_str(&tolerance, "1e-15");
    
    /* Initial guess: E = M */
    apf_copy(E, M);
    
    for (iter = 0; iter < 50; iter++) {
        /* f(E) = E - e*sin(E) - M */
        apfx_sin(&sin_E, E);
        apf_mul(&f, e, &sin_E);
        apf_sub(&f, E, &f);
        apf_sub(&f, &f, M);
        
        /* f'(E) = 1 - e*cos(E) */
        apfx_cos(&cos_E, E);
        apf_mul(&df, e, &cos_E);
        apf_sub(&df, &one, &df);
        
        /* E_new = E - f/f' */
        apf_div(&delta, &f, &df);
        apf_sub(&E_new, E, &delta);
        
        /* Check convergence */
        apf_abs(&delta, &delta);
        if (apf_cmp(&delta, &tolerance) < 0) {
            apf_copy(E, &E_new);
            return 1;
        }
        
        apf_copy(E, &E_new);
    }
    
    return 0;  /* Did not converge */
}

/* True anomaly from eccentric anomaly */
void true_anomaly(apf *nu, const apf *E, const apf *e) {
    apf one, one_plus_e, one_minus_e, ratio, sqrt_ratio;
    apf half_E, tan_half_E, tan_half_nu, two;
    
    apf_from_int(&one, 1);
    apf_from_int(&two, 2);
    
    /* (1+e)/(1-e) */
    apf_add(&one_plus_e, &one, e);
    apf_sub(&one_minus_e, &one, e);
    apf_div(&ratio, &one_plus_e, &one_minus_e);
    apf_sqrt(&sqrt_ratio, &ratio);
    
    /* tan(E/2) */
    apf_div(&half_E, E, &two);
    {
        apf sin_half, cos_half;
        apfx_sin(&sin_half, &half_E);
        apfx_cos(&cos_half, &half_E);
        apf_div(&tan_half_E, &sin_half, &cos_half);
    }
    
    /* tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2) */
    apf_mul(&tan_half_nu, &sqrt_ratio, &tan_half_E);
    
    /* nu = 2 * atan(tan_half_nu) */
    apfx_atan(nu, &tan_half_nu);
    apf_mul(nu, nu, &two);
}

/* Radius from true anomaly: r = a*(1-e^2)/(1+e*cos(nu)) */
void radius_from_anomaly(apf *r, const apf *a, const apf *e, const apf *nu) {
    apf one, e_sq, one_minus_e_sq, cos_nu, e_cos_nu, denom, numer;
    
    apf_from_int(&one, 1);
    
    /* 1 - e^2 */
    apf_mul(&e_sq, e, e);
    apf_sub(&one_minus_e_sq, &one, &e_sq);
    
    /* a * (1 - e^2) */
    apf_mul(&numer, a, &one_minus_e_sq);
    
    /* e * cos(nu) */
    apfx_cos(&cos_nu, nu);
    apf_mul(&e_cos_nu, e, &cos_nu);
    
    /* 1 + e*cos(nu) */
    apf_add(&denom, &one, &e_cos_nu);
    
    /* r = numer / denom */
    apf_div(r, &numer, &denom);
}

/* ========== Constants Access ========== */

int get_mu(apf *mu, const char *body) {
    if (strcmp(body, "sun") == 0 || strcmp(body, "Sun") == 0) {
        apf_from_str(mu, MU_SUN);
        return 1;
    }
    if (strcmp(body, "earth") == 0 || strcmp(body, "Earth") == 0) {
        apf_from_str(mu, MU_EARTH);
        return 1;
    }
    if (strcmp(body, "moon") == 0 || strcmp(body, "Moon") == 0) {
        apf_from_str(mu, MU_MOON);
        return 1;
    }
    if (strcmp(body, "mars") == 0 || strcmp(body, "Mars") == 0) {
        apf_from_str(mu, MU_MARS);
        return 1;
    }
    if (strcmp(body, "jupiter") == 0 || strcmp(body, "Jupiter") == 0) {
        apf_from_str(mu, MU_JUPITER);
        return 1;
    }
    return 0;
}

int get_radius(apf *r, const char *body) {
    if (strcmp(body, "earth") == 0 || strcmp(body, "Earth") == 0) {
        apf_from_str(r, R_EARTH);
        return 1;
    }
    if (strcmp(body, "moon") == 0 || strcmp(body, "Moon") == 0) {
        apf_from_str(r, R_MOON);
        return 1;
    }
    if (strcmp(body, "mars") == 0 || strcmp(body, "Mars") == 0) {
        apf_from_str(r, R_MARS);
        return 1;
    }
    if (strcmp(body, "sun") == 0 || strcmp(body, "Sun") == 0) {
        apf_from_str(r, R_SUN);
        return 1;
    }
    return 0;
}

void cmd_orbital(const char *args) {
    char buf[64];
    const char *p = args;
    int digits = display_digits ? display_digits : 10;
    
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (*p == '\0') {
        printf("Orbital Mechanics Commands:\n");
        printf("  orb v_circ <mu> <r>     - Circular orbital velocity\n");
        printf("  orb v_esc <mu> <r>      - Escape velocity\n");
        printf("  orb period <mu> <a>     - Orbital period\n");
        printf("  orb hohmann <mu> <r1> <r2> - Hohmann transfer delta-V\n");
        printf("  orb soi <a> <m> <M>     - Sphere of influence\n");
        printf("  orb mu <body>           - Get gravitational parameter\n");
        printf("  orb radius <body>       - Get body radius\n");
        printf("\nBodies: sun, earth, moon, mars, jupiter\n");
        printf("Units: km and km/s (mu in km^3/s^2)\n");
        return;
    }
    
    if (strncmp(p, "v_circ ", 7) == 0) {
        apf mu, r, v;
        p += 7;
        apf_from_str(&mu, p);
        p = strchr(p, ' ');
        if (p) {
            apf_from_str(&r, p + 1);
            orbital_velocity(&v, &mu, &r);
            apf_to_str(buf, sizeof(buf), &v, digits);
            printf("v_circ = %s km/s\n", buf);
        }
        return;
    }
    
    if (strncmp(p, "v_esc ", 6) == 0) {
        apf mu, r, v;
        p += 6;
        apf_from_str(&mu, p);
        p = strchr(p, ' ');
        if (p) {
            apf_from_str(&r, p + 1);
            escape_velocity(&v, &mu, &r);
            apf_to_str(buf, sizeof(buf), &v, digits);
            printf("v_escape = %s km/s\n", buf);
        }
        return;
    }
    
    if (strncmp(p, "period ", 7) == 0) {
        apf mu, a, T;
        p += 7;
        apf_from_str(&mu, p);
        p = strchr(p, ' ');
        if (p) {
            apf_from_str(&a, p + 1);
            orbital_period(&T, &mu, &a);
            apf_to_str(buf, sizeof(buf), &T, digits);
            printf("Period = %s seconds\n", buf);
            /* Also show in hours/days */
            {
                apf hours, days;
                apf sixty, twentyfour;
                apf_from_int(&sixty, 60);
                apf_from_int(&twentyfour, 24);
                apf_div(&hours, &T, &sixty);
                apf_div(&hours, &hours, &sixty);
                apf_div(&days, &hours, &twentyfour);
                apf_to_str(buf, sizeof(buf), &hours, 4);
                printf("       = %s hours\n", buf);
                apf_to_str(buf, sizeof(buf), &days, 4);
                printf("       = %s days\n", buf);
            }
        }
        return;
    }
    
    if (strncmp(p, "hohmann ", 8) == 0) {
        apf mu, r1, r2, dv, t;
        p += 8;
        apf_from_str(&mu, p);
        p = strchr(p, ' ');
        if (p) {
            p++;
            apf_from_str(&r1, p);
            p = strchr(p, ' ');
            if (p) {
                apf_from_str(&r2, p + 1);
                hohmann_deltav(&dv, &mu, &r1, &r2);
                hohmann_time(&t, &mu, &r1, &r2);
                apf_to_str(buf, sizeof(buf), &dv, digits);
                printf("Hohmann delta-V = %s km/s\n", buf);
                {
                    apf hours, days;
                    apf sixty, twentyfour;
                    apf_from_int(&sixty, 60);
                    apf_from_int(&twentyfour, 24);
                    apf_div(&hours, &t, &sixty);
                    apf_div(&hours, &hours, &sixty);
                    apf_div(&days, &hours, &twentyfour);
                    apf_to_str(buf, sizeof(buf), &days, 4);
                    printf("Transfer time = %s days\n", buf);
                }
            }
        }
        return;
    }
    
    if (strncmp(p, "mu ", 3) == 0) {
        apf mu;
        p += 3;
        if (get_mu(&mu, p)) {
            apf_to_str(buf, sizeof(buf), &mu, digits);
            printf("mu(%s) = %s km^3/s^2\n", p, buf);
        } else {
            printf("Unknown body: %s\n", p);
        }
        return;
    }
    
    if (strncmp(p, "radius ", 7) == 0) {
        apf r;
        p += 7;
        if (get_radius(&r, p)) {
            apf_to_str(buf, sizeof(buf), &r, digits);
            printf("R(%s) = %s km\n", p, buf);
        } else {
            printf("Unknown body: %s\n", p);
        }
        return;
    }
    
    printf("Unknown orbital command. Type 'orb' for help.\n");
}

#endif /* HAVE_ORBITAL */
