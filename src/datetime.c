/* datetime.c - DateTime and timetable operations for scalc
 * C89 compliant for Watcom C / DOS
 * NO double/float - all math via apf for FPU-less systems
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "apfx.h"

/* ========== DateTime Parsing ========== */

/* Parse datetime string to Unix timestamp (seconds since 1970-01-01 00:00:00 UTC) 
 * Formats supported:
 *   'YYYY-MM-DD HH:MM:SS'
 *   'YYYY-MM-DD HH:MM'
 *   'YYYY-MM-DD'
 * Result stored in apf
 * Returns 1 on success, 0 on error
 */
int datetime_parse(apf *result, const char *str)
{
    int year, month, day, hour = 0, min = 0, sec = 0;
    long days;
    int i;
    apf tmp, tmp2;
    
    /* Days in each month (non-leap year) */
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    
    if (!str) return 0;
    
    /* Try full datetime format: YYYY-MM-DD HH:MM:SS */
    if (sscanf(str, "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &min, &sec) >= 3) {
        /* Valid */
    } else {
        return 0;
    }
    
    /* Validate */
    if (year < 1900 || year > 2100) return 0;
    if (month < 1 || month > 12) return 0;
    if (day < 1 || day > 31) return 0;
    if (hour < 0 || hour > 23) return 0;
    if (min < 0 || min > 59) return 0;
    if (sec < 0 || sec > 59) return 0;
    
    /* Calculate days since 1970-01-01 (Unix epoch) */
    days = 0;
    
    /* Years */
    for (i = 1970; i < year; i++) {
        days += 365;
        if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) {
            days += 1;  /* Leap year */
        }
    }
    for (i = year; i < 1970; i++) {
        days -= 365;
        if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) {
            days -= 1;  /* Leap year */
        }
    }
    
    /* Months */
    for (i = 1; i < month; i++) {
        days += days_in_month[i-1];
        if (i == 2 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) {
            days += 1;  /* Feb in leap year */
        }
    }
    
    /* Days */
    days += (day - 1);
    
    /* Convert to seconds: days * 86400 + hour * 3600 + min * 60 + sec */
    /* All arithmetic via apf */
    apf_from_int(result, days);
    apf_from_int(&tmp, 86400);
    apf_mul(&tmp2, result, &tmp);
    
    apf_from_int(&tmp, hour);
    apf_from_int(result, 3600);
    apf_mul(&tmp, &tmp, result);
    apf_add(&tmp2, &tmp2, &tmp);
    
    apf_from_int(&tmp, min);
    apf_from_int(result, 60);
    apf_mul(&tmp, &tmp, result);
    apf_add(&tmp2, &tmp2, &tmp);
    
    apf_from_int(&tmp, sec);
    apf_add(result, &tmp2, &tmp);
    
    return 1;
}

/* Format seconds since epoch to datetime string 
 * Uses apf integer operations only
 */
void datetime_format(char *buf, int bufsize, const apf *seconds)
{
    long days_long, rem;
    int year, month, day, hour, min, sec;
    int i;
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int leap;
    apf days_apf, secs_per_day, time_of_day;
    apf tmp, tmp2;
    
    /* Divide by 86400 to get days */
    apf_from_int(&secs_per_day, 86400);
    apf_div(&days_apf, seconds, &secs_per_day);
    apf_floor(&tmp, &days_apf);
    
    /* Get time of day = seconds - days * 86400 */
    apf_mul(&tmp2, &tmp, &secs_per_day);
    apf_sub(&time_of_day, seconds, &tmp2);
    
    /* Convert to long integers using apf_to_str and parsing */
    {
        char sbuf[64];
        apf_to_str(sbuf, sizeof(sbuf), &tmp, 0);
        days_long = atol(sbuf);
        apf_to_str(sbuf, sizeof(sbuf), &time_of_day, 0);
        rem = atol(sbuf);
    }
    
    if (rem < 0) { rem += 86400; days_long--; }
    
    hour = (int)(rem / 3600);
    rem %= 3600;
    min = (int)(rem / 60);
    sec = (int)(rem % 60);
    
    /* Convert days since 1970-01-01 (Unix epoch) to year/month/day */
    year = 1970;
    
    if (days_long >= 0) {
        while (1) {
            int year_days;
            leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
            year_days = 365 + leap;
            if (days_long < year_days) break;
            days_long -= year_days;
            year++;
        }
    } else {
        while (days_long < 0) {
            year--;
            leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
            days_long += 365 + leap;
        }
    }
    
    leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
    
    /* Month */
    month = 1;
    for (i = 0; i < 12; i++) {
        int dm = days_in_month[i];
        if (i == 1 && leap) dm++;
        if (days_long < dm) break;
        days_long -= dm;
        month++;
    }
    day = (int)days_long + 1;
    
    sprintf(buf, "%04d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, min, sec);
    (void)bufsize;  /* Suppress unused warning in C89 mode */
}

/* ========== Duration Functions ========== */
/* Convert values to seconds using apf arithmetic */

void duration_seconds(apf *result, const apf *val)
{
    apf_copy(result, val);
}

void duration_minutes(apf *result, const apf *val)
{
    apf sixty;
    apf_from_int(&sixty, 60);
    apf_mul(result, val, &sixty);
}

void duration_hours(apf *result, const apf *val)
{
    apf factor;
    apf_from_int(&factor, 3600);
    apf_mul(result, val, &factor);
}

void duration_days(apf *result, const apf *val)
{
    apf factor;
    apf_from_int(&factor, 86400);
    apf_mul(result, val, &factor);
}

void duration_milliseconds(apf *result, const apf *val)
{
    apf factor;
    apf_from_int(&factor, 1000);
    apf_div(result, val, &factor);
}

/* ========== Interpolation (APF only) ========== */

/* Linear interpolation between two points 
 * result = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
 */
void interp_linear_point(apfc *result, 
                         const apf *x1, const apfc *y1,
                         const apf *x2, const apfc *y2,
                         const apf *x)
{
    apf dx, t;
    apfc dy, scaled;
    
    /* dx = x2 - x1 */
    apf_sub(&dx, x2, x1);
    
    if (apf_iszero(&dx)) {
        /* Points at same x, return y1 */
        apf_copy(&result->re, &y1->re);
        apf_copy(&result->im, &y1->im);
        return;
    }
    
    /* t = (x - x1) / dx */
    apf_sub(&t, x, x1);
    apf_div(&t, &t, &dx);
    
    /* dy = y2 - y1 */
    apf_sub(&dy.re, &y2->re, &y1->re);
    apf_sub(&dy.im, &y2->im, &y1->im);
    
    /* result = y1 + dy * t */
    apf_mul(&scaled.re, &dy.re, &t);
    apf_mul(&scaled.im, &dy.im, &t);
    apf_add(&result->re, &y1->re, &scaled.re);
    apf_add(&result->im, &y1->im, &scaled.im);
}

/* Binary search for interval containing x in sorted array */
int find_interval(const matrix_t *times, const apf *x)
{
    int lo = 0, hi = times->rows - 1, mid;
    
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (apf_gt(&MAT_AT(times, mid, 0).re, x)) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return lo;
}

/* Linear interpolation for a column of data at new times */
void interp_linear_col(matrix_t *result, int col,
                       const matrix_t *times, const matrix_t *data,
                       const matrix_t *new_times)
{
    int i, idx;
    int n_old = times->rows;
    int n_new = new_times->rows;
    
    for (i = 0; i < n_new; i++) {
        apf x = MAT_AT(new_times, i, 0).re;
        
        /* Check bounds */
        if (apf_le(&x, &MAT_AT(times, 0, 0).re)) {
            /* Before first point - extrapolate from first two */
            if (n_old >= 2) {
                interp_linear_point(&MAT_AT(result, i, col),
                    &MAT_AT(times, 0, 0).re, &MAT_AT(data, 0, col),
                    &MAT_AT(times, 1, 0).re, &MAT_AT(data, 1, col),
                    &x);
            } else {
                MAT_AT(result, i, col) = MAT_AT(data, 0, col);
            }
        } else if (apf_ge(&x, &MAT_AT(times, n_old-1, 0).re)) {
            /* After last point - extrapolate from last two */
            if (n_old >= 2) {
                interp_linear_point(&MAT_AT(result, i, col),
                    &MAT_AT(times, n_old-2, 0).re, &MAT_AT(data, n_old-2, col),
                    &MAT_AT(times, n_old-1, 0).re, &MAT_AT(data, n_old-1, col),
                    &x);
            } else {
                MAT_AT(result, i, col) = MAT_AT(data, n_old-1, col);
            }
        } else {
            /* Interior - find interval and interpolate */
            idx = find_interval(times, &x);
            interp_linear_point(&MAT_AT(result, i, col),
                &MAT_AT(times, idx, 0).re, &MAT_AT(data, idx, col),
                &MAT_AT(times, idx+1, 0).re, &MAT_AT(data, idx+1, col),
                &x);
        }
    }
}

/* ========== Spline Interpolation (APF only) ========== */

/* Solve tridiagonal system for spline coefficients
 * Using Thomas algorithm - all in apf
 */
int spline_solve_tridiag(apf *y2, const matrix_t *times, const matrix_t *data, int col)
{
    int n = times->rows;
    int i;
    apf *c_prime;  /* Modified coefficients */
    apf *d_prime;  /* Modified RHS */
    apf tmp, tmp2;
    
    if (n < 2) return 0;
    
    c_prime = (apf *)mat_arena_alloc(n * sizeof(apf));
    d_prime = (apf *)mat_arena_alloc(n * sizeof(apf));
    if (!c_prime || !d_prime) return 0;
    
    /* Natural spline: y2[0] = y2[n-1] = 0 */
    apf_zero(&y2[0]);
    apf_zero(&y2[n-1]);
    
    if (n == 2) return 1;  /* Linear case */
    
    /* Forward sweep */
    /* For natural spline, first row is: y2[0] = 0 */
    apf_zero(&c_prime[0]);
    apf_zero(&d_prime[0]);
    
    for (i = 1; i < n - 1; i++) {
        apf h_prev, h_curr;
        apf a_i, b_i, c_i, d_i;  /* Tridiagonal coefficients */
        apf denom;
        apf slope_prev, slope_curr;
        
        /* h[i-1] = x[i] - x[i-1] */
        apf_sub(&h_prev, &MAT_AT(times, i, 0).re, &MAT_AT(times, i-1, 0).re);
        /* h[i] = x[i+1] - x[i] */
        apf_sub(&h_curr, &MAT_AT(times, i+1, 0).re, &MAT_AT(times, i, 0).re);
        
        /* a[i] = h[i-1] */
        apf_copy(&a_i, &h_prev);
        
        /* b[i] = 2 * (h[i-1] + h[i]) */
        apf_add(&b_i, &h_prev, &h_curr);
        apf_from_int(&tmp, 2);
        apf_mul(&b_i, &b_i, &tmp);
        
        /* c[i] = h[i] */
        apf_copy(&c_i, &h_curr);
        
        /* d[i] = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1]) */
        apf_sub(&slope_curr, &MAT_AT(data, i+1, col).re, &MAT_AT(data, i, col).re);
        apf_div(&slope_curr, &slope_curr, &h_curr);
        apf_sub(&slope_prev, &MAT_AT(data, i, col).re, &MAT_AT(data, i-1, col).re);
        apf_div(&slope_prev, &slope_prev, &h_prev);
        apf_sub(&d_i, &slope_curr, &slope_prev);
        apf_from_int(&tmp, 6);
        apf_mul(&d_i, &d_i, &tmp);
        
        /* Modified coefficients for Thomas algorithm */
        /* c'[i] = c[i] / (b[i] - a[i] * c'[i-1]) */
        apf_mul(&tmp, &a_i, &c_prime[i-1]);
        apf_sub(&denom, &b_i, &tmp);
        apf_div(&c_prime[i], &c_i, &denom);
        
        /* d'[i] = (d[i] - a[i] * d'[i-1]) / (b[i] - a[i] * c'[i-1]) */
        apf_mul(&tmp, &a_i, &d_prime[i-1]);
        apf_sub(&tmp2, &d_i, &tmp);
        apf_div(&d_prime[i], &tmp2, &denom);
    }
    
    /* Back substitution */
    /* y2[n-1] = 0 (natural spline) */
    for (i = n - 2; i >= 1; i--) {
        /* y2[i] = d'[i] - c'[i] * y2[i+1] */
        apf_mul(&tmp, &c_prime[i], &y2[i+1]);
        apf_sub(&y2[i], &d_prime[i], &tmp);
    }
    
    return 1;
}

/* Evaluate cubic spline at point x */
void spline_eval_point(apfc *result, 
                       const matrix_t *times, const matrix_t *data,
                       const apf *y2, int col, const apf *x)
{
    int n = times->rows;
    int lo;
    apf h, a, b;
    apf tmp, tmp2, tmp3;
    apf six;
    
    /* Find interval */
    lo = find_interval(times, x);
    if (lo >= n - 1) lo = n - 2;
    if (lo < 0) lo = 0;
    
    /* h = x[hi] - x[lo] */
    apf_sub(&h, &MAT_AT(times, lo+1, 0).re, &MAT_AT(times, lo, 0).re);
    
    if (apf_iszero(&h)) {
        apf_copy(&result->re, &MAT_AT(data, lo, col).re);
        apf_zero(&result->im);
        return;
    }
    
    /* a = (x[hi] - x) / h */
    apf_sub(&a, &MAT_AT(times, lo+1, 0).re, x);
    apf_div(&a, &a, &h);
    
    /* b = (x - x[lo]) / h */
    apf_sub(&b, x, &MAT_AT(times, lo, 0).re);
    apf_div(&b, &b, &h);
    
    /* result = a*y[lo] + b*y[hi] + ((a^3-a)*y2[lo] + (b^3-b)*y2[hi]) * h^2/6 */
    apf_from_int(&six, 6);
    
    /* a * y[lo] */
    apf_mul(&tmp, &a, &MAT_AT(data, lo, col).re);
    
    /* b * y[hi] */
    apf_mul(&tmp2, &b, &MAT_AT(data, lo+1, col).re);
    apf_add(&tmp, &tmp, &tmp2);
    
    /* (a^3 - a) */
    apf_mul(&tmp2, &a, &a);
    apf_mul(&tmp2, &tmp2, &a);
    apf_sub(&tmp2, &tmp2, &a);
    
    /* (a^3 - a) * y2[lo] */
    apf_mul(&tmp2, &tmp2, &y2[lo]);
    
    /* (b^3 - b) */
    apf_mul(&tmp3, &b, &b);
    apf_mul(&tmp3, &tmp3, &b);
    apf_sub(&tmp3, &tmp3, &b);
    
    /* (b^3 - b) * y2[hi] */
    apf_mul(&tmp3, &tmp3, &y2[lo+1]);
    
    /* Sum */
    apf_add(&tmp2, &tmp2, &tmp3);
    
    /* * h^2 / 6 */
    apf_mul(&tmp3, &h, &h);
    apf_div(&tmp3, &tmp3, &six);
    apf_mul(&tmp2, &tmp2, &tmp3);
    
    /* Final result */
    apf_add(&result->re, &tmp, &tmp2);
    apf_zero(&result->im);
}

/* Spline interpolation for a column */
void interp_spline_col(matrix_t *result, int col,
                       const matrix_t *times, const matrix_t *data,
                       const matrix_t *new_times)
{
    int n_old = times->rows;
    int n_new = new_times->rows;
    int i;
    apf *y2;
    
    /* Allocate second derivatives */
    y2 = (apf *)mat_arena_alloc(n_old * sizeof(apf));
    if (!y2) return;
    
    /* Compute spline coefficients */
    if (!spline_solve_tridiag(y2, times, data, col)) {
        /* Fall back to linear */
        interp_linear_col(result, col, times, data, new_times);
        return;
    }
    
    /* Evaluate at each new time */
    for (i = 0; i < n_new; i++) {
        spline_eval_point(&MAT_AT(result, i, col),
                         times, data, y2, col,
                         &MAT_AT(new_times, i, 0).re);
    }
}

/* ========== Aggregation Functions ========== */

/* Mean of values in a bin (start <= t < end) */
void aggregate_mean(apfc *result, const matrix_t *times, const matrix_t *data,
                    int col, const apf *start, const apf *end)
{
    int i, count = 0;
    apf sum;
    apf_zero(&sum);
    
    for (i = 0; i < times->rows; i++) {
        if (apf_ge(&MAT_AT(times, i, 0).re, start) &&
            apf_lt(&MAT_AT(times, i, 0).re, end)) {
            apf tmp;
            apf_add(&tmp, &sum, &MAT_AT(data, i, col).re);
            apf_copy(&sum, &tmp);
            count++;
        }
    }
    
    if (count > 0) {
        apf cnt;
        apf_from_int(&cnt, count);
        apf_div(&result->re, &sum, &cnt);
    } else {
        apf_zero(&result->re);  /* Or could set NaN */
    }
    apf_zero(&result->im);
}

/* ========== Retime Main Function ========== */

/* Generate regular time grid */
int retime_generate_times(matrix_t *result, const matrix_t *times, const apf *step)
{
    int n = times->rows;
    int i, count;
    apf tmin, tmax, t;
    apf start, end;
    apf tmp;
    
    /* Find min and max times */
    apf_copy(&tmin, &MAT_AT(times, 0, 0).re);
    apf_copy(&tmax, &MAT_AT(times, 0, 0).re);
    for (i = 1; i < n; i++) {
        if (apf_lt(&MAT_AT(times, i, 0).re, &tmin)) {
            apf_copy(&tmin, &MAT_AT(times, i, 0).re);
        }
        if (apf_gt(&MAT_AT(times, i, 0).re, &tmax)) {
            apf_copy(&tmax, &MAT_AT(times, i, 0).re);
        }
    }
    
    /* Round start down to step boundary */
    apf_div(&tmp, &tmin, step);
    apf_floor(&tmp, &tmp);
    apf_mul(&start, &tmp, step);
    
    /* Round end up to step boundary */
    apf_div(&tmp, &tmax, step);
    apf_ceil(&tmp, &tmp);
    apf_mul(&end, &tmp, step);
    
    /* Count points */
    count = 0;
    apf_copy(&t, &start);
    while (apf_le(&t, &end)) {
        count++;
        apf_add(&tmp, &t, step);
        apf_copy(&t, &tmp);
    }
    
    /* Allocate */
    mat_zero(result, count, 1);
    if (!result->data) return 0;
    
    /* Fill times */
    apf_copy(&t, &start);
    for (i = 0; i < count; i++) {
        apf_copy(&MAT_AT(result, i, 0).re, &t);
        apf_zero(&MAT_AT(result, i, 0).im);
        apf_add(&tmp, &t, step);
        apf_copy(&t, &tmp);
    }
    
    return 1;
}

/* Main retime function 
 * method: "linear", "spline", "mean"
 */
int timetable_retime(matrix_t *new_times_out, matrix_t *new_data_out,
                     const matrix_t *times, const matrix_t *data,
                     const char *time_spec, const char *method)
{
    int j;
    matrix_t new_times = {0, 0, NULL};
    apf step;
    
    /* Parse time specification */
    if (strcmp(time_spec, "hourly") == 0) {
        apf_from_int(&step, 3600);
    } else if (strcmp(time_spec, "minutely") == 0) {
        apf_from_int(&step, 60);
    } else if (strcmp(time_spec, "secondly") == 0) {
        apf_from_int(&step, 1);
    } else if (strcmp(time_spec, "daily") == 0) {
        apf_from_int(&step, 86400);
    } else {
        /* Assume numeric string for custom step in seconds */
        /* Check if it looks like a number */
        if (time_spec[0] == '\0' || 
            (time_spec[0] != '-' && time_spec[0] != '+' && 
             time_spec[0] != '.' && (time_spec[0] < '0' || time_spec[0] > '9'))) {
            printf("Error: Unknown time specification '%s'\n", time_spec);
            return 0;
        }
        apf_from_str(&step, time_spec);
    }
    
    /* Generate new time grid */
    if (!retime_generate_times(&new_times, times, &step)) {
        return 0;
    }
    
    /* Allocate output data */
    mat_zero(new_data_out, new_times.rows, data->cols);
    if (!new_data_out->data) return 0;
    
    /* Interpolate each column */
    for (j = 0; j < data->cols; j++) {
        if (strcmp(method, "linear") == 0) {
            interp_linear_col(new_data_out, j, times, data, &new_times);
        } else if (strcmp(method, "spline") == 0) {
            interp_spline_col(new_data_out, j, times, data, &new_times);
        } else if (strcmp(method, "mean") == 0) {
            /* Aggregate mean in each bin */
            int i;
            for (i = 0; i < new_times.rows; i++) {
                apf bin_start, bin_end, tmp;
                apf_copy(&bin_start, &MAT_AT(&new_times, i, 0).re);
                apf_add(&tmp, &bin_start, &step);
                apf_copy(&bin_end, &tmp);
                aggregate_mean(&MAT_AT(new_data_out, i, j), times, data, j,
                             &bin_start, &bin_end);
            }
        } else {
            printf("Error: Unknown interpolation method '%s'\n", method);
            return 0;
        }
    }
    
    /* Copy new times to output */
    mat_copy(new_times_out, &new_times);
    
    return 1;
}

/* ============================================================
 * Unix timestamp conversion functions for SaaS metrics
 * ============================================================ */

void unix_to_datetime(long ts, int *year, int *month, int *day, 
                      int *hour, int *min, int *sec)
{
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long days;
    int i, leap, dm, year_days;
    
    *sec = (int)(ts % 60);
    ts /= 60;
    *min = (int)(ts % 60);
    ts /= 60;
    *hour = (int)(ts % 24);
    days = ts / 24;
    
    /* Start from 1970 */
    *year = 1970;
    while (1) {
        leap = (*year % 4 == 0 && (*year % 100 != 0 || *year % 400 == 0));
        year_days = leap ? 366 : 365;
        if (days < year_days) break;
        days -= year_days;
        (*year)++;
    }
    
    leap = (*year % 4 == 0 && (*year % 100 != 0 || *year % 400 == 0));
    *month = 1;
    for (i = 0; i < 12; i++) {
        dm = days_in_month[i];
        if (i == 1 && leap) dm = 29;
        if (days < dm) break;
        days -= dm;
        (*month)++;
    }
    
    *day = (int)days + 1;
}

long datetime_to_unix(int year, int month, int day, int hour, int min, int sec)
{
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long days = 0;
    int y, i, leap;
    
    /* Days from 1970 to year */
    for (y = 1970; y < year; y++) {
        leap = (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0));
        days += leap ? 366 : 365;
    }
    
    /* Days in months */
    leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    for (i = 0; i < month - 1; i++) {
        days += days_in_month[i];
        if (i == 1 && leap) days++;
    }
    
    days += day - 1;
    
    return days * 86400 + hour * 3600 + min * 60 + sec;
}
