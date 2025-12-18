/* datetime.h - DateTime and timetable operations for scalc
 * C89 compliant for Watcom C / DOS
 */
#ifndef DATETIME_H
#define DATETIME_H

#include "matrix.h"
#include "apf.h"

/* DateTime parsing and formatting */
int datetime_parse(apf *result, const char *str);
void datetime_format(char *buf, int bufsize, const apf *seconds);

/* Duration conversions (to seconds) */
void duration_seconds(apf *result, const apf *val);
void duration_minutes(apf *result, const apf *val);
void duration_hours(apf *result, const apf *val);
void duration_days(apf *result, const apf *val);
void duration_milliseconds(apf *result, const apf *val);

/* Interpolation */
void interp_linear_col(matrix_t *result, int col,
                       const matrix_t *times, const matrix_t *data,
                       const matrix_t *new_times);
void interp_spline_col(matrix_t *result, int col,
                       const matrix_t *times, const matrix_t *data,
                       const matrix_t *new_times);

/* Timetable operations */
int retime_generate_times(matrix_t *result, const matrix_t *times, const apf *step);
int timetable_retime(matrix_t *new_times_out, matrix_t *new_data_out,
                     const matrix_t *times, const matrix_t *data,
                     const char *time_spec, const char *method);

#endif /* DATETIME_H */
