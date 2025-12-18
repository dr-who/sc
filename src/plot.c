/* plot.c - Function plotting for scalc
 * C89 compliant for Watcom C / DOS
 * 100% soft-float APF-based - NO double or float types
 */
#define _XOPEN_SOURCE 500  /* For usleep */
#define _POSIX_C_SOURCE 200809L
#include "sc.h"
#include <time.h>

#ifdef HAVE_CONIO
/* DOS VGA Graphics Mode 13h: 320x200, 256 colors */

#define GFX_WIDTH  320
#define GFX_HEIGHT 200
#define COLOR_YELLOW 14
#define COLOR_CYAN   3

static void set_video_mode(unsigned char mode)
{
    union REGS regs;
    regs.h.ah = 0x00;
    regs.h.al = mode;
    int86(0x10, &regs, &regs);
}

static void bios_put_pixel(int x, int y, unsigned char color)
{
    union REGS regs;
    if (x < 0 || x >= GFX_WIDTH || y < 0 || y >= GFX_HEIGHT) return;
    regs.h.ah = 0x0C;
    regs.h.al = color;
    regs.h.bh = 0;
    regs.x.cx = x;
    regs.x.dx = y;
    int86(0x10, &regs, &regs);
}

static void draw_hline(int x1, int x2, int y, unsigned char color)
{
    int x;
    if (y < 0 || y >= GFX_HEIGHT) return;
    if (x1 > x2) { int t = x1; x1 = x2; x2 = t; }
    if (x1 < 0) x1 = 0;
    if (x2 >= GFX_WIDTH) x2 = GFX_WIDTH - 1;
    for (x = x1; x <= x2; x++) {
        bios_put_pixel(x, y, color);
    }
}

static void draw_vline(int x, int y1, int y2, unsigned char color)
{
    int y;
    if (x < 0 || x >= GFX_WIDTH) return;
    if (y1 > y2) { int t = y1; y1 = y2; y2 = t; }
    if (y1 < 0) y1 = 0;
    if (y2 >= GFX_HEIGHT) y2 = GFX_HEIGHT - 1;
    for (y = y1; y <= y2; y++) {
        bios_put_pixel(x, y, color);
    }
}

/* Draw line between two points */
static void draw_line(int x1, int y1, int x2, int y2, unsigned char color)
{
    int dx, dy, steps, s;
    dx = x2 - x1;
    dy = y2 - y1;
    if (dx < 0) dx = -dx;
    if (dy < 0) dy = -dy;
    steps = dx > dy ? dx : dy;
    if (steps == 0) {
        bios_put_pixel(x1, y1, color);
        return;
    }
    for (s = 0; s <= steps; s++) {
        int px = x1 + (x2 - x1) * s / steps;
        int py = y1 + (y2 - y1) * s / steps;
        bios_put_pixel(px, py, color);
    }
}

int do_plot(const char *expr, char var, apf *xmin, apf *xmax)
{
    int func_idx, x;
    apf x_range, x_step, temp;
    apf ymin, ymax, y_range;
    apfc arg, result;
    int axis_x, axis_y;
    int prev_x, prev_y, first_point;
    int use_expr = 0;
    
    func_idx = get_func_index(expr);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        use_expr = 1;
    }
    
    apf_sub(&x_range, xmax, xmin);
    apf_from_int(&temp, GFX_WIDTH - 1);
    apf_div(&x_step, &x_range, &temp);
    
    /* Quick y-range scan */
    {
        int first_valid = 0;
        int eval_ok;
        for (x = 0; x < GFX_WIDTH; x += 40) {
            apf x_val, idx_apf;
            apf_from_int(&idx_apf, x);
            apf_mul(&temp, &idx_apf, &x_step);
            apf_add(&x_val, xmin, &temp);
            
            apf_copy(&arg.re, &x_val);
            apf_zero(&arg.im);
            
            if (use_expr) {
                eval_ok = eval_expr_with_var(&result, expr, var, &arg);
            } else {
                eval_ok = eval_user_func(&result, func_idx, &arg);
            }
            
            if (eval_ok && !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
                apf_iszero(&result.im)) {
                if (!first_valid) {
                    apf_copy(&ymin, &result.re);
                    apf_copy(&ymax, &result.re);
                    first_valid = 1;
                } else {
                    if (apf_cmp(&result.re, &ymin) < 0) apf_copy(&ymin, &result.re);
                    if (apf_cmp(&result.re, &ymax) > 0) apf_copy(&ymax, &result.re);
                }
            }
        }
        if (!first_valid) {
            printf("No valid points in range\n");
            return 0;
        }
    }
    
    /* Add 10% margin */
    apf_sub(&y_range, &ymax, &ymin);
    if (apf_iszero(&y_range)) apf_from_int(&y_range, 1);
    {
        apf margin, ten;
        apf_from_int(&ten, 10);
        apf_div(&margin, &y_range, &ten);
        apf_sub(&ymin, &ymin, &margin);
        apf_add(&ymax, &ymax, &margin);
        apf_sub(&y_range, &ymax, &ymin);
    }
    
    /* Calculate axis positions */
    axis_x = -1;
    axis_y = -1;
    {
        apf zero_apf, norm;
        apf_zero(&zero_apf);
        
        if (apf_cmp(&ymin, &zero_apf) <= 0 && apf_cmp(&ymax, &zero_apf) >= 0) {
            apf_sub(&temp, &zero_apf, &ymin);
            apf_div(&norm, &temp, &y_range);
            apf_from_int(&temp, GFX_HEIGHT - 1);
            apf_mul(&norm, &norm, &temp);
            axis_y = (GFX_HEIGHT - 1) - (int)apf_to_long(&norm);
        }
        
        if (apf_cmp(xmin, &zero_apf) <= 0 && apf_cmp(xmax, &zero_apf) >= 0) {
            apf_sub(&temp, &zero_apf, xmin);
            apf_div(&norm, &temp, &x_range);
            apf_from_int(&temp, GFX_WIDTH - 1);
            apf_mul(&norm, &norm, &temp);
            axis_x = (int)apf_to_long(&norm);
        }
    }
    
    set_video_mode(0x13);
    
    if (axis_y >= 0 && axis_y < GFX_HEIGHT) {
        draw_hline(0, GFX_WIDTH - 1, axis_y, COLOR_CYAN);
    }
    if (axis_x >= 0 && axis_x < GFX_WIDTH) {
        draw_vline(axis_x, 0, GFX_HEIGHT - 1, COLOR_CYAN);
    }
    
    prev_x = -1;
    prev_y = -1;
    first_point = 1;
    
    for (x = 0; x < GFX_WIDTH; x++) {
        apf x_val, idx_apf, norm_y;
        int py, eval_ok;
        
        if (kbhit()) {
            getch();
            break;
        }
        
        apf_from_int(&idx_apf, x);
        apf_mul(&temp, &idx_apf, &x_step);
        apf_add(&x_val, xmin, &temp);
        
        apf_copy(&arg.re, &x_val);
        apf_zero(&arg.im);
        
        if (use_expr) {
            eval_ok = eval_expr_with_var(&result, expr, var, &arg);
        } else {
            eval_ok = eval_user_func(&result, func_idx, &arg);
        }
        
        if (eval_ok && !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
            apf_iszero(&result.im)) {
            
            apf_sub(&temp, &result.re, &ymin);
            apf_div(&norm_y, &temp, &y_range);
            apf_from_int(&temp, GFX_HEIGHT - 1);
            apf_mul(&norm_y, &norm_y, &temp);
            py = (GFX_HEIGHT - 1) - (int)apf_to_long(&norm_y);
            
            if (apf_cmp(&result.re, &ymin) < 0 || apf_cmp(&result.re, &ymax) > 0) {
                if (py < 0) py = -1;
                if (py >= GFX_HEIGHT) py = -1;
            }
            
            if (py >= 0 && py < GFX_HEIGHT) {
                if (!first_point && prev_y >= 0) {
                    draw_line(prev_x, prev_y, x, py, COLOR_YELLOW);
                } else {
                    bios_put_pixel(x, py, COLOR_YELLOW);
                }
                prev_x = x;
                prev_y = py;
                first_point = 0;
            } else {
                first_point = 1;
            }
        } else {
            first_point = 1;
        }
    }
    
    while (kbhit()) getch();
    getch();
    
    set_video_mode(0x03);
    dos_init_screen();
    
    return 1;
}

/* Text mode plot for DOS */
static int text_plot_impl(const char *expr, char var, apf *xmin, apf *xmax)
{
    int func_idx, x, row;
    apf x_range, x_step, x_val, temp;
    apfc arg, result;
    apf ymin, ymax, y_range;
    static apf y_vals[80];
    static int valid[80];
    static int y_screen[80];
    char line[82];
    int width = 70, height = 20, axis_row = -1;
    int first_valid = 0;
    int use_expr = 0;
    
    func_idx = get_func_index(expr);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        use_expr = 1;
    }
    
    apf_sub(&x_range, xmax, xmin);
    apf_from_int(&temp, width);
    apf_div(&x_step, &x_range, &temp);
    
    for (x = 0; x < width; x++) {
        apf screen_x;
        int eval_ok;
        apf_from_int(&screen_x, x);
        apf_mul(&temp, &screen_x, &x_step);
        apf_add(&x_val, xmin, &temp);
        
        apf_copy(&arg.re, &x_val);
        apf_zero(&arg.im);
        
        if (use_expr) {
            eval_ok = eval_expr_with_var(&result, expr, var, &arg);
        } else {
            eval_ok = eval_user_func(&result, func_idx, &arg);
        }
        
        if (eval_ok &&
            !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
            apf_iszero(&result.im)) {
            apf_copy(&y_vals[x], &result.re);
            valid[x] = 1;
            if (!first_valid || apf_cmp(&result.re, &ymin) < 0) apf_copy(&ymin, &result.re);
            if (!first_valid || apf_cmp(&result.re, &ymax) > 0) apf_copy(&ymax, &result.re);
            first_valid = 1;
        } else {
            valid[x] = 0;
        }
    }
    
    if (!first_valid) {
        printf("No valid points in range\n");
        return 0;
    }
    
    apf_sub(&y_range, &ymax, &ymin);
    if (apf_iszero(&y_range)) apf_from_int(&y_range, 1);
    {
        apf margin, ten;
        apf_from_int(&ten, 10);
        apf_div(&margin, &y_range, &ten);
        apf_sub(&ymin, &ymin, &margin);
        apf_add(&ymax, &ymax, &margin);
        apf_sub(&y_range, &ymax, &ymin);
    }
    
    for (x = 0; x < width; x++) {
        if (valid[x]) {
            apf norm_y;
            apf_sub(&temp, &y_vals[x], &ymin);
            apf_div(&norm_y, &temp, &y_range);
            apf_from_int(&temp, height - 1);
            apf_mul(&norm_y, &norm_y, &temp);
            y_screen[x] = (height - 1) - (int)apf_to_long(&norm_y);
        }
    }
    
    {
        apf zero_apf, norm;
        apf_zero(&zero_apf);
        if (apf_cmp(&ymin, &zero_apf) <= 0 && apf_cmp(&ymax, &zero_apf) >= 0) {
            apf_sub(&temp, &zero_apf, &ymin);
            apf_div(&norm, &temp, &y_range);
            apf_from_int(&temp, height - 1);
            apf_mul(&norm, &norm, &temp);
            axis_row = (height - 1) - (int)apf_to_long(&norm);
        }
    }
    
    printf("\n");
    for (row = 0; row < height; row++) {
        int col;
        for (col = 0; col < width; col++) {
            if (valid[col] && y_screen[col] == row) line[col] = '*';
            else if (row == axis_row) line[col] = '-';
            else line[col] = ' ';
        }
        line[width] = '\0';
        if (row == 0) printf("%8ld |%s\n", apf_to_long(&ymax), line);
        else if (row == height - 1) printf("%8ld |%s\n", apf_to_long(&ymin), line);
        else if (row == axis_row) printf("%8d |%s\n", 0, line);
        else printf("         |%s\n", line);
    }
    
    printf("         +");
    for (x = 0; x < width; x++) printf("-");
    printf("\n         %-8ld", apf_to_long(xmin));
    for (x = 8; x < width - 8; x++) printf(" ");
    printf("%8ld\n\n", apf_to_long(xmax));
    
    return 1;
}

int do_textplot(const char *expr, char var, apf *xmin, apf *xmax)
{
    return text_plot_impl(expr, var, xmin, xmax);
}

#else
/* Non-DOS: ASCII plot only */

int do_plot(const char *expr, char var, apf *xmin, apf *xmax)
{
    int func_idx, x, row;
    apf x_range, x_step, x_val, temp;
    apfc arg, result;
    apf ymin, ymax, y_range;
    static apf y_vals[80];
    static int valid[80];
    static int y_screen[80];
    char line[82];
    int width = 70, height = 20, axis_row = -1;
    int first_valid = 0;
    int use_expr = 0;
    
    func_idx = get_func_index(expr);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        use_expr = 1;
    }
    
    apf_sub(&x_range, xmax, xmin);
    apf_from_int(&temp, width);
    apf_div(&x_step, &x_range, &temp);
    
    for (x = 0; x < width; x++) {
        apf screen_x;
        apf_from_int(&screen_x, x);
        apf_mul(&temp, &screen_x, &x_step);
        apf_add(&x_val, xmin, &temp);
        
        apf_copy(&arg.re, &x_val);
        apf_zero(&arg.im);
        
        if (use_expr) {
            if (eval_expr_with_var(&result, expr, var, &arg) &&
                !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
                apf_iszero(&result.im)) {
                apf_copy(&y_vals[x], &result.re);
                valid[x] = 1;
                if (!first_valid || apf_cmp(&result.re, &ymin) < 0) apf_copy(&ymin, &result.re);
                if (!first_valid || apf_cmp(&result.re, &ymax) > 0) apf_copy(&ymax, &result.re);
                first_valid = 1;
            } else {
                valid[x] = 0;
            }
        } else {
            if (eval_user_func(&result, func_idx, &arg) &&
                !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
                apf_iszero(&result.im)) {
                apf_copy(&y_vals[x], &result.re);
                valid[x] = 1;
                if (!first_valid || apf_cmp(&result.re, &ymin) < 0) apf_copy(&ymin, &result.re);
                if (!first_valid || apf_cmp(&result.re, &ymax) > 0) apf_copy(&ymax, &result.re);
                first_valid = 1;
            } else {
                valid[x] = 0;
            }
        }
    }
    
    if (!first_valid) {
        printf("No valid points in range\n");
        return 0;
    }
    
    apf_sub(&y_range, &ymax, &ymin);
    if (apf_iszero(&y_range)) apf_from_int(&y_range, 1);
    {
        apf margin, ten;
        apf_from_int(&ten, 10);
        apf_div(&margin, &y_range, &ten);
        apf_sub(&ymin, &ymin, &margin);
        apf_add(&ymax, &ymax, &margin);
        apf_sub(&y_range, &ymax, &ymin);
    }
    
    for (x = 0; x < width; x++) {
        if (valid[x]) {
            apf norm_y;
            apf_sub(&temp, &y_vals[x], &ymin);
            apf_div(&norm_y, &temp, &y_range);
            apf_from_int(&temp, height - 1);
            apf_mul(&norm_y, &norm_y, &temp);
            y_screen[x] = (height - 1) - (int)apf_to_long(&norm_y);
        }
    }
    
    {
        apf zero_apf, norm;
        apf_zero(&zero_apf);
        if (apf_cmp(&ymin, &zero_apf) <= 0 && apf_cmp(&ymax, &zero_apf) >= 0) {
            apf_sub(&temp, &zero_apf, &ymin);
            apf_div(&norm, &temp, &y_range);
            apf_from_int(&temp, height - 1);
            apf_mul(&norm, &norm, &temp);
            axis_row = (height - 1) - (int)apf_to_long(&norm);
        }
    }
    
    printf("\n");
    for (row = 0; row < height; row++) {
        int col;
        for (col = 0; col < width; col++) {
            if (valid[col] && y_screen[col] == row) line[col] = '*';
            else if (row == axis_row) line[col] = '-';
            else line[col] = ' ';
        }
        line[width] = '\0';
        if (row == 0) printf("%8ld |%s\n", apf_to_long(&ymax), line);
        else if (row == height - 1) printf("%8ld |%s\n", apf_to_long(&ymin), line);
        else if (row == axis_row) printf("%8d |%s\n", 0, line);
        else printf("         |%s\n", line);
    }
    
    printf("         +");
    for (x = 0; x < width; x++) printf("-");
    printf("\n         %-8ld", apf_to_long(xmin));
    for (x = 8; x < width - 8; x++) printf(" ");
    printf("%8ld\n\n", apf_to_long(xmax));
    
    return 1;
}

int do_textplot(const char *expr, char var, apf *xmin, apf *xmax)
{
    return do_plot(expr, var, xmin, xmax);
}

#endif /* HAVE_CONIO */

/* ========== Chaotic Systems - 100% APF implementation ========== */

/*
 * Lorenz Attractor - the classic chaotic system
 * dx/dt = sigma * (y - x)
 * dy/dt = x * (rho - z) - y
 * dz/dt = x * y - beta * z
 * Classic parameters: sigma=10, rho=28, beta=8/3
 */
int do_lorenz_text(long sigma_l, long rho_l, long beta_num, long beta_den, int steps)
{
    apf x, y, z, sigma, rho, beta;
    apf dt, half_dt;
    apf xmin_a, xmax_a, zmin_a, zmax_a;
    int width = 70, height = 24;
    static char screen[30][80];
    int i, row, col;
    
    /* Initialize parameters */
    apf_from_int(&sigma, sigma_l);
    apf_from_int(&rho, rho_l);
    {
        apf num, den;
        apf_from_int(&num, beta_num);
        apf_from_int(&den, beta_den);
        apf_div(&beta, &num, &den);
    }
    
    /* dt = 0.01 */
    {
        apf num, den;
        apf_from_int(&num, 1);
        apf_from_int(&den, 100);
        apf_div(&dt, &num, &den);
        apf_from_int(&den, 2);
        apf_div(&half_dt, &dt, &den);
    }
    
    /* Initial conditions */
    apf_from_int(&x, 1);
    apf_from_int(&y, 1);
    apf_from_int(&z, 1);
    
    /* Plot bounds */
    apf_from_int(&xmin_a, -25);
    apf_from_int(&xmax_a, 25);
    apf_from_int(&zmin_a, 0);
    apf_from_int(&zmax_a, 55);
    
    /* Clear screen buffer */
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
        screen[row][width] = '\0';
    }
    
    printf("\nLorenz Attractor (sigma=%ld, rho=%ld, beta=%ld/%ld)\n", 
           sigma_l, rho_l, beta_num, beta_den);
    printf("Projecting x-z plane, %d iterations\n\n", steps);
    
    /* Integrate using RK4 */
    for (i = 0; i < steps; i++) {
        apf k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
        apf tx, ty, tz, tmp, tmp2;
        apf dx_final, dy_final, dz_final;
        int px, pz;
        
        /* k1 = f(x, y, z) */
        apf_sub(&tmp, &y, &x);
        apf_mul(&k1x, &sigma, &tmp);            /* k1x = sigma * (y - x) */
        
        apf_sub(&tmp, &rho, &z);
        apf_mul(&tmp2, &x, &tmp);
        apf_sub(&k1y, &tmp2, &y);               /* k1y = x * (rho - z) - y */
        
        apf_mul(&tmp, &x, &y);
        apf_mul(&tmp2, &beta, &z);
        apf_sub(&k1z, &tmp, &tmp2);             /* k1z = x * y - beta * z */
        
        /* k2 = f(x + 0.5*dt*k1, ...) */
        apf_mul(&tmp, &half_dt, &k1x);
        apf_add(&tx, &x, &tmp);
        apf_mul(&tmp, &half_dt, &k1y);
        apf_add(&ty, &y, &tmp);
        apf_mul(&tmp, &half_dt, &k1z);
        apf_add(&tz, &z, &tmp);
        
        apf_sub(&tmp, &ty, &tx);
        apf_mul(&k2x, &sigma, &tmp);
        apf_sub(&tmp, &rho, &tz);
        apf_mul(&tmp2, &tx, &tmp);
        apf_sub(&k2y, &tmp2, &ty);
        apf_mul(&tmp, &tx, &ty);
        apf_mul(&tmp2, &beta, &tz);
        apf_sub(&k2z, &tmp, &tmp2);
        
        /* k3 = f(x + 0.5*dt*k2, ...) */
        apf_mul(&tmp, &half_dt, &k2x);
        apf_add(&tx, &x, &tmp);
        apf_mul(&tmp, &half_dt, &k2y);
        apf_add(&ty, &y, &tmp);
        apf_mul(&tmp, &half_dt, &k2z);
        apf_add(&tz, &z, &tmp);
        
        apf_sub(&tmp, &ty, &tx);
        apf_mul(&k3x, &sigma, &tmp);
        apf_sub(&tmp, &rho, &tz);
        apf_mul(&tmp2, &tx, &tmp);
        apf_sub(&k3y, &tmp2, &ty);
        apf_mul(&tmp, &tx, &ty);
        apf_mul(&tmp2, &beta, &tz);
        apf_sub(&k3z, &tmp, &tmp2);
        
        /* k4 = f(x + dt*k3, ...) */
        apf_mul(&tmp, &dt, &k3x);
        apf_add(&tx, &x, &tmp);
        apf_mul(&tmp, &dt, &k3y);
        apf_add(&ty, &y, &tmp);
        apf_mul(&tmp, &dt, &k3z);
        apf_add(&tz, &z, &tmp);
        
        apf_sub(&tmp, &ty, &tx);
        apf_mul(&k4x, &sigma, &tmp);
        apf_sub(&tmp, &rho, &tz);
        apf_mul(&tmp2, &tx, &tmp);
        apf_sub(&k4y, &tmp2, &ty);
        apf_mul(&tmp, &tx, &ty);
        apf_mul(&tmp2, &beta, &tz);
        apf_sub(&k4z, &tmp, &tmp2);
        
        /* dx = (k1 + 2*k2 + 2*k3 + k4) / 6 */
        {
            apf two, six;
            apf_from_int(&two, 2);
            apf_from_int(&six, 6);
            
            apf_mul(&tmp, &two, &k2x);
            apf_add(&dx_final, &k1x, &tmp);
            apf_mul(&tmp, &two, &k3x);
            apf_add(&dx_final, &dx_final, &tmp);
            apf_add(&dx_final, &dx_final, &k4x);
            apf_div(&dx_final, &dx_final, &six);
            
            apf_mul(&tmp, &two, &k2y);
            apf_add(&dy_final, &k1y, &tmp);
            apf_mul(&tmp, &two, &k3y);
            apf_add(&dy_final, &dy_final, &tmp);
            apf_add(&dy_final, &dy_final, &k4y);
            apf_div(&dy_final, &dy_final, &six);
            
            apf_mul(&tmp, &two, &k2z);
            apf_add(&dz_final, &k1z, &tmp);
            apf_mul(&tmp, &two, &k3z);
            apf_add(&dz_final, &dz_final, &tmp);
            apf_add(&dz_final, &dz_final, &k4z);
            apf_div(&dz_final, &dz_final, &six);
        }
        
        /* Update: x += dt * dx */
        apf_mul(&tmp, &dt, &dx_final);
        apf_add(&x, &x, &tmp);
        apf_mul(&tmp, &dt, &dy_final);
        apf_add(&y, &y, &tmp);
        apf_mul(&tmp, &dt, &dz_final);
        apf_add(&z, &z, &tmp);
        
        /* Map to screen (x-z projection) */
        {
            apf range, norm;
            apf_sub(&range, &xmax_a, &xmin_a);
            apf_sub(&tmp, &x, &xmin_a);
            apf_div(&norm, &tmp, &range);
            apf_from_int(&tmp, width - 1);
            apf_mul(&norm, &norm, &tmp);
            px = (int)apf_to_long(&norm);
            
            apf_sub(&range, &zmax_a, &zmin_a);
            apf_sub(&tmp, &z, &zmin_a);
            apf_div(&norm, &tmp, &range);
            apf_from_int(&tmp, height - 1);
            apf_mul(&norm, &norm, &tmp);
            pz = (height - 1) - (int)apf_to_long(&norm);
        }
        
        if (px >= 0 && px < width && pz >= 0 && pz < height) {
            if (screen[pz][px] == ' ') screen[pz][px] = '.';
            else if (screen[pz][px] == '.') screen[pz][px] = ':';
            else if (screen[pz][px] == ':') screen[pz][px] = '+';
            else if (screen[pz][px] == '+') screen[pz][px] = '*';
            else if (screen[pz][px] == '*') screen[pz][px] = '#';
        }
    }
    
    /* Print result */
    printf("z=55 |");
    for (col = 0; col < width; col++) printf("%c", screen[0][col]);
    printf("\n");
    
    for (row = 1; row < height - 1; row++) {
        printf("     |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("z=0  |");
    for (col = 0; col < width; col++) printf("%c", screen[height-1][col]);
    printf("\n");
    
    printf("     +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("     x=-25");
    for (col = 0; col < width - 12; col++) printf(" ");
    printf("x=25\n\n");
    
    return 1;
}

/* Rossler Attractor using APF
 * dx/dt = -y - z
 * dy/dt = x + a*y
 * dz/dt = b + z*(x - c)
 */
int do_rossler_text(long a_num, long a_den, long b_num, long b_den, 
                    long c_num, long c_den, int steps)
{
    apf x, y, z, a, b, c;
    apf dt;
    apf xmin_a, xmax_a, ymin_a, ymax_a;
    int width = 70, height = 24;
    static char screen[30][80];
    int i, row, col;
    
    /* Initialize parameters */
    {
        apf num, den;
        apf_from_int(&num, a_num);
        apf_from_int(&den, a_den);
        apf_div(&a, &num, &den);
        apf_from_int(&num, b_num);
        apf_from_int(&den, b_den);
        apf_div(&b, &num, &den);
        apf_from_int(&num, c_num);
        apf_from_int(&den, c_den);
        apf_div(&c, &num, &den);
    }
    
    /* dt = 0.02 */
    {
        apf num, den;
        apf_from_int(&num, 1);
        apf_from_int(&den, 50);
        apf_div(&dt, &num, &den);
    }
    
    apf_from_int(&x, 1);
    apf_from_int(&y, 1);
    apf_from_int(&z, 1);
    
    apf_from_int(&xmin_a, -15);
    apf_from_int(&xmax_a, 15);
    apf_from_int(&ymin_a, -15);
    apf_from_int(&ymax_a, 15);
    
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
        screen[row][width] = '\0';
    }
    
    printf("\nRossler Attractor (a=%ld/%ld, b=%ld/%ld, c=%ld/%ld)\n",
           a_num, a_den, b_num, b_den, c_num, c_den);
    printf("Projecting x-y plane, %d iterations\n\n", steps);
    
    /* Simple Euler integration */
    for (i = 0; i < steps; i++) {
        apf dx, dy, dz, tmp, tmp2;
        int px, py;
        
        /* dx = -y - z */
        apf_add(&tmp, &y, &z);
        apf_neg(&dx, &tmp);
        
        /* dy = x + a*y */
        apf_mul(&tmp, &a, &y);
        apf_add(&dy, &x, &tmp);
        
        /* dz = b + z*(x - c) */
        apf_sub(&tmp, &x, &c);
        apf_mul(&tmp2, &z, &tmp);
        apf_add(&dz, &b, &tmp2);
        
        /* Update */
        apf_mul(&tmp, &dt, &dx);
        apf_add(&x, &x, &tmp);
        apf_mul(&tmp, &dt, &dy);
        apf_add(&y, &y, &tmp);
        apf_mul(&tmp, &dt, &dz);
        apf_add(&z, &z, &tmp);
        
        /* Map to screen */
        {
            apf range, norm;
            apf_sub(&range, &xmax_a, &xmin_a);
            apf_sub(&tmp, &x, &xmin_a);
            apf_div(&norm, &tmp, &range);
            apf_from_int(&tmp, width - 1);
            apf_mul(&norm, &norm, &tmp);
            px = (int)apf_to_long(&norm);
            
            apf_sub(&range, &ymax_a, &ymin_a);
            apf_sub(&tmp, &y, &ymin_a);
            apf_div(&norm, &tmp, &range);
            apf_from_int(&tmp, height - 1);
            apf_mul(&norm, &norm, &tmp);
            py = (height - 1) - (int)apf_to_long(&norm);
        }
        
        if (px >= 0 && px < width && py >= 0 && py < height) {
            if (screen[py][px] == ' ') screen[py][px] = '.';
            else if (screen[py][px] == '.') screen[py][px] = ':';
            else if (screen[py][px] == ':') screen[py][px] = '+';
            else if (screen[py][px] == '+') screen[py][px] = '*';
            else if (screen[py][px] == '*') screen[py][px] = '#';
        }
    }
    
    /* Print */
    printf("y=15 |");
    for (col = 0; col < width; col++) printf("%c", screen[0][col]);
    printf("\n");
    
    for (row = 1; row < height - 1; row++) {
        printf("     |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("y=-15|");
    for (col = 0; col < width; col++) printf("%c", screen[height-1][col]);
    printf("\n");
    
    printf("     +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("     x=-15");
    for (col = 0; col < width - 12; col++) printf(" ");
    printf("x=15\n\n");
    
    return 1;
}

/* Parametric curve plotting using APF */
int do_parametric_text(const char *xfunc, const char *yfunc, 
                       apf *tmin, apf *tmax, int steps)
{
    int xidx, yidx;
    apf t, dt, t_range;
    apf xmin_a, xmax_a, ymin_a, ymax_a;
    int width = 70, height = 24;
    static char screen[30][80];
    /* Store points using integers (screen coords) to save memory */
    static int px_arr[500], py_arr[500];
    int i, row, col, valid_count = 0;
    int first_valid = 0;
    
    if (steps > 500) steps = 500;
    
    xidx = get_func_index(xfunc);
    yidx = get_func_index(yfunc);
    
    if (xidx < 0 || !user_funcs[xidx].defined) {
        printf("Error: function '%s' not defined\n", xfunc);
        return 0;
    }
    if (yidx < 0 || !user_funcs[yidx].defined) {
        printf("Error: function '%s' not defined\n", yfunc);
        return 0;
    }
    
    apf_sub(&t_range, tmax, tmin);
    {
        apf steps_apf;
        apf_from_int(&steps_apf, steps);
        apf_div(&dt, &t_range, &steps_apf);
    }
    
    /* First pass: find range */
    apf_copy(&t, tmin);
    for (i = 0; i < steps; i++) {
        apfc tval, xres, yres;
        
        apf_copy(&tval.re, &t);
        apf_zero(&tval.im);
        
        if (eval_user_func(&xres, xidx, &tval) && 
            eval_user_func(&yres, yidx, &tval) &&
            !apf_isnan(&xres.re) && !apf_isnan(&yres.re) &&
            !apf_isinf(&xres.re) && !apf_isinf(&yres.re)) {
            
            if (!first_valid) {
                apf_copy(&xmin_a, &xres.re);
                apf_copy(&xmax_a, &xres.re);
                apf_copy(&ymin_a, &yres.re);
                apf_copy(&ymax_a, &yres.re);
                first_valid = 1;
            } else {
                if (apf_cmp(&xres.re, &xmin_a) < 0) apf_copy(&xmin_a, &xres.re);
                if (apf_cmp(&xres.re, &xmax_a) > 0) apf_copy(&xmax_a, &xres.re);
                if (apf_cmp(&yres.re, &ymin_a) < 0) apf_copy(&ymin_a, &yres.re);
                if (apf_cmp(&yres.re, &ymax_a) > 0) apf_copy(&ymax_a, &yres.re);
            }
            px_arr[i] = 1;  /* Mark valid */
            valid_count++;
        } else {
            px_arr[i] = -1;  /* Mark invalid */
        }
        
        apf_add(&t, &t, &dt);
    }
    
    if (valid_count == 0) {
        printf("No valid points computed\n");
        return 0;
    }
    
    /* Add margin */
    {
        apf margin, twenty, x_range, y_range;
        apf_from_int(&twenty, 20);
        
        apf_sub(&x_range, &xmax_a, &xmin_a);
        apf_div(&margin, &x_range, &twenty);
        apf_sub(&xmin_a, &xmin_a, &margin);
        apf_add(&xmax_a, &xmax_a, &margin);
        
        apf_sub(&y_range, &ymax_a, &ymin_a);
        apf_div(&margin, &y_range, &twenty);
        apf_sub(&ymin_a, &ymin_a, &margin);
        apf_add(&ymax_a, &ymax_a, &margin);
    }
    
    /* Second pass: compute screen coords */
    apf_copy(&t, tmin);
    for (i = 0; i < steps; i++) {
        if (px_arr[i] == 1) {
            apfc tval, xres, yres;
            apf tmp, norm, x_range, y_range;
            
            apf_copy(&tval.re, &t);
            apf_zero(&tval.im);
            
            eval_user_func(&xres, xidx, &tval);
            eval_user_func(&yres, yidx, &tval);
            
            apf_sub(&x_range, &xmax_a, &xmin_a);
            apf_sub(&tmp, &xres.re, &xmin_a);
            apf_div(&norm, &tmp, &x_range);
            apf_from_int(&tmp, width - 1);
            apf_mul(&norm, &norm, &tmp);
            px_arr[i] = (int)apf_to_long(&norm);
            
            apf_sub(&y_range, &ymax_a, &ymin_a);
            apf_sub(&tmp, &yres.re, &ymin_a);
            apf_div(&norm, &tmp, &y_range);
            apf_from_int(&tmp, height - 1);
            apf_mul(&norm, &norm, &tmp);
            py_arr[i] = (height - 1) - (int)apf_to_long(&norm);
        } else {
            px_arr[i] = -1;
            py_arr[i] = -1;
        }
        
        apf_add(&t, &t, &dt);
    }
    
    /* Clear screen and plot */
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
    }
    
    printf("\nParametric: %s(t), %s(t)  t=[%ld, %ld]\n\n", 
           xfunc, yfunc, apf_to_long(tmin), apf_to_long(tmax));
    
    for (i = 0; i < steps; i++) {
        int px = px_arr[i];
        int py = py_arr[i];
        if (px >= 0 && px < width && py >= 0 && py < height) {
            if (screen[py][px] == ' ') screen[py][px] = '.';
            else if (screen[py][px] == '.') screen[py][px] = '*';
            else screen[py][px] = '#';
        }
    }
    
    /* Print */
    for (row = 0; row < height; row++) {
        if (row == 0) printf("%7ld |", apf_to_long(&ymax_a));
        else if (row == height - 1) printf("%7ld |", apf_to_long(&ymin_a));
        else printf("        |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("        +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("        %-7ld", apf_to_long(&xmin_a));
    for (col = 0; col < width - 14; col++) printf(" ");
    printf("%7ld\n\n", apf_to_long(&xmax_a));
    
    return 1;
}

/* ========================================================================
 * Interactive Plot (iplot) - Full-screen ANSI terminal plotting
 * ======================================================================== */


/* ========================================================================
 * Interactive Plot (iplot) and Mandelbrot Viewer
 * Full-screen ANSI terminal plotting with pan/zoom
 * ======================================================================== */

#if !defined(HAVE_CONIO) && !defined(_WIN32)
/* Unix terminal handling */
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <signal.h>

static struct termios iplot_orig_termios;
static int iplot_term_setup = 0;

static void iplot_restore_term(void)
{
    if (iplot_term_setup) {
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &iplot_orig_termios);
        iplot_term_setup = 0;
    }
    /* Show cursor, exit alternate screen */
    printf("\033[?25h\033[?1049l");
    fflush(stdout);
}

static void iplot_setup_term(void)
{
    struct termios raw;
    
    if (!isatty(STDIN_FILENO)) return;
    
    tcgetattr(STDIN_FILENO, &iplot_orig_termios);
    raw = iplot_orig_termios;
    raw.c_lflag &= ~(ECHO | ICANON);
    raw.c_cc[VMIN] = 1;
    raw.c_cc[VTIME] = 0;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
    iplot_term_setup = 1;
    
    /* Enter alternate screen, hide cursor */
    printf("\033[?1049h\033[?25l");
    fflush(stdout);
}

static void iplot_get_terminal_size(int *width, int *height)
{
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0) {
        *width = w.ws_col;
        *height = w.ws_row;
    } else {
        *width = 80;
        *height = 24;
    }
}

static int iplot_read_key(void)
{
    unsigned char c;
    if (read(STDIN_FILENO, &c, 1) != 1) return -1;
    
    if (c == 27) {  /* Escape sequence */
        unsigned char seq[3];
        if (read(STDIN_FILENO, &seq[0], 1) != 1) return 27;  /* Just ESC */
        if (seq[0] != '[') return 27;
        if (read(STDIN_FILENO, &seq[1], 1) != 1) return 27;
        
        switch (seq[1]) {
            case 'A': return 1001;  /* Up */
            case 'B': return 1002;  /* Down */
            case 'C': return 1003;  /* Right */
            case 'D': return 1004;  /* Left */
        }
        return 27;
    }
    return c;
}

int do_iplot(const char *expr, char var, apf *xmin_init, apf *xmax_init)
{
    int width, height;
    int plot_width, plot_height;
    int x, y, i;
    apf xmin, xmax, ymin, ymax;
    apf ymin_fixed, ymax_fixed;  /* User-controlled y range */
    apf x_range, y_range, x_step;
    apf *y_values;
    int *valid;
    char *screen;
    int running = 1;
    int need_redraw = 1;
    int func_idx;
    int use_expr = 0;
    int auto_y = 1;  /* Auto-scale y axis */
    const int label_width = 8;  /* Fixed width for Y axis labels */
    
    /* Check function */
    func_idx = get_func_index(expr);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        use_expr = 1;
    }
    
    /* Copy initial range - default to -10:10 if not specified */
    if (apf_eq(xmin_init, xmax_init)) {
        apf_from_int(&xmin, -10);
        apf_from_int(&xmax, 10);
    } else {
        apf_copy(&xmin, xmin_init);
        apf_copy(&xmax, xmax_init);
    }
    
    /* Default y range (will be auto-scaled initially) */
    apf_from_int(&ymin_fixed, -10);
    apf_from_int(&ymax_fixed, 10);
    
    /* Setup terminal */
    iplot_setup_term();
    iplot_get_terminal_size(&width, &height);
    
    /* Reserve space for axes labels (label_width + 1 for border) */
    plot_width = width - label_width - 2;
    plot_height = height - 4;
    
    /* Adjust for character aspect ratio: chars are ~2x tall as wide
     * So we need to scale x by 0.5 (or equivalently, use half the width for same visual range)
     * This is handled by making the x_range effectively cover more "visual" space */
    
    if (plot_width < 20 || plot_height < 10) {
        iplot_restore_term();
        printf("Terminal too small for iplot\n");
        return 0;
    }
    
    /* Allocate buffers */
    y_values = (apf *)malloc(plot_width * sizeof(apf));
    valid = (int *)malloc(plot_width * sizeof(int));
    screen = (char *)malloc(plot_width * plot_height);
    
    if (!y_values || !valid || !screen) {
        iplot_restore_term();
        printf("Out of memory\n");
        free(y_values); free(valid); free(screen);
        return 0;
    }
    
    while (running) {
        if (need_redraw) {
            int axis_x = -1, axis_y = -1;
            int first_valid = 1;
            apf temp, idx_apf;
            apf auto_ymin, auto_ymax;
            
            /* Calculate x step - scale for aspect ratio (chars ~2x tall as wide) */
            apf_sub(&x_range, &xmax, &xmin);
            apf_from_int(&temp, plot_width - 1);
            apf_div(&x_step, &x_range, &temp);
            
            /* Evaluate function at each x */
            for (x = 0; x < plot_width; x++) {
                apf x_val;
                apfc arg, result;
                int eval_ok;
                
                apf_from_int(&idx_apf, x);
                apf_mul(&temp, &idx_apf, &x_step);
                apf_add(&x_val, &xmin, &temp);
                
                apf_copy(&arg.re, &x_val);
                apf_zero(&arg.im);
                
                if (use_expr) {
                    eval_ok = eval_expr_with_var(&result, expr, var, &arg);
                } else {
                    eval_ok = eval_user_func(&result, func_idx, &arg);
                }
                
                if (eval_ok && !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
                    apf_iszero(&result.im)) {
                    valid[x] = 1;
                    apf_copy(&y_values[x], &result.re);
                    
                    if (first_valid) {
                        apf_copy(&auto_ymin, &result.re);
                        apf_copy(&auto_ymax, &result.re);
                        first_valid = 0;
                    } else {
                        if (apf_cmp(&result.re, &auto_ymin) < 0) apf_copy(&auto_ymin, &result.re);
                        if (apf_cmp(&result.re, &auto_ymax) > 0) apf_copy(&auto_ymax, &result.re);
                    }
                } else {
                    valid[x] = 0;
                }
            }
            
            /* Determine y range */
            if (auto_y) {
                if (first_valid) {
                    apf_from_int(&ymin, -1);
                    apf_from_int(&ymax, 1);
                } else {
                    apf_copy(&ymin, &auto_ymin);
                    apf_copy(&ymax, &auto_ymax);
                }
                
                /* Add 5% margin to y */
                apf_sub(&y_range, &ymax, &ymin);
                if (apf_iszero(&y_range)) apf_from_int(&y_range, 2);
                {
                    apf margin, twenty;
                    apf_from_int(&twenty, 20);
                    apf_div(&margin, &y_range, &twenty);
                    apf_sub(&ymin, &ymin, &margin);
                    apf_add(&ymax, &ymax, &margin);
                }
            } else {
                apf_copy(&ymin, &ymin_fixed);
                apf_copy(&ymax, &ymax_fixed);
            }
            
            apf_sub(&y_range, &ymax, &ymin);
            
            /* Find axis positions */
            {
                apf zero_apf, norm;
                apf_zero(&zero_apf);
                
                if (apf_cmp(&ymin, &zero_apf) <= 0 && apf_cmp(&ymax, &zero_apf) >= 0) {
                    apf_sub(&temp, &zero_apf, &ymin);
                    apf_div(&norm, &temp, &y_range);
                    apf_from_int(&temp, plot_height - 1);
                    apf_mul(&norm, &norm, &temp);
                    axis_y = (plot_height - 1) - (int)apf_to_long(&norm);
                }
                
                if (apf_cmp(&xmin, &zero_apf) <= 0 && apf_cmp(&xmax, &zero_apf) >= 0) {
                    apf_sub(&temp, &zero_apf, &xmin);
                    apf_div(&norm, &temp, &x_range);
                    apf_from_int(&temp, plot_width - 1);
                    apf_mul(&norm, &norm, &temp);
                    axis_x = (int)apf_to_long(&norm);
                }
            }
            
            /* Clear screen buffer */
            memset(screen, ' ', plot_width * plot_height);
            
            /* Draw axes */
            if (axis_y >= 0 && axis_y < plot_height) {
                for (x = 0; x < plot_width; x++) {
                    screen[axis_y * plot_width + x] = '-';
                }
            }
            if (axis_x >= 0 && axis_x < plot_width) {
                for (y = 0; y < plot_height; y++) {
                    char c = screen[y * plot_width + axis_x];
                    screen[y * plot_width + axis_x] = (c == '-') ? '+' : '|';
                }
            }
            
            /* Plot points */
            for (x = 0; x < plot_width; x++) {
                if (valid[x]) {
                    apf norm;
                    int py;
                    
                    apf_sub(&temp, &y_values[x], &ymin);
                    apf_div(&norm, &temp, &y_range);
                    apf_from_int(&temp, plot_height - 1);
                    apf_mul(&norm, &norm, &temp);
                    py = (plot_height - 1) - (int)apf_to_long(&norm);
                    
                    if (py >= 0 && py < plot_height) {
                        screen[py * plot_width + x] = '*';
                    }
                }
            }
            
            /* Render to terminal */
            printf("\033[H\033[2J");  /* Clear screen, home cursor */
            
            /* Title */
            printf("\033[1;36m iplot: \033[1;33m%s\033[0m", expr);
            printf("  \033[90m[arrows=pan +/-=zoom a=auto-y r=reset q=quit]\033[0m\n");
            
            /* Top border with y-max */
            {
                double ymax_d = apf_to_double(&ymax);
                printf("\033[33m%7.2g\033[0m \033[90m+", ymax_d);
                for (x = 0; x < plot_width; x++) printf("-");
                printf("+\033[0m\n");
            }
            
            /* Plot area */
            for (y = 0; y < plot_height; y++) {
                if (y == plot_height / 2) {
                    double ymid = apf_to_double(&ymin) + apf_to_double(&y_range) / 2;
                    printf("\033[33m%7.2g\033[0m \033[90m|\033[0m", ymid);
                } else {
                    printf("        \033[90m|\033[0m");
                }
                
                for (x = 0; x < plot_width; x++) {
                    char c = screen[y * plot_width + x];
                    if (c == '*') {
                        printf("\033[1;32m*\033[0m");
                    } else if (c == '-' || c == '|' || c == '+') {
                        printf("\033[90m%c\033[0m", c);
                    } else {
                        putchar(c);
                    }
                }
                printf("\033[90m|\033[0m\n");
            }
            
            /* Bottom border with y-min */
            {
                double ymin_d = apf_to_double(&ymin);
                printf("\033[33m%7.2g\033[0m \033[90m+", ymin_d);
                for (x = 0; x < plot_width; x++) printf("-");
                printf("+\033[0m\n");
            }
            
            /* X axis labels - compact */
            {
                double xmin_d = apf_to_double(&xmin);
                double xmax_d = apf_to_double(&xmax);
                double xmid = (xmin_d + xmax_d) / 2;
                int label_space = (plot_width - 21) / 2;  /* Space between labels */
                printf("        \033[33m%-7.2g", xmin_d);
                for (i = 0; i < label_space; i++) printf(" ");
                printf("%7.2g", xmid);
                for (i = 0; i < label_space; i++) printf(" ");
                printf("%7.2g\033[0m\n", xmax_d);
            }
            
            fflush(stdout);
            need_redraw = 0;
        }
        
        /* Read key */
        {
            int key = iplot_read_key();
            apf shift, temp_key;
            
            switch (key) {
                case 'q':
                case 'Q':
                case 27:   /* ESC */
                    running = 0;
                    break;
                    
                case 1004:  /* Left - pan left */
                    apf_sub(&x_range, &xmax, &xmin);
                    apf_from_int(&temp_key, 20);
                    apf_div(&shift, &x_range, &temp_key);
                    apf_sub(&xmin, &xmin, &shift);
                    apf_sub(&xmax, &xmax, &shift);
                    need_redraw = 1;
                    break;
                    
                case 1003:  /* Right - pan right */
                    apf_sub(&x_range, &xmax, &xmin);
                    apf_from_int(&temp_key, 20);
                    apf_div(&shift, &x_range, &temp_key);
                    apf_add(&xmin, &xmin, &shift);
                    apf_add(&xmax, &xmax, &shift);
                    need_redraw = 1;
                    break;
                    
                case 1001:  /* Up - pan up (increase y) */
                    if (auto_y) {
                        /* Switch to fixed y mode based on current view */
                        apf_copy(&ymin_fixed, &ymin);
                        apf_copy(&ymax_fixed, &ymax);
                        auto_y = 0;
                    }
                    apf_sub(&y_range, &ymax_fixed, &ymin_fixed);
                    apf_from_int(&temp_key, 20);
                    apf_div(&shift, &y_range, &temp_key);
                    apf_add(&ymin_fixed, &ymin_fixed, &shift);
                    apf_add(&ymax_fixed, &ymax_fixed, &shift);
                    need_redraw = 1;
                    break;
                    
                case 1002:  /* Down - pan down (decrease y) */
                    if (auto_y) {
                        apf_copy(&ymin_fixed, &ymin);
                        apf_copy(&ymax_fixed, &ymax);
                        auto_y = 0;
                    }
                    apf_sub(&y_range, &ymax_fixed, &ymin_fixed);
                    apf_from_int(&temp_key, 20);
                    apf_div(&shift, &y_range, &temp_key);
                    apf_sub(&ymin_fixed, &ymin_fixed, &shift);
                    apf_sub(&ymax_fixed, &ymax_fixed, &shift);
                    need_redraw = 1;
                    break;
                    
                case '+':
                case '=':  /* Zoom in */
                    apf_sub(&x_range, &xmax, &xmin);
                    apf_from_int(&temp_key, 20);
                    apf_div(&shift, &x_range, &temp_key);
                    apf_add(&xmin, &xmin, &shift);
                    apf_sub(&xmax, &xmax, &shift);
                    if (!auto_y) {
                        apf_sub(&y_range, &ymax_fixed, &ymin_fixed);
                        apf_div(&shift, &y_range, &temp_key);
                        apf_add(&ymin_fixed, &ymin_fixed, &shift);
                        apf_sub(&ymax_fixed, &ymax_fixed, &shift);
                    }
                    need_redraw = 1;
                    break;
                    
                case '-':
                case '_':  /* Zoom out */
                    apf_sub(&x_range, &xmax, &xmin);
                    apf_from_int(&temp_key, 18);
                    apf_div(&shift, &x_range, &temp_key);
                    apf_sub(&xmin, &xmin, &shift);
                    apf_add(&xmax, &xmax, &shift);
                    if (!auto_y) {
                        apf_sub(&y_range, &ymax_fixed, &ymin_fixed);
                        apf_div(&shift, &y_range, &temp_key);
                        apf_sub(&ymin_fixed, &ymin_fixed, &shift);
                        apf_add(&ymax_fixed, &ymax_fixed, &shift);
                    }
                    need_redraw = 1;
                    break;
                    
                case 'a':
                case 'A':  /* Toggle auto-y */
                    auto_y = !auto_y;
                    need_redraw = 1;
                    break;
                    
                case 'r':
                case 'R':  /* Reset */
                    if (apf_eq(xmin_init, xmax_init)) {
                        apf_from_int(&xmin, -10);
                        apf_from_int(&xmax, 10);
                    } else {
                        apf_copy(&xmin, xmin_init);
                        apf_copy(&xmax, xmax_init);
                    }
                    auto_y = 1;
                    need_redraw = 1;
                    break;
            }
        }
    }
    
    /* Cleanup */
    free(y_values);
    free(valid);
    free(screen);
    iplot_restore_term();
    
    return 1;
}

/* ========================================================================
 * Interactive Mandelbrot Viewer
 * ======================================================================== */

/* ANSI color palette for Mandelbrot - 256 color mode */
static const char *mandel_colors[] = {
    "\033[38;5;17m",   /* Deep blue */
    "\033[38;5;18m",
    "\033[38;5;19m",
    "\033[38;5;20m",
    "\033[38;5;21m",
    "\033[38;5;27m",
    "\033[38;5;33m",
    "\033[38;5;39m",
    "\033[38;5;45m",
    "\033[38;5;51m",   /* Cyan */
    "\033[38;5;50m",
    "\033[38;5;49m",
    "\033[38;5;48m",
    "\033[38;5;47m",
    "\033[38;5;46m",   /* Green */
    "\033[38;5;82m",
    "\033[38;5;118m",
    "\033[38;5;154m",
    "\033[38;5;190m",
    "\033[38;5;226m",  /* Yellow */
    "\033[38;5;220m",
    "\033[38;5;214m",
    "\033[38;5;208m",
    "\033[38;5;202m",
    "\033[38;5;196m",  /* Red */
    "\033[38;5;197m",
    "\033[38;5;198m",
    "\033[38;5;199m",
    "\033[38;5;200m",
    "\033[38;5;201m",  /* Magenta */
    "\033[38;5;165m",
    "\033[38;5;129m",
};
#define MANDEL_NUM_COLORS 32

/* Characters for density */
static const char mandel_chars[] = " .:-=+*#%@";
#define MANDEL_NUM_CHARS 10

static int mandelbrot_iter(double cx, double cy, int max_iter)
{
    double zx = 0, zy = 0;
    double zx2, zy2;
    int iter;
    
    for (iter = 0; iter < max_iter; iter++) {
        zx2 = zx * zx;
        zy2 = zy * zy;
        if (zx2 + zy2 > 4.0) return iter;
        zy = 2 * zx * zy + cy;
        zx = zx2 - zy2 + cx;
    }
    return max_iter;
}

int do_mandelbrot(void)
{
    int width, height;
    int plot_width, plot_height;
    double xmin = -2.5, xmax = 1.0, ymin, ymax;
    int max_iter = 100;
    int running = 1;
    int need_redraw = 1;
    int x, y;
    int *iter_buf;
    int color_mode = 1;  /* 0 = ASCII, 1 = color */
    double aspect_ratio = 0.5;  /* Chars are ~2x tall as wide */
    
    /* Setup terminal */
    iplot_setup_term();
    iplot_get_terminal_size(&width, &height);
    
    plot_width = width;
    plot_height = height - 2;  /* Leave room for status line */
    
    /* Calculate y range to maintain aspect ratio */
    {
        double x_range = xmax - xmin;
        double y_range = x_range * plot_height / (plot_width * aspect_ratio);
        double y_center = 0.0;  /* Mandelbrot is roughly centered at y=0 */
        ymin = y_center - y_range / 2;
        ymax = y_center + y_range / 2;
    }
    
    if (plot_width < 40 || plot_height < 20) {
        iplot_restore_term();
        printf("Terminal too small for mandelbrot\n");
        return 0;
    }
    
    iter_buf = (int *)malloc(plot_width * plot_height * sizeof(int));
    if (!iter_buf) {
        iplot_restore_term();
        printf("Out of memory\n");
        return 0;
    }
    
    while (running) {
        if (need_redraw) {
            double dx = (xmax - xmin) / plot_width;
            double dy = (ymax - ymin) / plot_height;
            int max_found = 0;
            
            /* Calculate all iterations */
            for (y = 0; y < plot_height; y++) {
                double cy = ymax - y * dy;
                for (x = 0; x < plot_width; x++) {
                    double cx = xmin + x * dx;
                    int iter = mandelbrot_iter(cx, cy, max_iter);
                    iter_buf[y * plot_width + x] = iter;
                    if (iter < max_iter && iter > max_found) max_found = iter;
                }
            }
            
            /* Render */
            printf("\033[H");  /* Home cursor */
            
            for (y = 0; y < plot_height; y++) {
                for (x = 0; x < plot_width; x++) {
                    int iter = iter_buf[y * plot_width + x];
                    
                    if (iter == max_iter) {
                        /* In set - black */
                        if (color_mode) {
                            printf("\033[38;5;0m ");
                        } else {
                            putchar(' ');
                        }
                    } else {
                        if (color_mode) {
                            int ci = (iter * MANDEL_NUM_COLORS / (max_found + 1)) % MANDEL_NUM_COLORS;
                            printf("%s#", mandel_colors[ci]);
                        } else {
                            int ci = (iter * MANDEL_NUM_CHARS / (max_found + 1)) % MANDEL_NUM_CHARS;
                            putchar(mandel_chars[ci]);
                        }
                    }
                }
                if (y < plot_height - 1) putchar('\n');
            }
            
            /* Status line */
            printf("\033[0m\n");
            printf("\033[K\033[1;36mMandelbrot\033[0m x:[%.6f,%.6f] y:[%.6f,%.6f] iter:%d  ", 
                   xmin, xmax, ymin, ymax, max_iter);
            printf("\033[90m[arrows=pan +/-=zoom i/I=iter c=color r=reset q=quit]\033[0m");
            fflush(stdout);
            
            need_redraw = 0;
        }
        
        /* Read key */
        {
            int key = iplot_read_key();
            double xrange = xmax - xmin;
            double yrange = ymax - ymin;
            double shift;
            
            switch (key) {
                case 'q':
                case 'Q':
                case 27:
                    running = 0;
                    break;
                    
                case 1004:  /* Left */
                    shift = xrange * 0.1;
                    xmin -= shift;
                    xmax -= shift;
                    need_redraw = 1;
                    break;
                    
                case 1003:  /* Right */
                    shift = xrange * 0.1;
                    xmin += shift;
                    xmax += shift;
                    need_redraw = 1;
                    break;
                    
                case 1001:  /* Up */
                    shift = yrange * 0.1;
                    ymin += shift;
                    ymax += shift;
                    need_redraw = 1;
                    break;
                    
                case 1002:  /* Down */
                    shift = yrange * 0.1;
                    ymin -= shift;
                    ymax -= shift;
                    need_redraw = 1;
                    break;
                    
                case '+':
                case '=':  /* Zoom in */
                    {
                        double cx = (xmin + xmax) / 2;
                        double cy = (ymin + ymax) / 2;
                        xrange *= 0.8;
                        yrange *= 0.8;
                        xmin = cx - xrange / 2;
                        xmax = cx + xrange / 2;
                        ymin = cy - yrange / 2;
                        ymax = cy + yrange / 2;
                        /* Increase iterations when zooming in */
                        if (max_iter < 1000) max_iter = (int)(max_iter * 1.1);
                    }
                    need_redraw = 1;
                    break;
                    
                case '-':
                case '_':  /* Zoom out */
                    {
                        double cx = (xmin + xmax) / 2;
                        double cy = (ymin + ymax) / 2;
                        xrange *= 1.25;
                        yrange *= 1.25;
                        xmin = cx - xrange / 2;
                        xmax = cx + xrange / 2;
                        ymin = cy - yrange / 2;
                        ymax = cy + yrange / 2;
                        if (max_iter > 50) max_iter = (int)(max_iter * 0.9);
                    }
                    need_redraw = 1;
                    break;
                    
                case 'i':  /* Decrease iterations */
                    if (max_iter > 20) max_iter -= 10;
                    need_redraw = 1;
                    break;
                    
                case 'I':  /* Increase iterations */
                    if (max_iter < 2000) max_iter += 10;
                    need_redraw = 1;
                    break;
                    
                case 'c':
                case 'C':  /* Toggle color mode */
                    color_mode = !color_mode;
                    need_redraw = 1;
                    break;
                    
                case 'r':
                case 'R':  /* Reset */
                    xmin = -2.5; xmax = 1.0;
                    {
                        double x_range = xmax - xmin;
                        double y_range = x_range * plot_height / (plot_width * aspect_ratio);
                        ymin = -y_range / 2;
                        ymax = y_range / 2;
                    }
                    max_iter = 100;
                    need_redraw = 1;
                    break;
                    
                case 'w':  /* Zoom to specific area - cardioid */
                    {
                        double cx = -0.75, cy = 0.15;
                        double x_range = 0.1;
                        double y_range = x_range * plot_height / (plot_width * aspect_ratio);
                        xmin = cx - x_range / 2; xmax = cx + x_range / 2;
                        ymin = cy - y_range / 2; ymax = cy + y_range / 2;
                    }
                    max_iter = 200;
                    need_redraw = 1;
                    break;
                    
                case 'e':  /* Elephant valley */
                    {
                        double cx = 0.27, cy = 0.01;
                        double x_range = 0.02;
                        double y_range = x_range * plot_height / (plot_width * aspect_ratio);
                        xmin = cx - x_range / 2; xmax = cx + x_range / 2;
                        ymin = cy - y_range / 2; ymax = cy + y_range / 2;
                    }
                    max_iter = 300;
                    need_redraw = 1;
                    break;
                    
                case 's':  /* Seahorse valley */
                    {
                        double cx = -0.745, cy = 0.105;
                        double x_range = 0.01;
                        double y_range = x_range * plot_height / (plot_width * aspect_ratio);
                        xmin = cx - x_range / 2; xmax = cx + x_range / 2;
                        ymin = cy - y_range / 2; ymax = cy + y_range / 2;
                    }
                    max_iter = 300;
                    need_redraw = 1;
                    break;
            }
        }
    }
    
    free(iter_buf);
    iplot_restore_term();
    
    return 1;
}

/* Mandelbrot iteration function callable from calculator */
void fn_mandelbrot_iter(apfc *result, const apfc *cx, const apfc *cy, const apfc *max_iter_apf)
{
    double x = apf_to_double(&cx->re);
    double y = apf_to_double(&cy->re);
    int max_iter = (int)apf_to_double(&max_iter_apf->re);
    int iter;
    
    if (max_iter < 1) max_iter = 100;
    if (max_iter > 10000) max_iter = 10000;
    
    iter = mandelbrot_iter(x, y, max_iter);
    
    apf_from_int(&result->re, iter);
    apf_zero(&result->im);
}

#elif defined(_WIN32)
/* Windows stubs */
int do_iplot(const char *expr, char var, apf *xmin_init, apf *xmax_init)
{
    (void)expr; (void)var; (void)xmin_init; (void)xmax_init;
    printf("iplot not implemented on Windows yet - use plot instead\n");
    return 0;
}

int do_mandelbrot(void)
{
    printf("mandelbrot viewer not implemented on Windows yet\n");
    return 0;
}

void fn_mandelbrot_iter(apfc *result, const apfc *cx, const apfc *cy, const apfc *max_iter_apf)
{
    (void)cx; (void)cy; (void)max_iter_apf;
    apf_zero(&result->re);
    apf_zero(&result->im);
}

#else
/* Stubs for DOS/other */
int do_iplot(const char *expr, char var, apf *xmin_init, apf *xmax_init)
{
    (void)expr; (void)var; (void)xmin_init; (void)xmax_init;
    printf("iplot requires Unix terminal\n");
    return 0;
}

int do_mandelbrot(void)
{
    printf("mandelbrot viewer requires Unix terminal\n");
    return 0;
}

void fn_mandelbrot_iter(apfc *result, const apfc *cx, const apfc *cy, const apfc *max_iter_apf)
{
    (void)cx; (void)cy; (void)max_iter_apf;
    apf_zero(&result->re);
    apf_zero(&result->im);
}
#endif

/* ========================================================================
 * Planet Ephemeris - Simple circular orbit model
 * ======================================================================== */

#include <math.h>
#include <strings.h>  /* For strcasecmp */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Comprehensive planet data structure */
typedef struct {
    const char *name;
    double a;           /* Semi-major axis in AU */
    double e;           /* Eccentricity */
    double i;           /* Inclination to ecliptic (degrees) */
    double omega;       /* Longitude of ascending node (degrees) */
    double w;           /* Argument of perihelion (degrees) */
    double L0;          /* Mean longitude at epoch J2000 (degrees) */
    double period;      /* Orbital period in days */
    double radius;      /* Physical radius in AU */
    double mass;        /* Mass in Earth masses */
    double rot_period;  /* Rotation period in hours (negative = retrograde) */
} planet_data_t;

/* Accurate orbital elements (J2000 epoch) */
/* Sources: NASA JPL, IAU */
static const planet_data_t planet_db[] = {
    /* Sun - stationary at origin */
    {"sun",     0.0,      0.0,     0.0,    0.0,    0.0,    0.0,
                1.0,      0.00465047, 332946.0, 609.12},
    
    /* Mercury */
    {"mercury", 0.38710,  0.20563, 7.005,  48.331, 29.124, 252.251,
                87.969,   0.0000163,  0.0553,   1407.6},
    
    /* Venus */
    {"venus",   0.72333,  0.00677, 3.395,  76.680, 54.884, 181.980,
                224.701,  0.0000404,  0.815,    -5832.5},  /* Retrograde */
    
    /* Earth */
    {"earth",   1.00000,  0.01671, 0.000,  -11.26, 102.95, 100.464,
                365.256,  0.0000426,  1.000,    23.934},
    
    /* Mars */
    {"mars",    1.52368,  0.09341, 1.850,  49.558, 286.50, 355.453,
                686.980,  0.0000227,  0.107,    24.623},
    
    /* Jupiter */
    {"jupiter", 5.20260,  0.04839, 1.303,  100.46, 273.87, 34.404,
                4332.59,  0.000467,   317.83,   9.925},
    
    /* Saturn */
    {"saturn",  9.55491,  0.05415, 2.489,  113.67, 339.39, 50.077,
                10759.22, 0.000389,   95.16,    10.656},
    
    /* Uranus */
    {"uranus",  19.2184,  0.04717, 0.773,  74.006, 96.999, 314.055,
                30685.4,  0.000170,   14.54,    -17.24},  /* Retrograde */
    
    /* Neptune */
    {"neptune", 30.1104,  0.00859, 1.770,  131.78, 276.34, 304.349,
                60189.0,  0.000165,   17.15,    16.11},
    
    /* Moon - orbits Earth, inclination to ecliptic ~5.145 degrees */
    {"moon",    0.00257,  0.0549,  5.145,  125.08, 318.15, 218.32,
                27.32166, 0.0000116,  0.0123,   655.73},
    
    {NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};

static const planet_data_t *find_planet_data(const char *name)
{
    int i;
    for (i = 0; planet_db[i].name; i++) {
        if (strcasecmp(planet_db[i].name, name) == 0) {
            return &planet_db[i];
        }
    }
    return NULL;
}

/* Calculate planet position with full Keplerian elements */
/* Returns x, y (ecliptic plane) and z (height above ecliptic) */
static void planet_pos_3d(const char *name, double t_days, 
                          double *x, double *y, double *z)
{
    const planet_data_t *p = find_planet_data(name);
    double M, E, E_prev, nu, r;
    double cos_nu, sin_nu;
    double x_orb, y_orb;
    double cos_w, sin_w, cos_O, sin_O, cos_i, sin_i;
    int iter;
    
    if (!p || p->a == 0) {
        *x = *y = *z = 0;
        return;
    }
    
    /* Mean anomaly */
    M = (p->L0 + 360.0 * t_days / p->period) * M_PI / 180.0;
    M = fmod(M, 2 * M_PI);
    if (M < 0) M += 2 * M_PI;
    
    /* Solve Kepler's equation: E - e*sin(E) = M */
    E = M;  /* Initial guess */
    for (iter = 0; iter < 10; iter++) {
        E_prev = E;
        E = M + p->e * sin(E);
        if (fabs(E - E_prev) < 1e-8) break;
    }
    
    /* True anomaly */
    cos_nu = (cos(E) - p->e) / (1 - p->e * cos(E));
    sin_nu = sqrt(1 - p->e * p->e) * sin(E) / (1 - p->e * cos(E));
    nu = atan2(sin_nu, cos_nu);
    
    /* Heliocentric distance */
    r = p->a * (1 - p->e * cos(E));
    
    /* Position in orbital plane */
    x_orb = r * cos(nu);
    y_orb = r * sin(nu);
    
    /* Convert to ecliptic coordinates */
    cos_w = cos(p->w * M_PI / 180.0);
    sin_w = sin(p->w * M_PI / 180.0);
    cos_O = cos(p->omega * M_PI / 180.0);
    sin_O = sin(p->omega * M_PI / 180.0);
    cos_i = cos(p->i * M_PI / 180.0);
    sin_i = sin(p->i * M_PI / 180.0);
    
    /* Rotation matrices combined */
    *x = (cos_w * cos_O - sin_w * sin_O * cos_i) * x_orb +
         (-sin_w * cos_O - cos_w * sin_O * cos_i) * y_orb;
    *y = (cos_w * sin_O + sin_w * cos_O * cos_i) * x_orb +
         (-sin_w * sin_O + cos_w * cos_O * cos_i) * y_orb;
    *z = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb;
    
    /* Special case: Moon orbits Earth */
    if (strcmp(name, "moon") == 0) {
        double ex, ey, ez;
        planet_pos_3d("earth", t_days, &ex, &ey, &ez);
        *x += ex;
        *y += ey;
        *z += ez;
    }
}

/* Legacy 2D wrapper */
static void planet_pos(const char *name, double t_days, double *x, double *y)
{
    double z;
    planet_pos_3d(name, t_days, x, y, &z);
}

#if !defined(HAVE_CONIO) && !defined(_WIN32)

/* ========================================================================
 * Interactive Scatter Plot (iscatter)
 * ======================================================================== */

/* Density characters for ASCII mode */
static const char density_chars[] = " .'`^\",:;Il!i><~+_-?][}{1)(|/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
#define NUM_DENSITY_CHARS 70

int do_iscatter(double *xs, double *ys, int *iters, int n, int max_iter)
{
    int width, height;
    int plot_width, plot_height;
    double xmin, xmax, ymin, ymax;
    int running = 1;
    int need_redraw = 1;
    int x, y, i;
    int *screen_iter;  /* Iteration count per cell */
    int *screen_count; /* Point count per cell */
    int color_mode = 1;
    
    /* Find bounds */
    xmin = xmax = xs[0];
    ymin = ymax = ys[0];
    for (i = 1; i < n; i++) {
        if (xs[i] < xmin) xmin = xs[i];
        if (xs[i] > xmax) xmax = xs[i];
        if (ys[i] < ymin) ymin = ys[i];
        if (ys[i] > ymax) ymax = ys[i];
    }
    
    /* Add margin */
    {
        double mx = (xmax - xmin) * 0.05;
        double my = (ymax - ymin) * 0.05;
        if (mx < 0.001) mx = 0.1;
        if (my < 0.001) my = 0.1;
        xmin -= mx; xmax += mx;
        ymin -= my; ymax += my;
    }
    
    iplot_setup_term();
    iplot_get_terminal_size(&width, &height);
    
    plot_width = width;
    plot_height = height - 2;
    
    screen_iter = (int *)malloc(plot_width * plot_height * sizeof(int));
    screen_count = (int *)malloc(plot_width * plot_height * sizeof(int));
    
    if (!screen_iter || !screen_count) {
        iplot_restore_term();
        free(screen_iter); free(screen_count);
        return 0;
    }
    
    while (running) {
        if (need_redraw) {
            double dx = (xmax - xmin) / plot_width;
            double dy = (ymax - ymin) / plot_height;
            
            memset(screen_iter, 0, plot_width * plot_height * sizeof(int));
            memset(screen_count, 0, plot_width * plot_height * sizeof(int));
            
            /* Bin points */
            for (i = 0; i < n; i++) {
                int px = (int)((xs[i] - xmin) / dx);
                int py = (int)((ymax - ys[i]) / dy);
                if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                    int idx = py * plot_width + px;
                    screen_iter[idx] += iters[i];
                    screen_count[idx]++;
                }
            }
            
            /* Render */
            printf("\033[H");
            
            for (y = 0; y < plot_height; y++) {
                for (x = 0; x < plot_width; x++) {
                    int idx = y * plot_width + x;
                    int count = screen_count[idx];
                    
                    if (count == 0) {
                        putchar(' ');
                    } else {
                        int avg_iter = screen_iter[idx] / count;
                        
                        if (color_mode) {
                            int ci = (avg_iter * MANDEL_NUM_COLORS / max_iter) % MANDEL_NUM_COLORS;
                            if (avg_iter >= max_iter) {
                                printf("\033[38;5;0m ");
                            } else {
                                printf("%s#", mandel_colors[ci]);
                            }
                        } else {
                            int di = (avg_iter * NUM_DENSITY_CHARS / max_iter);
                            if (di >= NUM_DENSITY_CHARS) di = NUM_DENSITY_CHARS - 1;
                            putchar(density_chars[di]);
                        }
                    }
                }
                if (y < plot_height - 1) putchar('\n');
            }
            
            printf("\033[0m\n");
            printf("\033[K\033[1;36miscatter\033[0m n=%d x:[%.4g,%.4g] y:[%.4g,%.4g]  ",
                   n, xmin, xmax, ymin, ymax);
            printf("\033[90m[arrows=pan +/-=zoom c=color q=quit]\033[0m");
            fflush(stdout);
            
            need_redraw = 0;
        }
        
        {
            int key = iplot_read_key();
            double xrange = xmax - xmin;
            double yrange = ymax - ymin;
            double shift;
            
            switch (key) {
                case 'q': case 'Q': case 27:
                    running = 0;
                    break;
                case 1004:  /* Left */
                    shift = xrange * 0.1;
                    xmin -= shift; xmax -= shift;
                    need_redraw = 1;
                    break;
                case 1003:  /* Right */
                    shift = xrange * 0.1;
                    xmin += shift; xmax += shift;
                    need_redraw = 1;
                    break;
                case 1001:  /* Up */
                    shift = yrange * 0.1;
                    ymin += shift; ymax += shift;
                    need_redraw = 1;
                    break;
                case 1002:  /* Down */
                    shift = yrange * 0.1;
                    ymin -= shift; ymax -= shift;
                    need_redraw = 1;
                    break;
                case '+': case '=':
                    {
                        double cx = (xmin + xmax) / 2;
                        double cy = (ymin + ymax) / 2;
                        xrange *= 0.8; yrange *= 0.8;
                        xmin = cx - xrange/2; xmax = cx + xrange/2;
                        ymin = cy - yrange/2; ymax = cy + yrange/2;
                    }
                    need_redraw = 1;
                    break;
                case '-': case '_':
                    {
                        double cx = (xmin + xmax) / 2;
                        double cy = (ymin + ymax) / 2;
                        xrange *= 1.25; yrange *= 1.25;
                        xmin = cx - xrange/2; xmax = cx + xrange/2;
                        ymin = cy - yrange/2; ymax = cy + yrange/2;
                    }
                    need_redraw = 1;
                    break;
                case 'c': case 'C':
                    color_mode = !color_mode;
                    need_redraw = 1;
                    break;
            }
        }
    }
    
    free(screen_iter);
    free(screen_count);
    iplot_restore_term();
    return 1;
}

/* ========================================================================
 * Interactive Time Scatter - Solar System Viewer (tscatter)
 * ======================================================================== */

#include <sys/select.h>
#include <unistd.h>

/* Planet display info - just visual properties */
typedef struct {
    const char *name;
    char symbol;
    const char *color;
} planet_display_t;

static const planet_display_t planet_displays[] = {
    {"sun",     'O', "\033[1;33m"},     /* Yellow */
    {"mercury", '.', "\033[38;5;250m"}, /* Gray */
    {"venus",   'v', "\033[38;5;229m"}, /* Pale yellow */
    {"earth",   'e', "\033[1;34m"},     /* Blue */
    {"mars",    'm', "\033[1;31m"},     /* Red */
    {"jupiter", 'J', "\033[38;5;208m"}, /* Orange */
    {"saturn",  'S', "\033[38;5;222m"}, /* Gold */
    {"uranus",  'U', "\033[1;36m"},     /* Cyan */
    {"neptune", 'N', "\033[1;34m"},     /* Blue */
    {"moon",    'o', "\033[38;5;252m"}, /* Light gray */
    {NULL, 0, NULL}
};

static const planet_display_t *get_planet_display(const char *name)
{
    int i;
    for (i = 0; planet_displays[i].name; i++) {
        if (strcasecmp(planet_displays[i].name, name) == 0) {
            return &planet_displays[i];
        }
    }
    return NULL;
}

/* Get planet radius from planet_db */
static double get_planet_radius(const char *name)
{
    const planet_data_t *p = find_planet_data(name);
    return p ? p->radius : 0.00001;
}

/* Draw a filled circle using Bresenham-style algorithm
 * Scales from single dot to full screen
 * radius_x: radius in character cells (horizontal)
 * Uses 'O' for interior, '.' for edge pixels when large enough
 */
static void draw_filled_circle(char *screen, char *colors, int screen_w, int screen_h,
                               int cx, int cy, double radius_x, char symbol, int color_idx)
{
    int rx, ry;
    int dy, dx;
    
    /* Radius in x (chars) and y (account for char aspect ~2:1) */
    rx = (int)(radius_x + 0.5);
    ry = (int)(radius_x / 2.0 + 0.5);  /* Chars are ~2x tall as wide */
    
    /* Minimum size: single character */
    if (rx < 1) rx = 1;
    if (ry < 1) ry = 1;
    
    /* For tiny circles, just draw a single char */
    if (rx == 1 && ry == 1) {
        if (cx >= 0 && cx < screen_w && cy >= 0 && cy < screen_h) {
            int idx = cy * screen_w + cx;
            /* Size determines character: . for tiny, o for small, O for medium+ */
            char c = (radius_x < 0.3) ? '.' : ((radius_x < 0.7) ? 'o' : symbol);
            screen[idx] = c;
            colors[idx] = color_idx;
        }
        return;
    }
    
    /* Draw filled ellipse */
    for (dy = -ry; dy <= ry; dy++) {
        for (dx = -rx; dx <= rx; dx++) {
            /* Ellipse equation: (dx/rx)^2 + (dy/ry)^2 <= 1 */
            double ex = (double)dx / rx;
            double ey = (double)dy / ry;
            double dist_sq = ex*ex + ey*ey;
            
            if (dist_sq <= 1.0) {
                int sx = cx + dx;
                int sy = cy + dy;
                if (sx >= 0 && sx < screen_w && sy >= 0 && sy < screen_h) {
                    int idx = sy * screen_w + sx;
                    char c;
                    
                    /* Edge detection: use '.' for pixels near edge, 'O' for interior */
                    if (rx >= 3 && dist_sq > 0.7) {
                        c = '.';  /* Edge */
                    } else if (rx >= 2 && dist_sq > 0.85) {
                        c = 'o';  /* Near edge */
                    } else {
                        c = symbol;  /* Interior - use planet symbol */
                    }
                    
                    screen[idx] = c;
                    colors[idx] = color_idx;
                }
            }
        }
    }
}

int do_tscatter(const char **bodies, int num_bodies, const char *follow_init)
{
    int width, height;
    int plot_width, plot_height;
    double t_days;          /* Days since J2000 */
    double dt = 1.0;        /* Days per step */
    int running = 1;
    int need_redraw = 1;
    int playing = 0;
    int x, y, i;
    time_t now;
    char *screen;
    char *colors;          /* Color index per cell */
    double xmin, xmax, ymin, ymax;
    double view_scale = 2.0;  /* AU - initial view radius */
    double center_x = 0, center_y = 0;
    double aspect_ratio;   /* Character aspect correction */
    int view_3d = 0;       /* 0 = top-down (X-Y), 1 = tilted 3D view */
    double tilt_angle = 30.0;  /* Tilt angle in degrees (0=top-down, 90=edge-on) */
    double size_scale = 1.0;   /* Planet size exaggeration factor */
    
    /* Initialize to current time: days since J2000 (2000-01-01 12:00:00 UTC) */
    /* Follow mode: 0=sun, 1=earth, 2=earth-moon barycenter, 3+=cycle through bodies */
    int follow_mode = 0;
    
    /* Initialize to current time: days since J2000 (2000-01-01 12:00:00 UTC) */
    now = time(NULL);
    t_days = (double)(now - 946728000) / 86400.0;
    
    /* Set initial follow mode based on parameter */
    if (follow_init) {
        if (strcasecmp(follow_init, "earth") == 0) follow_mode = 1;
        else if (strcasecmp(follow_init, "earth,moon") == 0 || 
                 strcasecmp(follow_init, "earth+moon") == 0 ||
                 strcasecmp(follow_init, "moon") == 0) {
            follow_mode = 2;
            view_scale = 0.006;  /* Zoom in to see moon orbit (0.00257 AU) */
        }
    }
    
    iplot_setup_term();;
    iplot_get_terminal_size(&width, &height);
    
    plot_width = width;
    plot_height = height - 3;
    
    /* Characters are ~2x tall as wide, so scale x by 0.5 for square aspect */
    aspect_ratio = 0.5;
    
    screen = (char *)malloc(plot_width * plot_height);
    colors = (char *)malloc(plot_width * plot_height);
    
    if (!screen || !colors) {
        iplot_restore_term();
        free(screen); free(colors);
        return 0;
    }
    
    while (running) {
        if (need_redraw) {
            double positions[20][3];  /* x, y, z for each body */
            double follow_x = 0, follow_y = 0, follow_z = 0;
            double x_scale, y_scale;
            const char *follow_label = "sun";
            char view_label[32];
            double cos_tilt, sin_tilt;
            
            /* Calculate tilt transform */
            cos_tilt = cos(tilt_angle * M_PI / 180.0);
            sin_tilt = sin(tilt_angle * M_PI / 180.0);
            
            if (view_3d) {
                sprintf(view_label, "3D %.0f°", tilt_angle);
            } else {
                strcpy(view_label, "2D");
            }
            
            /* Get 3D positions */
            for (i = 0; i < num_bodies && i < 20; i++) {
                planet_pos_3d(bodies[i], t_days, &positions[i][0], &positions[i][1], &positions[i][2]);
            }
            
            /* Calculate follow center based on follow_mode */
            if (follow_mode == 0) {
                /* Sun - origin */
                follow_x = follow_y = follow_z = 0;
                follow_label = "sun";
            } else if (follow_mode == 1) {
                /* Earth */
                planet_pos_3d("earth", t_days, &follow_x, &follow_y, &follow_z);
                follow_label = "earth";
            } else if (follow_mode == 2) {
                /* Earth-Moon barycenter (mass-weighted) */
                /* Earth mass = 1.0, Moon mass = 0.0123 Earth masses */
                double ex, ey, ez, mx, my, mz;
                double earth_mass = 1.0, moon_mass = 0.0123;
                double total_mass = earth_mass + moon_mass;
                planet_pos_3d("earth", t_days, &ex, &ey, &ez);
                planet_pos_3d("moon", t_days, &mx, &my, &mz);
                follow_x = (ex * earth_mass + mx * moon_mass) / total_mass;
                follow_y = (ey * earth_mass + my * moon_mass) / total_mass;
                follow_z = (ez * earth_mass + mz * moon_mass) / total_mass;
                follow_label = "E-M bary";
            } else if (follow_mode >= 3 && follow_mode - 3 < num_bodies) {
                /* Follow specific body from list */
                int bi = follow_mode - 3;
                follow_x = positions[bi][0];
                follow_y = positions[bi][1];
                follow_z = positions[bi][2];
                follow_label = bodies[bi];
            }
            
            /* Project follow center based on view mode */
            if (view_3d) {
                /* 3D: Y on screen = Y*cos(tilt) + Z*sin(tilt) */
                center_x = follow_x;
                center_y = follow_y * cos_tilt + follow_z * sin_tilt;
            } else {
                center_x = follow_x;
                center_y = follow_y;
            }
            
            /* Set view bounds with aspect ratio correction */
            x_scale = view_scale * (plot_width * aspect_ratio) / plot_height;
            y_scale = view_scale;
            
            xmin = center_x - x_scale;
            xmax = center_x + x_scale;
            ymin = center_y - y_scale;
            ymax = center_y + y_scale;
            
            /* Clear screen buffer */
            memset(screen, ' ', plot_width * plot_height);
            memset(colors, 0, plot_width * plot_height);
            
            /* Draw orbit traces (faint) - only in 2D mode */
            if (!view_3d) {
                for (i = 0; i < num_bodies; i++) {
                    const planet_data_t *pd = find_planet_data(bodies[i]);
                    if (pd && pd->a > 0 && strcmp(bodies[i], "moon") != 0) {
                        /* Draw elliptical orbit using orbital elements */
                        int steps = 100;
                        int s;
                        for (s = 0; s < steps; s++) {
                            double theta = 2 * M_PI * s / steps;
                            /* Simplified ellipse in orbital plane */
                            double r = pd->a * (1 - pd->e * pd->e) / (1 + pd->e * cos(theta));
                            double ox = r * cos(theta) - follow_x;
                            double oy = r * sin(theta) - follow_y;
                            
                            int px = (int)((ox - (xmin - follow_x)) / (xmax - xmin) * plot_width);
                            int py = (int)(((ymax - follow_y) - oy) / (ymax - ymin) * plot_height);
                            
                            if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                                int idx = py * plot_width + px;
                                if (screen[idx] == ' ') {
                                    screen[idx] = '.';
                                    colors[idx] = 10;  /* Dim */
                                }
                            }
                        }
                    }
                }
            } else {
                /* In 3D mode, draw a line at z=0 (ecliptic) */
                int px;
                int py = (int)((ymax - 0) / (ymax - ymin) * plot_height);
                for (px = 0; px < plot_width; px++) {
                    /* Draw ecliptic plane line */
                    if (py >= 0 && py < plot_height) {
                        int idx = py * plot_width + px;
                        if (screen[idx] == ' ') {
                            screen[idx] = '-';
                            colors[idx] = 10;
                        }
                    }
                }
            }
            
            /* Draw bodies - smallest first so larger ones overdraw */
            for (i = num_bodies - 1; i >= 0; i--) {
                double bx, by, bz, screen_y;
                double radius_au, x_per_pixel, radius_chars;
                const planet_display_t *disp;
                int cx, cy;
                
                bx = positions[i][0];
                by = positions[i][1];
                bz = positions[i][2];
                disp = get_planet_display(bodies[i]);
                radius_au = get_planet_radius(bodies[i]);
                
                /* Apply size exaggeration */
                radius_au *= size_scale;
                
                /* Project position based on view mode */
                if (view_3d) {
                    /* 3D projection: rotate around X axis by tilt angle */
                    screen_y = by * cos_tilt + bz * sin_tilt;
                } else {
                    screen_y = by;
                }
                
                /* Convert radius to character cells */
                x_per_pixel = (xmax - xmin) / plot_width;
                radius_chars = radius_au / x_per_pixel;
                
                /* Screen coordinates */
                cx = (int)((bx - xmin) / (xmax - xmin) * plot_width);
                cy = (int)((ymax - screen_y) / (ymax - ymin) * plot_height);
                
                /* Draw using filled circle function */
                draw_filled_circle(screen, colors, plot_width, plot_height,
                                   cx, cy, radius_chars, 
                                   disp ? disp->symbol : '*', i + 1);
            }
            
            /* Render */
            printf("\033[H");
            
            for (y = 0; y < plot_height; y++) {
                for (x = 0; x < plot_width; x++) {
                    int idx = y * plot_width + x;
                    char c = screen[idx];
                    int ci = colors[idx];
                    
                    if (ci == 0 || ci == 10) {
                        /* Empty or orbit trace */
                        if (c == '.') {
                            printf("\033[38;5;240m.\033[0m");
                        } else {
                            putchar(' ');
                        }
                    } else {
                        /* Planet */
                        const planet_display_t *disp = get_planet_display(bodies[ci - 1]);
                        if (disp && disp->color) {
                            printf("%s%c\033[0m", disp->color, c);
                        } else {
                            putchar(c);
                        }
                    }
                }
                if (y < plot_height - 1) putchar('\n');
            }
            
            /* Legend */
            printf("\033[0m\n\033[K");
            for (i = 0; i < num_bodies && i < 10; i++) {
                const planet_display_t *disp = get_planet_display(bodies[i]);
                if (disp && disp->color) {
                    printf("%s%c\033[0m=%s ", disp->color, disp->symbol, bodies[i]);
                }
            }
            
            /* Status - show date instead of days */
            {
                time_t display_time = 946728000 + (time_t)(t_days * 86400);
                struct tm *tm_info = gmtime(&display_time);
                char date_buf[32];
                char size_buf[32];
                strftime(date_buf, sizeof(date_buf), "%Y-%m-%d", tm_info);
                if (size_scale >= 1e9) {
                    sprintf(size_buf, "%.0fB×", size_scale / 1e9);
                } else if (size_scale >= 1e6) {
                    sprintf(size_buf, "%.0fM×", size_scale / 1e6);
                } else if (size_scale >= 1e3) {
                    sprintf(size_buf, "%.0fK×", size_scale / 1e3);
                } else if (size_scale > 1.0) {
                    sprintf(size_buf, "%.0f×", size_scale);
                } else {
                    size_buf[0] = '\0';  /* Don't show 1× */
                }
                printf("\n\033[K\033[1;36msolar\033[0m \033[33m%s\033[0m dt=%.1fd %s @\033[36m%s\033[0m %s%s%s  ",
                       date_buf, dt, playing ? "\033[32m▶\033[0m" : "\033[90m⏸\033[0m", 
                       follow_label, 
                       size_buf[0] ? "\033[35m" : "", size_buf, size_buf[0] ? "\033[0m " : "");
                printf("\033[90m%s\033[0m", view_label);
            }
            if (view_3d) {
                printf(" \033[90m[←→=time +/-=zoom s/S=size w/e=tilt d=2D]\033[0m");
            } else {
                printf(" \033[90m[←→=time ↑↓=dt +/-=zoom s/S=size f=follow d=3D]\033[0m");
            }
            fflush(stdout);
            
            need_redraw = 0;
        }
        
        /* Check for auto-advance if playing */
        if (playing) {
            struct timespec ts = {0, 50000000};  /* 50ms = 20 fps */
            nanosleep(&ts, NULL);
            t_days += dt;
            need_redraw = 1;
        }
        
        /* Non-blocking key check when playing */
        {
            fd_set fds;
            struct timeval tv;
            int key = -1;
            
            FD_ZERO(&fds);
            FD_SET(STDIN_FILENO, &fds);
            tv.tv_sec = 0;
            tv.tv_usec = playing ? 0 : 100000;  /* Block if not playing */
            
            if (select(STDIN_FILENO + 1, &fds, NULL, NULL, &tv) > 0) {
                key = iplot_read_key();
            }
            
            if (key > 0) {
                switch (key) {
                    case 'q': case 'Q': case 27:
                        running = 0;
                        break;
                    case 1004:  /* Left - go back in time */
                        t_days -= dt * 10;
                        need_redraw = 1;
                        break;
                    case 1003:  /* Right - go forward */
                        t_days += dt * 10;
                        need_redraw = 1;
                        break;
                    case 1001:  /* Up - increase dt */
                        dt *= 2;
                        if (dt > 365) dt = 365;
                        need_redraw = 1;
                        break;
                    case 1002:  /* Down - decrease dt */
                        dt /= 2;
                        if (dt < 0.01) dt = 0.01;
                        need_redraw = 1;
                        break;
                    case '+': case '=':  /* Zoom in */
                        view_scale *= 0.8;
                        if (view_scale < 0.001) view_scale = 0.001;
                        need_redraw = 1;
                        break;
                    case '-': case '_':  /* Zoom out */
                        view_scale *= 1.25;
                        if (view_scale > 50) view_scale = 50;
                        need_redraw = 1;
                        break;
                    case ' ':  /* Play/pause */
                        playing = !playing;
                        need_redraw = 1;
                        break;
                    case 'r': case 'R':  /* Reset to NOW */
                        now = time(NULL);
                        t_days = (double)(now - 946728000) / 86400.0;
                        dt = 1.0;
                        view_scale = 2.0;
                        follow_mode = 0;
                        need_redraw = 1;
                        break;
                    case 'f':  /* Cycle follow mode forward */
                        follow_mode++;
                        /* Modes: 0=sun, 1=earth, 2=earth+moon, 3+=bodies */
                        if (follow_mode > 2 + num_bodies) follow_mode = 0;
                        need_redraw = 1;
                        break;
                    case 'F':  /* Cycle follow mode backward */
                        follow_mode--;
                        if (follow_mode < 0) follow_mode = 2 + num_bodies;
                        need_redraw = 1;
                        break;
                    case '1':  /* Inner planets view */
                        view_scale = 2.0;
                        follow_mode = 0;  /* Follow sun */
                        need_redraw = 1;
                        break;
                    case '2':  /* Outer planets view */
                        view_scale = 35.0;
                        follow_mode = 0;  /* Follow sun */
                        need_redraw = 1;
                        break;
                    case '3':  /* Moon view - follow earth+moon */
                        view_scale = 0.006;  /* Moon orbit is 0.00257 AU radius */
                        follow_mode = 2;  /* Earth+moon barycenter */
                        need_redraw = 1;
                        break;
                    case '4':  /* Follow earth */
                        follow_mode = 1;
                        need_redraw = 1;
                        break;
                    case 'd': case 'D':  /* Toggle 2D/3D view */
                        view_3d = !view_3d;
                        if (view_3d && tilt_angle < 5) tilt_angle = 30;  /* Default tilt */
                        need_redraw = 1;
                        break;
                    case 's':  /* Increase size scale (×3) */
                        size_scale *= 3.0;
                        if (size_scale > 1e12) size_scale = 1e12;
                        need_redraw = 1;
                        break;
                    case 'S':  /* Decrease size scale (÷3) */
                        size_scale /= 3.0;
                        if (size_scale < 1.0) size_scale = 1.0;
                        need_redraw = 1;
                        break;
                    case 'w':  /* Tilt up (more top-down) */
                        if (view_3d) {
                            tilt_angle -= 10;
                            if (tilt_angle < 0) tilt_angle = 0;
                        }
                        need_redraw = 1;
                        break;
                    case 'e':  /* Tilt down (more edge-on) */
                        if (view_3d) {
                            tilt_angle += 10;
                            if (tilt_angle > 90) tilt_angle = 90;
                        }
                        need_redraw = 1;
                        break;
                }
            }
        }
    }
    
    free(screen);
    free(colors);
    iplot_restore_term();
    return 1;
}

#else
/* Stubs for non-Unix */
int do_iscatter(double *xs, double *ys, int *iters, int n, int max_iter)
{
    (void)xs; (void)ys; (void)iters; (void)n; (void)max_iter;
    printf("iscatter requires Unix terminal\n");
    return 0;
}

int do_tscatter(const char **bodies, int num_bodies, const char *follow)
{
    (void)bodies; (void)num_bodies; (void)follow;
    printf("tscatter requires Unix terminal\n");
    return 0;
}
#endif

/* Planet data function for calculator - returns array of 10 values:
 * [0] x - heliocentric X position (AU)
 * [1] y - heliocentric Y position (AU, in ecliptic plane)
 * [2] z - height above ecliptic (AU)
 * [3] radius - physical radius (AU)
 * [4] mass - mass (Earth masses)
 * [5] period - orbital period (days)
 * [6] a - semi-major axis (AU)
 * [7] e - eccentricity
 * [8] i - inclination (degrees)
 * [9] rot_period - rotation period (hours, negative=retrograde)
 */
#define PLANET_DATA_SIZE 10

void fn_planet_data(const char *name, double t_days, double *data)
{
    const planet_data_t *p = find_planet_data(name);
    double x, y, z;
    int i;
    
    /* Initialize to zero */
    for (i = 0; i < PLANET_DATA_SIZE; i++) data[i] = 0;
    
    if (!p) return;
    
    /* Get 3D position */
    planet_pos_3d(name, t_days, &x, &y, &z);
    
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = p->radius;
    data[4] = p->mass;
    data[5] = p->period;
    data[6] = p->a;
    data[7] = p->e;
    data[8] = p->i;
    data[9] = p->rot_period;
}

/* Legacy 2D wrapper for backward compatibility */
void fn_planet(const char *name, double t_days, double *x, double *y)
{
    planet_pos(name, t_days, x, y);
}
