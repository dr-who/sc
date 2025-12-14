/* plot.c - Function plotting for DOS VGA graphics
 * C89 compliant for Watcom C / DOS
 * Uses bisection to progressively fill plot - press key to stop
 */
#include "sc.h"

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

int do_plot(const char *func_name, apf *xmin, apf *xmax)
{
    int func_idx, x;
    apf x_range, x_step, temp;
    apf ymin, ymax, y_range;
    apfc arg, result;
    int axis_x, axis_y;
    int prev_x, prev_y, first_point;
    
    /* Find function */
    func_idx = get_func_index(func_name);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        printf("Error: function '%s' not defined\n", func_name);
        return 0;
    }
    
    /* Compute x step */
    apf_sub(&x_range, xmax, xmin);
    apf_from_int(&temp, GFX_WIDTH - 1);
    apf_div(&x_step, &x_range, &temp);
    
    /* Quick y-range scan: just 8 sample points */
    {
        int first_valid = 0;
        for (x = 0; x < GFX_WIDTH; x += 40) {
            apf x_val, idx_apf;
            apf_from_int(&idx_apf, x);
            apf_mul(&temp, &idx_apf, &x_step);
            apf_add(&x_val, xmin, &temp);
            
            apf_copy(&arg.re, &x_val);
            apf_zero(&arg.im);
            
            if (eval_user_func(&result, func_idx, &arg) &&
                !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
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
    
    /* Add 10% margin to y range */
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
    
    /* Switch to graphics immediately */
    set_video_mode(0x13);
    
    /* Draw axes */
    if (axis_y >= 0 && axis_y < GFX_HEIGHT) {
        draw_hline(0, GFX_WIDTH - 1, axis_y, COLOR_CYAN);
    }
    if (axis_x >= 0 && axis_x < GFX_WIDTH) {
        draw_vline(axis_x, 0, GFX_HEIGHT - 1, COLOR_CYAN);
    }
    
    /* Plot points left to right */
    prev_x = -1;
    prev_y = -1;
    first_point = 1;
    
    for (x = 0; x < GFX_WIDTH; x++) {
        apf x_val, idx_apf, norm_y;
        int py;
        
        /* Check for cancel */
        if (kbhit()) {
            getch();
            break;
        }
        
        /* Compute x value */
        apf_from_int(&idx_apf, x);
        apf_mul(&temp, &idx_apf, &x_step);
        apf_add(&x_val, xmin, &temp);
        
        apf_copy(&arg.re, &x_val);
        apf_zero(&arg.im);
        
        /* Evaluate function */
        if (eval_user_func(&result, func_idx, &arg) &&
            !apf_isnan(&result.re) && !apf_isinf(&result.re) &&
            apf_iszero(&result.im)) {
            
            /* Map to screen Y */
            apf_sub(&temp, &result.re, &ymin);
            apf_div(&norm_y, &temp, &y_range);
            apf_from_int(&temp, GFX_HEIGHT - 1);
            apf_mul(&norm_y, &norm_y, &temp);
            py = (GFX_HEIGHT - 1) - (int)apf_to_long(&norm_y);
            
            /* Update y range if needed (auto-scale) */
            if (apf_cmp(&result.re, &ymin) < 0 || apf_cmp(&result.re, &ymax) > 0) {
                /* Point out of range - just skip drawing */
                if (py < 0) py = -1;
                if (py >= GFX_HEIGHT) py = -1;
            }
            
            if (py >= 0 && py < GFX_HEIGHT) {
                /* Draw line from previous point */
                if (!first_point && prev_y >= 0) {
                    draw_line(prev_x, prev_y, x, py, COLOR_YELLOW);
                } else {
                    bios_put_pixel(x, py, COLOR_YELLOW);
                }
                prev_x = x;
                prev_y = py;
                first_point = 0;
            } else {
                /* Gap in curve */
                first_point = 1;
            }
        } else {
            /* Invalid point - gap in curve */
            first_point = 1;
        }
    }
    
    /* Wait for key to exit */
    while (kbhit()) getch();
    getch();
    
    set_video_mode(0x03);
    dos_init_screen();
    
    return 1;
}

/* Text mode plot for DOS - shared code */
static int text_plot_impl(const char *func_name, apf *xmin, apf *xmax)
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
    
    func_idx = get_func_index(func_name);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        printf("Error: function '%s' not defined\n", func_name);
        return 0;
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

int do_textplot(const char *func_name, apf *xmin, apf *xmax)
{
    return text_plot_impl(func_name, xmin, xmax);
}

#else
/* Non-DOS: ASCII plot only */

int do_plot(const char *func_name, apf *xmin, apf *xmax)
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
    
    func_idx = get_func_index(func_name);
    if (func_idx < 0 || !user_funcs[func_idx].defined) {
        printf("Error: function '%s' not defined\n", func_name);
        return 0;
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

int do_textplot(const char *func_name, apf *xmin, apf *xmax)
{
    return do_plot(func_name, xmin, xmax);
}

/*
 * Lorenz Attractor - the classic chaotic system
 * dx/dt = sigma * (y - x)
 * dy/dt = x * (rho - z) - y
 * dz/dt = x * y - beta * z
 *
 * Classic parameters: sigma=10, rho=28, beta=8/3
 * Projects onto x-z plane for the butterfly shape
 */

/* Text-mode Lorenz attractor */
int do_lorenz_text(double sigma, double rho, double beta, int steps)
{
    double x = 1.0, y = 1.0, z = 1.0;
    double dt = 0.01;
    double xmin = -25, xmax = 25, zmin = 0, zmax = 55;
    int width = 70, height = 24;
    static char screen[30][80];
    int i, row, col;
    int px, pz;
    
    /* Clear screen buffer */
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
        screen[row][width] = '\0';
    }
    
    printf("\nLorenz Attractor (sigma=%.1f, rho=%.1f, beta=%.3f)\n", sigma, rho, beta);
    printf("Projecting x-z plane, %d iterations\n\n", steps);
    
    /* Integrate and plot */
    for (i = 0; i < steps; i++) {
        double dx, dy, dz;
        
        /* RK4 integration step */
        {
            double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
            double tx, ty, tz;
            
            /* k1 */
            k1x = sigma * (y - x);
            k1y = x * (rho - z) - y;
            k1z = x * y - beta * z;
            
            /* k2 */
            tx = x + 0.5 * dt * k1x;
            ty = y + 0.5 * dt * k1y;
            tz = z + 0.5 * dt * k1z;
            k2x = sigma * (ty - tx);
            k2y = tx * (rho - tz) - ty;
            k2z = tx * ty - beta * tz;
            
            /* k3 */
            tx = x + 0.5 * dt * k2x;
            ty = y + 0.5 * dt * k2y;
            tz = z + 0.5 * dt * k2z;
            k3x = sigma * (ty - tx);
            k3y = tx * (rho - tz) - ty;
            k3z = tx * ty - beta * tz;
            
            /* k4 */
            tx = x + dt * k3x;
            ty = y + dt * k3y;
            tz = z + dt * k3z;
            k4x = sigma * (ty - tx);
            k4y = tx * (rho - tz) - ty;
            k4z = tx * ty - beta * tz;
            
            dx = (k1x + 2*k2x + 2*k3x + k4x) / 6.0;
            dy = (k1y + 2*k2y + 2*k3y + k4y) / 6.0;
            dz = (k1z + 2*k2z + 2*k3z + k4z) / 6.0;
        }
        
        x += dt * dx;
        y += dt * dy;
        z += dt * dz;
        
        /* Map to screen (x-z projection) */
        px = (int)((x - xmin) / (xmax - xmin) * (width - 1));
        pz = (int)((z - zmin) / (zmax - zmin) * (height - 1));
        pz = (height - 1) - pz;  /* Flip Y */
        
        if (px >= 0 && px < width && pz >= 0 && pz < height) {
            /* Use different chars for density */
            if (screen[pz][px] == ' ') screen[pz][px] = '.';
            else if (screen[pz][px] == '.') screen[pz][px] = ':';
            else if (screen[pz][px] == ':') screen[pz][px] = '+';
            else if (screen[pz][px] == '+') screen[pz][px] = '*';
            else if (screen[pz][px] == '*') screen[pz][px] = '#';
        }
    }
    
    /* Print the result */
    printf("z=%.0f |", zmax);
    for (col = 0; col < width; col++) printf("%c", screen[0][col]);
    printf("\n");
    
    for (row = 1; row < height - 1; row++) {
        printf("      |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("z=%.0f  |", zmin);
    for (col = 0; col < width; col++) printf("%c", screen[height-1][col]);
    printf("\n");
    
    printf("      +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("      x=%.0f", xmin);
    for (col = 0; col < width - 16; col++) printf(" ");
    printf("x=%.0f\n\n", xmax);
    
    return 1;
}

/* Rossler Attractor - another classic strange attractor
 * dx/dt = -y - z
 * dy/dt = x + a*y
 * dz/dt = b + z*(x - c)
 * Classic parameters: a=0.2, b=0.2, c=5.7
 */
int do_rossler_text(double a, double b, double c, int steps)
{
    double x = 1.0, y = 1.0, z = 1.0;
    double dt = 0.02;
    double xmin = -15, xmax = 15, ymin = -15, ymax = 15;
    int width = 70, height = 24;
    static char screen[30][80];
    int i, row, col;
    int px, py;
    
    /* Clear screen buffer */
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
        screen[row][width] = '\0';
    }
    
    printf("\nRossler Attractor (a=%.1f, b=%.1f, c=%.1f)\n", a, b, c);
    printf("Projecting x-y plane, %d iterations\n\n", steps);
    
    /* Integrate and plot */
    for (i = 0; i < steps; i++) {
        double dx, dy, dz;
        
        /* Simple Euler - Rossler is less stiff */
        dx = -y - z;
        dy = x + a * y;
        dz = b + z * (x - c);
        
        x += dt * dx;
        y += dt * dy;
        z += dt * dz;
        
        /* Map to screen (x-y projection) */
        px = (int)((x - xmin) / (xmax - xmin) * (width - 1));
        py = (int)((y - ymin) / (ymax - ymin) * (height - 1));
        py = (height - 1) - py;
        
        if (px >= 0 && px < width && py >= 0 && py < height) {
            if (screen[py][px] == ' ') screen[py][px] = '.';
            else if (screen[py][px] == '.') screen[py][px] = ':';
            else if (screen[py][px] == ':') screen[py][px] = '+';
            else if (screen[py][px] == '+') screen[py][px] = '*';
            else if (screen[py][px] == '*') screen[py][px] = '#';
        }
    }
    
    /* Print */
    printf("y=%.0f |", ymax);
    for (col = 0; col < width; col++) printf("%c", screen[0][col]);
    printf("\n");
    
    for (row = 1; row < height - 1; row++) {
        printf("      |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("y=%.0f |", ymin);
    for (col = 0; col < width; col++) printf("%c", screen[height-1][col]);
    printf("\n");
    
    printf("      +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("      x=%.0f", xmin);
    for (col = 0; col < width - 16; col++) printf(" ");
    printf("x=%.0f\n\n", xmax);
    
    return 1;
}

/* Parametric curve plotting: x(t), y(t) 
 * For Lissajous figures, spirals, etc.
 */
int do_parametric_text(const char *xfunc, const char *yfunc, 
                       double tmin, double tmax, int steps)
{
    int xidx, yidx;
    double t, dt;
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    int width = 70, height = 24;
    static char screen[30][80];
    static double xs[2000], ys[2000];
    int i, row, col, valid_count = 0;
    
    /* Find function indices */
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
    
    if (steps > 2000) steps = 2000;
    dt = (tmax - tmin) / steps;
    
    /* First pass: compute values and find range */
    for (i = 0; i < steps; i++) {
        apfc tval, xres, yres;
        t = tmin + i * dt;
        
        apf_from_double(&tval.re, t);
        apf_zero(&tval.im);
        
        if (eval_user_func(&xres, xidx, &tval) && 
            eval_user_func(&yres, yidx, &tval) &&
            !apf_isnan(&xres.re) && !apf_isnan(&yres.re) &&
            !apf_isinf(&xres.re) && !apf_isinf(&yres.re)) {
            
            xs[i] = apf_to_double(&xres.re);
            ys[i] = apf_to_double(&yres.re);
            
            if (xs[i] < xmin) xmin = xs[i];
            if (xs[i] > xmax) xmax = xs[i];
            if (ys[i] < ymin) ymin = ys[i];
            if (ys[i] > ymax) ymax = ys[i];
            valid_count++;
        } else {
            xs[i] = 1e30;  /* Mark invalid */
            ys[i] = 1e30;
        }
    }
    
    if (valid_count == 0) {
        printf("No valid points computed\n");
        return 0;
    }
    
    /* Add margin */
    {
        double xmargin = (xmax - xmin) * 0.05;
        double ymargin = (ymax - ymin) * 0.05;
        if (xmargin < 0.1) xmargin = 0.1;
        if (ymargin < 0.1) ymargin = 0.1;
        xmin -= xmargin; xmax += xmargin;
        ymin -= ymargin; ymax += ymargin;
    }
    
    /* Clear screen */
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            screen[row][col] = ' ';
        }
    }
    
    printf("\nParametric: %s(t), %s(t)  t=[%.2f, %.2f]\n\n", xfunc, yfunc, tmin, tmax);
    
    /* Plot points */
    for (i = 0; i < steps; i++) {
        int px, py;
        if (xs[i] > 1e20) continue;  /* Skip invalid */
        
        px = (int)((xs[i] - xmin) / (xmax - xmin) * (width - 1));
        py = (int)((ys[i] - ymin) / (ymax - ymin) * (height - 1));
        py = (height - 1) - py;
        
        if (px >= 0 && px < width && py >= 0 && py < height) {
            if (screen[py][px] == ' ') screen[py][px] = '.';
            else if (screen[py][px] == '.') screen[py][px] = '*';
            else screen[py][px] = '#';
        }
    }
    
    /* Print */
    for (row = 0; row < height; row++) {
        if (row == 0) printf("%7.1f |", ymax);
        else if (row == height - 1) printf("%7.1f |", ymin);
        else printf("        |");
        for (col = 0; col < width; col++) printf("%c", screen[row][col]);
        printf("\n");
    }
    
    printf("        +");
    for (col = 0; col < width; col++) printf("-");
    printf("\n");
    printf("        %-7.1f", xmin);
    for (col = 0; col < width - 14; col++) printf(" ");
    printf("%7.1f\n\n", xmax);
    
    return 1;
}

#endif
