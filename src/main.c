/* main.c - sc - Scientific Calculator entry point
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"
#include "help.h"

/* Forward declaration for evaluating a single expression
 * Returns: 0 for success/true, 1 for boolean false */
int eval_expr_line(const char *line, int quiet);

/* Check if input is from a pipe/redirect (non-interactive) */
static int is_pipe_input(void)
{
#ifdef HAVE_TERMIOS
    extern int isatty(int);
    extern int fileno(FILE *);
    return !isatty(fileno(stdin));
#else
    return 0;  /* Assume interactive on DOS/Windows */
#endif
}

int main(int argc, char *argv[])
{
    static char input[MAX_INPUT];
    int len;
    int interactive;
    int last_bool_result = 1;  /* Track last boolean result for exit code (1=true, 0=false) */

    init_user_funcs();
    init_variables();
    
    /* Initialize last_ans */
    last_ans.type = VAL_SCALAR;
    apf_zero(&last_ans.v.scalar.re);
    apf_zero(&last_ans.v.scalar.im);
    last_ans_valid = 0;
    
    /* Command-line expression: sc "1+2" or sc 1+2 */
    if (argc > 1) {
        int i;
        int exit_code;
        static char expr[MAX_INPUT];
        expr[0] = '\0';
        
        /* Concatenate all arguments with spaces */
        for (i = 1; i < argc; i++) {
            if (i > 1) strcat(expr, " ");
            strncat(expr, argv[i], MAX_INPUT - strlen(expr) - 2);
        }
        
        exit_code = eval_expr_line(expr, 1);  /* quiet=1: just print result, no "= " prefix */
        return exit_code;  /* 0 for true/non-boolean, 1 for false */
    }
    
    /* Check for pipe input: echo "1+2" | sc */
    interactive = !is_pipe_input();
    
    /* Interactive mode (or pipe mode with same logic, just no banner/prompts) */
    if (interactive) {
        dos_init_screen();
        printf("sc - Scientific Calculator.  Type 'help' for commands, 'vars' to list variables.\n");
    }

    for (;;) {
        if (interactive) {
            print_prompt();
            fflush(stdout);
        }

        len = read_line(input, MAX_INPUT);
        if (len < 0) { 
            if (interactive) printf("\n"); 
            break; 
        }
        if (input[0] == '\0') continue;

        /* Commands */
        if (str_eq(input, "quit") || str_eq(input, "exit")) break;
        if (str_eq(input, "help")) { print_help(); continue; }
        if (str_starts(input, "help ")) {
            const char *fn = input + 5;
            while (*fn == ' ') fn++;
            if (*fn) show_function_help(fn);
            else print_help();
            continue;
        }
        if (str_starts(input, "demo ")) {
            const char *fn = input + 5;
            while (*fn == ' ') fn++;
            if (*fn) show_function_demo(fn);
            continue;
        }
        if (str_eq(input, "test")) { run_tests(); continue; }
        if (str_eq(input, "bench")) { run_bench(); continue; }
        if (str_starts(input, "mode")) { handle_mode(input); continue; }
        
        /* Angle mode commands */
        if (str_eq(input, "deg") || str_eq(input, "degrees")) {
            angle_mode = ANGLE_DEG;
            printf("Angle mode: degrees\n");
            continue;
        }
        if (str_eq(input, "rad") || str_eq(input, "radians")) {
            angle_mode = ANGLE_RAD;
            printf("Angle mode: radians\n");
            continue;
        }
        if (str_eq(input, "grad") || str_eq(input, "gradians") || str_eq(input, "gon")) {
            angle_mode = ANGLE_GRAD;
            printf("Angle mode: gradians\n");
            continue;
        }
        
        /* Digits command: set display precision */
        if (str_starts(input, "digits")) {
            const char *arg = input + 6;
            while (*arg == ' ') arg++;
            if (*arg == '\0') {
                /* Show current setting */
                printf("Display: %d digits (max %d, internal %d)\n", 
                       display_digits, MAX_DISPLAY_DIGITS, INTERNAL_DIGITS);
            } else {
                int n = 0;
                while (*arg >= '0' && *arg <= '9') {
                    n = n * 10 + (*arg - '0');
                    arg++;
                }
                if (n <= 0) {
                    /* digits 0 means show maximum safe digits */
                    display_digits = MAX_DISPLAY_DIGITS;
                    printf("Display: %d digits (maximum safe)\n", display_digits);
                } else {
                    if (n > MAX_DISPLAY_DIGITS) {
                        printf("Warning: max safe display is %d (using %d exposes guard digit)\n",
                               MAX_DISPLAY_DIGITS, INTERNAL_DIGITS);
                        if (n > INTERNAL_DIGITS) n = INTERNAL_DIGITS;
                    }
                    display_digits = n;
                    printf("Display: %d digits\n", n);
                }
            }
            continue;
        }
        
        /* Round command: set output rounding mode */
        /* Only match "round" alone or "round <mode>", not "round(x)" function */
        if (str_starts(input, "round") && (input[5] == '\0' || input[5] == ' ')) {
            const char *arg = input + 5;
            while (*arg == ' ') arg++;
            if (*arg == '\0') {
                /* Show current mode */
                int mode = apf_get_round_mode();
                printf("Rounding: %s\n", 
                       mode == APF_ROUND_NEAREST ? "nearest (ties to even)" :
                       mode == APF_ROUND_CEIL ? "ceil (toward +infinity)" :
                       mode == APF_ROUND_FLOOR ? "floor (toward -infinity)" :
                       mode == APF_ROUND_TRUNC ? "trunc (toward zero)" :
                       "away (ties away from zero)");
            } else if (str_eq(arg, "nearest") || str_eq(arg, "even")) {
                apf_set_round_mode(APF_ROUND_NEAREST);
                printf("Rounding: nearest (ties to even)\n");
            } else if (str_eq(arg, "ceil") || str_eq(arg, "up")) {
                apf_set_round_mode(APF_ROUND_CEIL);
                printf("Rounding: ceil (toward +infinity)\n");
            } else if (str_eq(arg, "floor") || str_eq(arg, "down")) {
                apf_set_round_mode(APF_ROUND_FLOOR);
                printf("Rounding: floor (toward -infinity)\n");
            } else if (str_eq(arg, "trunc") || str_eq(arg, "zero")) {
                apf_set_round_mode(APF_ROUND_TRUNC);
                printf("Rounding: trunc (toward zero)\n");
            } else if (str_eq(arg, "away")) {
                apf_set_round_mode(APF_ROUND_AWAY);
                printf("Rounding: away (ties away from zero)\n");
            } else {
                printf("Rounding modes: nearest, ceil, floor, trunc, away\n");
            }
            continue;
        }
        
        if (str_eq(input, "funcs") || str_eq(input, "functions")) {
            int i, count = 0;
            
            /* Show user-defined functions first */
            printf("User-defined functions:\n");
            for (i = 0; i < MAX_FUNCTIONS; i++) {
                if (user_funcs[i].defined) {
                    printf("  %c(%c) = %s\n", 'a'+i, user_funcs[i].param, user_funcs[i].body);
                    count++;
                }
            }
            if (count == 0) printf("  (none)\n");
            
            /* Show all built-in functions from help database */
            list_all_functions();
            continue;
        }
        
        if (str_eq(input, "vars") || str_eq(input, "variables")) {
            int i, count = 0;
            static char vbuf[256];
            /* Show scalar variables */
            for (i = 0; i < MAX_SCALAR_VARS; i++) {
                if (scalar_vars[i].defined) {
                    apfc_to_str(vbuf, sizeof(vbuf), &scalar_vars[i].val, display_digits);
                    printf("  %c = %s\n", 'a'+i, vbuf);
                    count++;
                }
            }
            /* Show matrix variables */
            for (i = 0; i < MAX_MATRIX_VARS; i++) {
                if (matrix_vars[i].defined) {
                    printf("  %c = [%dx%d matrix]\n", matrix_vars[i].name, 
                           matrix_vars[i].val.rows,
                           matrix_vars[i].val.cols);
                    count++;
                }
            }
            if (count == 0) printf("No variables defined.\n");
            continue;
        }
        
        if (str_eq(input, "date") || str_eq(input, "time")) {
            print_date_time();
            continue;
        }
        
#ifdef HAVE_RPN
        /* RPN mode commands */
        if (str_eq(input, "rpn")) {
            extern int rpn_mode;
            extern void rpn_init(void);
            rpn_mode = 1;
            rpn_init();
            printf("RPN mode enabled. Stack commands: drop swap dup clear .s\n");
            continue;
        }
        if (str_eq(input, "alg") || str_eq(input, "algebraic")) {
            extern int rpn_mode;
            rpn_mode = 0;
            printf("Algebraic (infix) mode enabled.\n");
            continue;
        }
#endif

#ifdef HAVE_SOLVER
        /* Quadratic solver */
        if (str_starts(input, "quad ")) {
            extern void cmd_quadratic(const char *args);
            cmd_quadratic(input + 5);
            continue;
        }
        if (str_starts(input, "solve ")) {
            extern void cmd_solve(const char *args);
            cmd_solve(input + 6);
            continue;
        }
#endif

#ifdef HAVE_STATS
        /* Statistics commands */
        if (str_starts(input, "stat")) {
            extern void cmd_stats(const char *args);
            const char *p = input + 4;
            while (*p == ' ') p++;
            cmd_stats(p);
            continue;
        }
        if (str_starts(input, "data+")) {
            extern void stats_add(const apf *x);
            extern void stats_add2(const apf *x, const apf *y);
            extern int stats_count(void);
            apf x, y;
            const char *p = input + 5;
            const char *sep;
            while (*p == ' ') p++;
            if (*p) {
                /* Parse first number */
                apf_from_str(&x, p);
                /* Skip first number */
                while (*p && *p != ' ' && *p != ',' && *p != '\t') p++;
                /* Check for separator (comma or space) */
                sep = p;
                while (*sep == ' ' || *sep == '\t') sep++;
                if (*sep == ',') {
                    sep++;
                    while (*sep == ' ' || *sep == '\t') sep++;
                }
                if (*sep && *sep != '\n' && *sep != '\r') {
                    apf_from_str(&y, sep);
                    stats_add2(&x, &y);
                } else {
                    stats_add(&x);
                }
                printf("n = %d\n", stats_count());
            }
            continue;
        }
        if (str_starts(input, "data-")) {
            extern void stats_sub(const apf *x);
            extern int stats_count(void);
            apf x;
            const char *p = input + 5;
            while (*p == ' ') p++;
            if (*p) {
                apf_from_str(&x, p);
                stats_sub(&x);
                printf("n = %d\n", stats_count());
            }
            continue;
        }
#endif

#ifdef HAVE_TVM
        /* Time Value of Money */
        if (str_starts(input, "tvm")) {
            extern void cmd_tvm(const char *args);
            const char *p = input + 3;
            while (*p == ' ') p++;
            cmd_tvm(p);
            continue;
        }
#endif

#ifdef HAVE_NEWTON
        /* Newton-Raphson solver */
        if (str_starts(input, "root ") || str_starts(input, "nsolve ")) {
            extern void cmd_newton(const char *args);
            const char *p = input;
            if (*p == 'r') p += 5; else p += 7;
            cmd_newton(p);
            continue;
        }
#endif

#ifdef HAVE_ORBITAL
        /* Orbital mechanics */
        if (str_starts(input, "orb")) {
            extern void cmd_orbital(const char *args);
            const char *p = input + 3;
            while (*p == ' ') p++;
            cmd_orbital(p);
            continue;
        }
#endif
        
        /* Lorenz attractor: lorenz [sigma rho beta] [steps] */
        if (str_starts(input, "lorenz")) {
            double sigma = 10.0, rho = 28.0, beta = 8.0/3.0;
            int steps = 5000;
            char *p = input + 6;
            
            /* Parse optional parameters */
            while (*p == ' ') p++;
            if (*p) {
                if (sscanf(p, "%lf %lf %lf %d", &sigma, &rho, &beta, &steps) < 3) {
                    if (sscanf(p, "%d", &steps) != 1) {
                        steps = 5000;
                    }
                }
            }
            if (steps < 100) steps = 100;
            if (steps > 50000) steps = 50000;
            
            do_lorenz_text((long)sigma, (long)rho, (long)(beta*3), 3, steps);
            continue;
        }
        
        /* Rossler attractor: rossler [a b c] [steps] */
        if (str_starts(input, "rossler")) {
            double a = 0.2, b = 0.2, c = 5.7;
            int steps = 5000;
            char *p = input + 7;
            
            while (*p == ' ') p++;
            if (*p) {
                if (sscanf(p, "%lf %lf %lf %d", &a, &b, &c, &steps) < 3) {
                    if (sscanf(p, "%d", &steps) != 1) {
                        steps = 5000;
                    }
                }
            }
            if (steps < 100) steps = 100;
            if (steps > 50000) steps = 50000;
            
            do_rossler_text((long)(a*10), 10, (long)(b*10), 10, (long)(c*10), 10, steps);
            continue;
        }
        
        /* Parametric plot: parametric xfunc yfunc [tmin:tmax] [steps] */
        if (str_starts(input, "parametric ")) {
            char xfunc[16], yfunc[16];
            double tmin = 0, tmax = 6.283185;  /* 0 to 2*pi default */
            int steps = 500;
            char *p = input + 11;
            int xi = 0, yi = 0;
            
            /* Skip whitespace and get first function name */
            while (*p == ' ') p++;
            while (*p && *p != ' ' && xi < 15) {
                xfunc[xi++] = *p++;
            }
            xfunc[xi] = '\0';
            
            /* Get second function name */
            while (*p == ' ') p++;
            while (*p && *p != ' ' && yi < 15) {
                yfunc[yi++] = *p++;
            }
            yfunc[yi] = '\0';
            
            /* Parse optional range and steps */
            while (*p == ' ') p++;
            if (*p) {
                double v1, v2;
                int n;
                if (sscanf(p, "%lf:%lf %d", &v1, &v2, &n) >= 2) {
                    tmin = v1;
                    tmax = v2;
                    if (n > 0) steps = n;
                } else if (sscanf(p, "%d", &n) == 1) {
                    steps = n;
                }
            }
            if (steps < 50) steps = 50;
            if (steps > 2000) steps = 2000;
            
            if (xi == 0 || yi == 0) {
                printf("Usage: parametric xfunc yfunc [tmin:tmax] [steps]\n");
                printf("Example: parametric f g 0:6.28 500\n");
            } else {
                apf tmin_apf, tmax_apf; { apf t; apf_from_int(&tmin_apf, (long)(tmin * 1000)); apf_from_int(&t, 1000); apf_div(&tmin_apf, &tmin_apf, &t); apf_from_int(&tmax_apf, (long)(tmax * 1000)); apf_div(&tmax_apf, &tmax_apf, &t); }
                do_parametric_text(xfunc, yfunc, &tmin_apf, &tmax_apf, steps);
            }
            continue;
        }
        
        /* Plot command: plot EXPR [range] - can plot any expression with x
         * Examples: plot sin(x)
         *           plot x^2+x -1:1
         *           plot f(x) -5:5
         *           plot isprime(x) 0:100
         */
        if (str_starts(input, "plot ") || str_starts(input, "textplot ")) {
            int text_mode = str_starts(input, "textplot ");
            static char expr_buf[256];
            char *p = input + (text_mode ? 9 : 5);
            char *range_start = NULL;
            int ei = 0;
            apf xmin, xmax;
            
            /* Skip leading whitespace */
            while (*p == ' ') p++;
            
            /* Find where range starts: look for " -N:" or " N:" pattern
             * Scan through to find space followed by digit/minus then eventually colon */
            {
                char *scan = p;
                while (*scan) {
                    if (*scan == ' ') {
                        char *after_space = scan + 1;
                        while (*after_space == ' ') after_space++;
                        /* Check if this looks like a range: starts with -, digit, or xrange */
                        if (str_starts(after_space, "xrange")) {
                            range_start = after_space;
                            break;
                        }
                        if (*after_space == '-' || (*after_space >= '0' && *after_space <= '9')) {
                            /* Scan forward to see if there's a colon */
                            char *colon_check = after_space;
                            while (*colon_check && *colon_check != ' ' && *colon_check != ':') {
                                colon_check++;
                            }
                            if (*colon_check == ':') {
                                range_start = after_space;
                                break;
                            }
                        }
                    }
                    scan++;
                }
            }
            
            /* Copy expression (everything before range) */
            if (range_start) {
                /* Copy up to the space before range */
                while (p < range_start - 1 && ei < 255) {
                    expr_buf[ei++] = *p++;
                }
                /* Trim trailing spaces */
                while (ei > 0 && expr_buf[ei-1] == ' ') ei--;
            } else {
                /* No range - copy whole expression */
                while (*p && ei < 255) {
                    expr_buf[ei++] = *p++;
                }
                /* Trim trailing spaces */
                while (ei > 0 && expr_buf[ei-1] == ' ') ei--;
            }
            expr_buf[ei] = '\0';
            
            /* If expression is just a single-letter function name followed by (x), strip the (x)
             * This is for user-defined functions like f(x), g(x), etc.
             * Don't strip from built-in functions like sin(x), cos(x), etc.
             */
            if (ei == 4 && expr_buf[0] >= 'a' && expr_buf[0] <= 'z' &&
                expr_buf[1] == '(' && expr_buf[2] == 'x' && expr_buf[3] == ')') {
                /* Single letter function like f(x) - strip to just "f" */
                expr_buf[1] = '\0';
            }
            
            /* Default range */
            apf_from_int(&xmin, -10);
            apf_from_int(&xmax, 10);
            
            /* Parse range if present */
            if (range_start && *range_start) {
                const char *saved = input_ptr;
                p = range_start;
                
                if (str_starts(p, "xrange(") || str_starts(p, "xrange (")) {
                    /* Old format: xrange(-5,5) */
                    while (*p && *p != '(') p++;
                    if (*p == '(') {
                        p++;
                        input_ptr = p;
                        next_token();
                        {
                            apfc xmin_c;
                            if (parse_expr(&xmin_c)) {
                                apf_copy(&xmin, &xmin_c.re);
                            }
                        }
                        if (current_token.type == TOK_COMMA) {
                            next_token();
                            {
                                apfc xmax_c;
                                if (parse_expr(&xmax_c)) {
                                    apf_copy(&xmax, &xmax_c.re);
                                }
                            }
                        }
                    }
                } else {
                    /* New format: -5:0.1:5 or -5:5 (xmin:step:xmax or xmin:xmax) */
                    input_ptr = p;
                    next_token();
                    {
                        apfc val1, val2, val3;
                        int have_val1 = 0, have_val2 = 0, have_val3 = 0;
                        
                        /* Parse first value (xmin) */
                        if (parse_expr(&val1)) {
                            have_val1 = 1;
                        }
                        
                        /* Check for colon */
                        if (current_token.type == TOK_COLON) {
                            next_token();
                            /* Parse second value */
                            if (parse_expr(&val2)) {
                                have_val2 = 1;
                            }
                            
                            /* Check for another colon */
                            if (current_token.type == TOK_COLON) {
                                next_token();
                                /* Parse third value (xmax) */
                                if (parse_expr(&val3)) {
                                    have_val3 = 1;
                                }
                            }
                        }
                        
                        /* Interpret values */
                        if (have_val1 && have_val2 && have_val3) {
                            /* -5:0.1:5 format: val1=xmin, val2=step (ignored), val3=xmax */
                            apf_copy(&xmin, &val1.re);
                            apf_copy(&xmax, &val3.re);
                        } else if (have_val1 && have_val2) {
                            /* -5:5 format: val1=xmin, val2=xmax */
                            apf_copy(&xmin, &val1.re);
                            apf_copy(&xmax, &val2.re);
                        } else if (have_val1) {
                            /* Single value: treat as xmax, xmin=0 */
                            apf_zero(&xmin);
                            apf_copy(&xmax, &val1.re);
                        }
                    }
                }
                input_ptr = saved;
            }
            
#ifdef HAVE_CONIO
            if (text_mode) {
                do_textplot(expr_buf, &xmin, &xmax);
            } else {
                do_plot(expr_buf, 'x', &xmin, &xmax);
            }
#else
            /* On Linux, both plot and textplot use text mode */
            do_plot(expr_buf, 'x', &xmin, &xmax);
#endif
            continue;
        }

        /* Function definition: f(x) = expr */
        if (len >= 5 && input[0] >= 'a' && input[0] <= 'z' &&
            input[1] == '(' && input[2] >= 'a' && input[2] <= 'z' && input[3] == ')') {
            char *eq = input + 4;
            while (*eq == ' ' || *eq == '\t') eq++;
            if (*eq == '=') {
                int fi = input[0] - 'a';
                char param = input[2];
                const char *body = eq + 1;
                while (*body == ' ' || *body == '\t') body++;
                if (*body == '\0') {
                    user_funcs[fi].defined = 0;
                    printf("Function %c deleted.\n", input[0]);
                } else {
                    user_funcs[fi].defined = 1;
                    user_funcs[fi].param = param;
                    strncpy(user_funcs[fi].body, body, MAX_FUNC_BODY - 1);
                    user_funcs[fi].body[MAX_FUNC_BODY - 1] = '\0';
                    printf("%c(%c) = %s\n", input[0], param, body);
                }
                continue;
            }
        }

#ifdef HAVE_RPN
        /* RPN mode - process as stack operations */
        {
            extern int rpn_mode;
            if (rpn_mode) {
                extern void rpn_eval_line(const char *line);
                rpn_eval_line(input);
                continue;
            }
        }
#endif

        /* Parse and evaluate expression(s) */
        input_ptr = input;
        next_token();

        /* Check for for loop: for v in start:end; body; end */
        if (current_token.type == TOK_FOR) {
            int var_idx;
            apfc start_val, end_val, step_val;
            long start_i, end_i, step_i, loop_i;
            const char *body_start;
            const char *body_end;
            
            next_token();
            /* Accept single-letter variable (TOK_FUNC) or 'i' (TOK_IMAG) */
            if (current_token.type == TOK_IMAG) {
                var_idx = 'i' - 'a';
            } else if (current_token.type == TOK_FUNC && strlen(current_token.func_name) == 1) {
                var_idx = get_var_index(current_token.func_name);
            } else {
                printf("Error: expected single-letter variable after 'for'\n");
                continue;
            }
            if (var_idx < 0) {
                printf("Error: invalid loop variable\n");
                continue;
            }
            
            next_token();
            if (current_token.type != TOK_IN) {
                printf("Error: expected 'in' after variable\n");
                continue;
            }
            
            next_token();
            if (!parse_expr(&start_val)) continue;
            
            if (current_token.type != TOK_COLON) {
                printf("Error: expected ':' in range\n");
                continue;
            }
            next_token();
            
            if (!parse_expr(&end_val)) continue;
            
            /* Check for step (start:step:end) */
            apf_from_int(&step_val.re, 1);
            apf_zero(&step_val.im);
            if (current_token.type == TOK_COLON) {
                /* end_val is actually step, need to read actual end */
                step_val = end_val;
                next_token();
                if (!parse_expr(&end_val)) continue;
            }
            
            if (current_token.type != TOK_SEMI) {
                printf("Error: expected ';' after range\n");
                continue;
            }
            
            /* Save body start position - it's right after the semicolon */
            /* input_ptr is currently pointing past the ';', at the body start */
            body_start = input_ptr;
            next_token();
            
            /* Find 'end' keyword */
            body_end = NULL;
            {
                int depth = 1;
                const char *p = input_ptr;
                while (*p && depth > 0) {
                    /* Skip to next keyword */
                    while (*p && *p != 'e' && *p != 'f') p++;
                    if (strncmp(p, "end", 3) == 0 && !isalpha(p[3])) {
                        depth--;
                        if (depth == 0) {
                            body_end = p;
                        }
                        p += 3;
                    } else if (strncmp(p, "for", 3) == 0 && !isalpha(p[3])) {
                        depth++;
                        p += 3;
                    } else if (*p) {
                        p++;
                    }
                }
            }
            
            if (!body_end) {
                printf("Error: missing 'end'\n");
                continue;
            }
            
            start_i = apf_to_long(&start_val.re);
            end_i = apf_to_long(&end_val.re);
            step_i = apf_to_long(&step_val.re);
            
            if (step_i == 0) {
                printf("Error: step cannot be zero\n");
                continue;
            }
            
            /* Execute loop */
            for (loop_i = start_i; (step_i > 0) ? (loop_i <= end_i) : (loop_i >= end_i); loop_i += step_i) {
                value_t loop_result;
                value_t loop_var_val;
                const char *body_ptr;
                
                /* Set loop variable */
                loop_var_val.type = VAL_SCALAR;
                apf_from_int(&loop_var_val.v.scalar.re, loop_i);
                apf_zero(&loop_var_val.v.scalar.im);
                set_var_value(var_idx, &loop_var_val);
                
                /* Execute body */
                body_ptr = body_start;
                while (body_ptr < body_end) {
                    input_ptr = body_ptr;
                    next_token();
                    
                    if (current_token.type == TOK_END || 
                        current_token.type == TOK_ENDFOR ||
                        input_ptr >= body_end) {
                        break;
                    }
                    
                    /* Handle assignment or expression (including 'i' which is TOK_IMAG) */
                    if (current_token.type == TOK_FUNC || current_token.type == TOK_IMAG) {
                        char var_name[16];
                        int vidx;
                        const char *saved = input_ptr;
                        token_t saved_tok = current_token;
                        
                        if (current_token.type == TOK_IMAG) {
                            strcpy(var_name, "i");
                        } else {
                            strcpy(var_name, current_token.func_name);
                        }
                        vidx = get_var_index(var_name);
                        next_token();
                        
                        if (vidx >= 0 && current_token.type == TOK_ASSIGN) {
                            next_token();
                            if (parse_value(&loop_result)) {
                                set_var_value(vidx, &loop_result);
                            }
                        } else {
                            input_ptr = saved;
                            current_token = saved_tok;
                            parse_value(&loop_result);
                        }
                    } else {
                        parse_value(&loop_result);
                    }
                    
                    /* Skip to next semicolon */
                    if (current_token.type == TOK_SEMI) {
                        body_ptr = input_ptr;
                    } else {
                        break;
                    }
                }
            }
            continue;
        }

        {
            int done = 0;
            value_t last_result;
            int has_result = 0;
            int is_equality_check = 0;
            
            result_is_boolean = 0;  /* Reset for this expression */
            last_result.type = VAL_SCALAR;
            apf_zero(&last_result.v.scalar.re);
            apf_zero(&last_result.v.scalar.im);
            
            while (!done) {
                /* Handle multi-value assignment: [A,B,C] = func(x) */
                if (current_token.type == TOK_LBRACKET) {
                    char mvnames[8][32]; int mvidx[8]; int mvnum = 0;
                    const char *mvsaved = input_ptr; token_t mvtok = current_token;
                    int ismv = 0;
                    next_token();
                    while (mvnum < 8 && current_token.type != TOK_RBRACKET && current_token.type != TOK_END) {
                        if (current_token.type == TOK_FUNC) { strcpy(mvnames[mvnum], current_token.func_name); mvidx[mvnum] = get_var_index(mvnames[mvnum]); mvnum++; }
                        else if (current_token.type == TOK_IMAG) { strcpy(mvnames[mvnum], "i"); mvidx[mvnum] = get_var_index("i"); mvnum++; }
                        else break;
                        next_token(); if (current_token.type == TOK_COMMA) next_token();
                    }
                    if (current_token.type == TOK_RBRACKET && mvnum > 0) { next_token(); if (current_token.type == TOK_ASSIGN) ismv = 1; }
                    if (ismv) {
                        next_token();
                        if (current_token.type == TOK_FUNC) {
                            char fn[32]; strcpy(fn, current_token.func_name); next_token();
                            if (current_token.type == TOK_LPAREN) {
                                next_token();
                                if (strcmp(fn, "svd") == 0 && mvnum == 3) {
                                    value_t arg; matrix_t U, S, V;
                                    if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                        if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                        if (mat_svd(&U, &S, &V, &arg.v.matrix)) {
                                            value_t v; v.type = VAL_MATRIX;
                                            mat_copy(&v.v.matrix, &U); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &S); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &V); if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix);
                                            has_result = 0;
                                        } else { printf("Error: SVD failed\n"); done = 1; }
                                    } else { printf("Error: svd syntax\n"); done = 1; }
                                } else if (strcmp(fn, "qr") == 0 && mvnum == 2) {
                                    value_t arg; matrix_t Q, R;
                                    if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                        if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                        if (mat_qr(&Q, &R, &arg.v.matrix)) {
                                            value_t v; v.type = VAL_MATRIX;
                                            mat_copy(&v.v.matrix, &Q); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &R); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                            has_result = 0;
                                        } else { printf("Error: QR failed\n"); done = 1; }
                                    } else { printf("Error: qr syntax\n"); done = 1; }
                                } else if (strcmp(fn, "lu") == 0 && mvnum >= 2) {
                                    value_t arg; matrix_t L, Um; int perm[MAT_MAX_ROWS];
                                    if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                        if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                        if (mat_lu(&L, &Um, perm, &arg.v.matrix)) {
                                            value_t v; int ii; v.type = VAL_MATRIX;
                                            mat_copy(&v.v.matrix, &L); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &Um); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                            if (mvnum == 3) { int nn = arg.v.matrix.rows; mat_zero(&v.v.matrix, nn, nn);
                                                for (ii = 0; ii < nn; ii++) apf_from_int(&MAT_AT(&v.v.matrix, ii, perm[ii]).re, 1);
                                                if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix); }
                                            has_result = 0;
                                        } else { printf("Error: LU failed\n"); done = 1; }
                                    } else { printf("Error: lu syntax\n"); done = 1; }
                                } else if (strcmp(fn, "pca") == 0 && mvnum == 3) {
                                    value_t arg; matrix_t coeff, score, latent;
                                    if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                        if (arg.type != VAL_MATRIX) { printf("Error: pca requires matrix\n"); done = 1; }
                                        else if (mat_pca(&coeff, &score, &latent, &arg.v.matrix)) {
                                            value_t v; v.type = VAL_MATRIX;
                                            mat_copy(&v.v.matrix, &coeff); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &score); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &latent); if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix);
                                            has_result = 0;
                                        } else { printf("Error: PCA failed\n"); done = 1; }
                                    } else { printf("Error: pca syntax\n"); done = 1; }
                                } else if (strcmp(fn, "kmeans") == 0 && mvnum == 2) {
                                    value_t arg, karg; matrix_t idx_m, cent; int k;
                                    if (parse_value(&arg) && current_token.type == TOK_COMMA) { next_token();
                                        if (parse_value(&karg) && current_token.type == TOK_RPAREN) { next_token();
                                            if (arg.type != VAL_MATRIX) { printf("Error: kmeans requires matrix\n"); done = 1; }
                                            else { k = (karg.type == VAL_SCALAR) ? apf_to_long(&karg.v.scalar.re) : 2;
                                                if (mat_kmeans(&idx_m, &cent, &arg.v.matrix, k)) {
                                                    value_t v; v.type = VAL_MATRIX;
                                                    mat_copy(&v.v.matrix, &idx_m); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                                    mat_copy(&v.v.matrix, &cent); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                                    has_result = 0;
                                                } else { printf("Error: kmeans failed\n"); done = 1; } }
                                        } else { printf("Error: kmeans(X,k) syntax\n"); done = 1; }
                                    } else { printf("Error: kmeans(X,k) syntax\n"); done = 1; }
                                } else { printf("Error: unknown [a,b]=func '%s' with %d outputs\n", fn, mvnum); done = 1; }
                            } else { printf("Error: expected '('\n"); done = 1; }
                        } else { printf("Error: expected function\n"); done = 1; }
                        if (!done) { if (current_token.type == TOK_SEMI) { next_token(); if (current_token.type == TOK_END) done = 1; } else if (current_token.type == TOK_END) done = 1; }
                        continue;
                    } else { input_ptr = mvsaved; current_token = mvtok; }
                }
                /* Check for variable assignment: x = expr (including 'i' which is TOK_IMAG) */
                if (current_token.type == TOK_FUNC || current_token.type == TOK_IMAG) {
                    char var_name[16];
                    int var_idx;
                    const char *saved_pos = input_ptr;
                    token_t saved_tok = current_token;
                    
                    if (current_token.type == TOK_IMAG) {
                        strcpy(var_name, "i");
                    } else {
                        strcpy(var_name, current_token.func_name);
                    }
                    var_idx = get_var_index(var_name);
                    
                    next_token();
                    if (var_idx >= 0 && current_token.type == TOK_ASSIGN) {
                        /* Check if this is assignment or equality check */
                        next_token();
                        if (current_token.type == TOK_END || current_token.type == TOK_SEMI) {
                            /* Just "x =" with nothing after - error */
                            printf("Error: expected expression after =\n");
                            done = 1;
                            has_result = 0;
                            continue;
                        }
                        /* It's an assignment */
                        if (!parse_value(&last_result)) {
                            done = 1;
                            continue;
                        }
                        set_var_value(var_idx, &last_result);
                        has_result = 1;
                    } else {
                        /* Not an assignment, restore and parse as expression */
                        int leading_not = 0;
                        input_ptr = saved_pos;
                        current_token = saved_tok;
                        /* Track leading 'not' - will apply after comparison */
                        while (current_token.type == TOK_NOT) {
                            leading_not = !leading_not;
                            next_token();
                        }
                        if (!parse_value(&last_result)) {
                            done = 1;
                            continue;
                        }
                        has_result = 1;
                        
                        /* Handle comparison: after value, check for comparison operator */
                        if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                            current_token.type == TOK_NE || current_token.type == TOK_LT ||
                            current_token.type == TOK_LE || current_token.type == TOK_GT ||
                            current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                            value_t rhs;
                            int cmp_result = 0;
                            token_type_t op = current_token.type;
                            next_token();
                            if (!parse_value(&rhs)) {
                                printf("Error: expected expression after comparison operator\n");
                                done = 1;
                                has_result = 0;
                                continue;
                            }
                            if (last_result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                                int cmp = apf_cmp(&last_result.v.scalar.re, &rhs.v.scalar.re);
                                int eq = apf_eq(&last_result.v.scalar.re, &rhs.v.scalar.re) &&
                                         apf_eq(&last_result.v.scalar.im, &rhs.v.scalar.im);
                                switch (op) {
                                    case TOK_ASSIGN:
                                    case TOK_EQUAL: cmp_result = eq; break;
                                    case TOK_NE:    cmp_result = !eq; break;
                                    case TOK_LT:    cmp_result = (cmp < 0); break;
                                    case TOK_LE:    cmp_result = (cmp <= 0); break;
                                    case TOK_GT:    cmp_result = (cmp > 0); break;
                                    case TOK_GE:    cmp_result = (cmp >= 0); break;
                                    case TOK_APPROX: {
                                        apf abs_lhs, abs_rhs, upper, lower, factor_hi, factor_lo;
                                        apf_abs(&abs_lhs, &last_result.v.scalar.re);
                                        apf_abs(&abs_rhs, &rhs.v.scalar.re);
                                        apf_from_str(&factor_hi, "1.0500001");
                                        apf_from_str(&factor_lo, "0.9499999");
                                        apf_mul(&upper, &abs_lhs, &factor_hi);
                                        apf_mul(&lower, &abs_lhs, &factor_lo);
                                        if (apf_is_zero(&last_result.v.scalar.re)) {
                                            cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                                        } else if (last_result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                                            cmp_result = 0;
                                        } else {
                                            cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                                          apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                                        }
                                        break;
                                    }
                                    default: cmp_result = 0; break;
                                }
                            }
                            last_result.type = VAL_SCALAR;
                            apf_from_int(&last_result.v.scalar.re, cmp_result);
                            apf_zero(&last_result.v.scalar.im);
                            is_equality_check = 1;
                        }
                        
                        /* Apply leading 'not' after comparison */
                        if (leading_not && last_result.type == VAL_SCALAR) {
                            int val = apf_is_zero(&last_result.v.scalar.re) ? 1 : 0;
                            apf_from_int(&last_result.v.scalar.re, val);
                            apf_zero(&last_result.v.scalar.im);
                            is_equality_check = 1;
                        }
                    }
                } else {
                    /* Handle leading 'not' - track, will apply after comparison */
                    int leading_not = 0;
                    while (current_token.type == TOK_NOT) {
                        leading_not = !leading_not;
                        next_token();
                    }

                    if (parse_value(&last_result)) {
                        has_result = 1;

                        
                        /* Handle comparison */
                        if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                            current_token.type == TOK_NE || current_token.type == TOK_LT ||
                            current_token.type == TOK_LE || current_token.type == TOK_GT ||
                            current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                            value_t rhs;
                            int cmp_result = 0;
                            token_type_t op = current_token.type;
                            next_token();
                            if (!parse_value(&rhs)) {
                                printf("Error: expected expression after comparison operator\n");
                                done = 1;
                                has_result = 0;
                                continue;
                            }
                            if (last_result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                                int cmp = apf_cmp(&last_result.v.scalar.re, &rhs.v.scalar.re);
                                int eq = apf_eq(&last_result.v.scalar.re, &rhs.v.scalar.re) &&
                                         apf_eq(&last_result.v.scalar.im, &rhs.v.scalar.im);
                                switch (op) {
                                    case TOK_ASSIGN:
                                    case TOK_EQUAL: cmp_result = eq; break;
                                    case TOK_NE:    cmp_result = !eq; break;
                                    case TOK_LT:    cmp_result = (cmp < 0); break;
                                    case TOK_LE:    cmp_result = (cmp <= 0); break;
                                    case TOK_GT:    cmp_result = (cmp > 0); break;
                                    case TOK_GE:    cmp_result = (cmp >= 0); break;
                                    case TOK_APPROX: {
                                        apf abs_lhs, abs_rhs, upper, lower, factor_hi, factor_lo;
                                        apf_abs(&abs_lhs, &last_result.v.scalar.re);
                                        apf_abs(&abs_rhs, &rhs.v.scalar.re);
                                        apf_from_str(&factor_hi, "1.0500001");
                                        apf_from_str(&factor_lo, "0.9499999");
                                        apf_mul(&upper, &abs_lhs, &factor_hi);
                                        apf_mul(&lower, &abs_lhs, &factor_lo);
                                        if (apf_is_zero(&last_result.v.scalar.re)) {
                                            cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                                        } else if (last_result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                                            cmp_result = 0;
                                        } else {
                                            cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                                          apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                                        }
                                        break;
                                    }
                                    default: cmp_result = 0; break;
                                }
                            }
                            last_result.type = VAL_SCALAR;
                            apf_from_int(&last_result.v.scalar.re, cmp_result);
                            apf_zero(&last_result.v.scalar.im);
                            is_equality_check = 1;
                        }
                        
                        /* Apply leading 'not' after comparison */
                        if (leading_not && last_result.type == VAL_SCALAR) {
                            int val = apf_is_zero(&last_result.v.scalar.re) ? 1 : 0;
                            apf_from_int(&last_result.v.scalar.re, val);
                            apf_zero(&last_result.v.scalar.im);
                            is_equality_check = 1;
                        }
                    } else {
                        done = 1;
                        continue;
                    }
                }
                
                /* Check for comparison: expr op expr where op is =, ==, <>, <, <=, >, >=, ~= */
                if (has_result && (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                                   current_token.type == TOK_NE || current_token.type == TOK_LT ||
                                   current_token.type == TOK_LE || current_token.type == TOK_GT ||
                                   current_token.type == TOK_GE || current_token.type == TOK_APPROX)) {
                    value_t rhs;
                    int cmp_result = 0;
                    token_type_t op = current_token.type;
                    next_token();
                    if (!parse_value(&rhs)) {
                        printf("Error: expected expression after comparison operator\n");
                        done = 1;
                        has_result = 0;
                        continue;
                    }
                    /* Compare the two values (scalars only for now) */
                    if (last_result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                        /* Only compare real parts for ordering */
                        int cmp = apf_cmp(&last_result.v.scalar.re, &rhs.v.scalar.re);
                        int eq = apf_eq(&last_result.v.scalar.re, &rhs.v.scalar.re) &&
                                 apf_eq(&last_result.v.scalar.im, &rhs.v.scalar.im);
                        switch (op) {
                            case TOK_ASSIGN:
                            case TOK_EQUAL: cmp_result = eq; break;
                            case TOK_NE:    cmp_result = !eq; break;
                            case TOK_LT:    cmp_result = (cmp < 0); break;
                            case TOK_LE:    cmp_result = (cmp <= 0); break;
                            case TOK_GT:    cmp_result = (cmp > 0); break;
                            case TOK_GE:    cmp_result = (cmp >= 0); break;
                            case TOK_APPROX: {
                                /* ~= : approximately equal within 5% of |lhs| */
                                /* Use slightly wider bounds to handle floating point rounding */
                                apf abs_lhs, abs_rhs, upper, lower, factor_hi, factor_lo;
                                apf_abs(&abs_lhs, &last_result.v.scalar.re);
                                apf_abs(&abs_rhs, &rhs.v.scalar.re);
                                apf_from_str(&factor_hi, "1.0500001");
                                apf_from_str(&factor_lo, "0.9499999");
                                apf_mul(&upper, &abs_lhs, &factor_hi);
                                apf_mul(&lower, &abs_lhs, &factor_lo);
                                if (apf_is_zero(&last_result.v.scalar.re)) {
                                    cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                                } else if (last_result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                                    cmp_result = 0;
                                } else {
                                    cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                                  apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                                }
                                break;
                            }
                            default:        cmp_result = 0; break;
                        }
                    } else {
                        /* Matrices: only equality is supported */
                        cmp_result = 0;
                    }
                    /* Store result as 1 or 0 */
                    last_result.type = VAL_SCALAR;
                    apf_from_int(&last_result.v.scalar.re, cmp_result);
                    apf_zero(&last_result.v.scalar.im);
                    is_equality_check = 1;
                    has_result = 1;
                }
                
                /* Handle boolean operators: and, or, xor */
                while (has_result && (current_token.type == TOK_AND || current_token.type == TOK_OR || 
                                      current_token.type == TOK_XOR)) {
                    token_type_t bool_op = current_token.type;
                    int lhs_bool, rhs_bool, bool_result;
                    value_t rhs;
                    
                    /* Get left side as boolean (non-zero = true) */
                    if (last_result.type != VAL_SCALAR) {
                        printf("Error: boolean operators require scalar values\n");
                        done = 1;
                        has_result = 0;
                        break;
                    }
                    lhs_bool = !apf_is_zero(&last_result.v.scalar.re);
                    
                    next_token();
                    
                    /* Handle 'not' prefix on right side */
                    {
                        int negate = 0;
                        while (current_token.type == TOK_NOT) {
                            negate = !negate;
                            next_token();
                        }
                        
                        if (!parse_value(&rhs)) {
                            printf("Error: expected expression after boolean operator\n");
                            done = 1;
                            has_result = 0;
                            break;
                        }
                        
                        /* Check for comparison on right side */
                        if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                            current_token.type == TOK_NE || current_token.type == TOK_LT ||
                            current_token.type == TOK_LE || current_token.type == TOK_GT ||
                            current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                            value_t rhs2;
                            int cmp_result = 0;
                            token_type_t op = current_token.type;
                            next_token();
                            if (!parse_value(&rhs2)) {
                                printf("Error: expected expression after comparison\n");
                                done = 1;
                                has_result = 0;
                                break;
                            }
                            if (rhs.type == VAL_SCALAR && rhs2.type == VAL_SCALAR) {
                                int cmp = apf_cmp(&rhs.v.scalar.re, &rhs2.v.scalar.re);
                                int eq = apf_eq(&rhs.v.scalar.re, &rhs2.v.scalar.re) &&
                                         apf_eq(&rhs.v.scalar.im, &rhs2.v.scalar.im);
                                switch (op) {
                                    case TOK_ASSIGN:
                                    case TOK_EQUAL: cmp_result = eq; break;
                                    case TOK_NE:    cmp_result = !eq; break;
                                    case TOK_LT:    cmp_result = (cmp < 0); break;
                                    case TOK_LE:    cmp_result = (cmp <= 0); break;
                                    case TOK_GT:    cmp_result = (cmp > 0); break;
                                    case TOK_GE:    cmp_result = (cmp >= 0); break;
                                    case TOK_APPROX: {
                                        apf abs_lhs2, abs_rhs2, upper, lower, fhi, flo;
                                        apf_abs(&abs_lhs2, &rhs.v.scalar.re);
                                        apf_abs(&abs_rhs2, &rhs2.v.scalar.re);
                                        apf_from_str(&fhi, "1.0500001");
                                        apf_from_str(&flo, "0.9499999");
                                        apf_mul(&upper, &abs_lhs2, &fhi);
                                        apf_mul(&lower, &abs_lhs2, &flo);
                                        if (apf_is_zero(&rhs.v.scalar.re)) {
                                            cmp_result = apf_is_zero(&rhs2.v.scalar.re) ? 1 : 0;
                                        } else if (rhs.v.scalar.re.sign != rhs2.v.scalar.re.sign) {
                                            cmp_result = 0;
                                        } else {
                                            cmp_result = (apf_cmp(&abs_rhs2, &lower) >= 0 && 
                                                          apf_cmp(&abs_rhs2, &upper) <= 0) ? 1 : 0;
                                        }
                                        break;
                                    }
                                    default: cmp_result = 0; break;
                                }
                            }
                            rhs.type = VAL_SCALAR;
                            apf_from_int(&rhs.v.scalar.re, cmp_result);
                            apf_zero(&rhs.v.scalar.im);
                        }
                        
                        if (rhs.type != VAL_SCALAR) {
                            printf("Error: boolean operators require scalar values\n");
                            done = 1;
                            has_result = 0;
                            break;
                        }
                        
                        rhs_bool = !apf_is_zero(&rhs.v.scalar.re);
                        if (negate) rhs_bool = !rhs_bool;
                    }
                    
                    /* Apply boolean operator */
                    switch (bool_op) {
                        case TOK_AND: bool_result = lhs_bool && rhs_bool; break;
                        case TOK_OR:  bool_result = lhs_bool || rhs_bool; break;
                        case TOK_XOR: bool_result = lhs_bool != rhs_bool; break;
                        default:      bool_result = 0; break;
                    }
                    
                    last_result.type = VAL_SCALAR;
                    apf_from_int(&last_result.v.scalar.re, bool_result);
                    apf_zero(&last_result.v.scalar.im);
                    is_equality_check = 1;
                }
                
                /* Update ans immediately after each expression (before semicolon check) */
                if (has_result) {
                    last_ans = last_result;
                    last_ans_valid = 1;
                }
                
                if (current_token.type == TOK_SEMI) {
                    next_token();
                    if (current_token.type == TOK_END) {
                        done = 1;
                    }
                    is_equality_check = 0;  /* Reset for next expression */
                } else if (current_token.type == TOK_END) {
                    done = 1;
                } else {
                    printf("Error: unexpected characters\n");
                    done = 1;
                    has_result = 0;
                }
            }
            
            if (has_result) {
                if ((is_equality_check || result_is_boolean) && last_result.type == VAL_SCALAR) {
                    /* Print true/false for equality checks and boolean functions */
                    int val = apf_is_zero(&last_result.v.scalar.re) ? 0 : 1;
                    last_bool_result = val;  /* Track for exit code */
                    if (interactive) {
                        printf("= %s\n", val ? "true" : "false");
                    } else {
                        printf("%s\n", val ? "true" : "false");
                    }
                    result_is_boolean = 0;  /* Reset for next expression */
                } else {
                    print_value_ex(&last_result, !interactive);
                    last_bool_result = 1;  /* Non-boolean: reset to true */
                }
            }
        }
    }

    if (interactive) printf("Goodbye!\n");
    return last_bool_result ? 0 : 1;  /* Exit 0 for true/non-boolean, 1 for false */
}

/* Evaluate a single expression line (for command-line/pipe mode) */
int eval_expr_line(const char *line, int quiet)
{
    value_t result;
    static char buf[512];
    int done = 0;
    int has_result = 0;
    int is_equality_check = 0;
    
    result_is_boolean = 0;  /* Reset for this expression */
    result.type = VAL_SCALAR;
    apf_zero(&result.v.scalar.re);
    apf_zero(&result.v.scalar.im);
    
    input_ptr = line;
    next_token();
    
    while (!done) {
        /* Handle multi-value assignment: [A,B,C] = func(x) */
        if (current_token.type == TOK_LBRACKET) {
            char mvnames[8][32]; int mvidx[8]; int mvnum = 0;
            const char *mvsaved = input_ptr; token_t mvtok = current_token;
            int ismv = 0;
            next_token();
            while (mvnum < 8 && current_token.type != TOK_RBRACKET && current_token.type != TOK_END) {
                if (current_token.type == TOK_FUNC) { strcpy(mvnames[mvnum], current_token.func_name); mvidx[mvnum] = get_var_index(mvnames[mvnum]); mvnum++; }
                else if (current_token.type == TOK_IMAG) { strcpy(mvnames[mvnum], "i"); mvidx[mvnum] = get_var_index("i"); mvnum++; }
                else break;
                next_token(); if (current_token.type == TOK_COMMA) next_token();
            }
            if (current_token.type == TOK_RBRACKET && mvnum > 0) { next_token(); if (current_token.type == TOK_ASSIGN) ismv = 1; }
            if (ismv) {
                next_token();
                if (current_token.type == TOK_FUNC) {
                    char fn[32]; strcpy(fn, current_token.func_name); next_token();
                    if (current_token.type == TOK_LPAREN) {
                        next_token();
                        if (strcmp(fn, "svd") == 0 && mvnum == 3) {
                            value_t arg; matrix_t U, S, V;
                            if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                if (mat_svd(&U, &S, &V, &arg.v.matrix)) {
                                    value_t v; v.type = VAL_MATRIX;
                                    mat_copy(&v.v.matrix, &U); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &S); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &V); if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix);
                                    has_result = 0;
                                } else { printf("Error: SVD failed\n"); done = 1; }
                            } else { printf("Error: svd syntax\n"); done = 1; }
                        } else if (strcmp(fn, "qr") == 0 && mvnum == 2) {
                            value_t arg; matrix_t Q, R;
                            if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                if (mat_qr(&Q, &R, &arg.v.matrix)) {
                                    value_t v; v.type = VAL_MATRIX;
                                    mat_copy(&v.v.matrix, &Q); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &R); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                    has_result = 0;
                                } else { printf("Error: QR failed\n"); done = 1; }
                            } else { printf("Error: qr syntax\n"); done = 1; }
                        } else if (strcmp(fn, "lu") == 0 && mvnum >= 2) {
                            value_t arg; matrix_t L, Um; int perm[MAT_MAX_ROWS];
                            if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                if (arg.type != VAL_MATRIX) { mat_zero(&arg.v.matrix, 1, 1); MAT_AT(&arg.v.matrix, 0, 0) = arg.v.scalar; }
                                if (mat_lu(&L, &Um, perm, &arg.v.matrix)) {
                                    value_t v; int ii; v.type = VAL_MATRIX;
                                    mat_copy(&v.v.matrix, &L); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &Um); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                    if (mvnum == 3) { int nn = arg.v.matrix.rows; mat_zero(&v.v.matrix, nn, nn);
                                        for (ii = 0; ii < nn; ii++) apf_from_int(&MAT_AT(&v.v.matrix, ii, perm[ii]).re, 1);
                                        if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix); }
                                    has_result = 0;
                                } else { printf("Error: LU failed\n"); done = 1; }
                            } else { printf("Error: lu syntax\n"); done = 1; }
                        } else if (strcmp(fn, "pca") == 0 && mvnum == 3) {
                            value_t arg; matrix_t coeff, score, latent;
                            if (parse_value(&arg) && current_token.type == TOK_RPAREN) { next_token();
                                if (arg.type != VAL_MATRIX) { printf("Error: pca requires matrix\n"); done = 1; }
                                else if (mat_pca(&coeff, &score, &latent, &arg.v.matrix)) {
                                    value_t v; v.type = VAL_MATRIX;
                                    mat_copy(&v.v.matrix, &coeff); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &score); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                    mat_copy(&v.v.matrix, &latent); if (mvidx[2] >= 0) set_var_value(mvidx[2], &v); else set_named_matrix(mvnames[2], &v.v.matrix);
                                    has_result = 0;
                                } else { printf("Error: PCA failed\n"); done = 1; }
                            } else { printf("Error: pca syntax\n"); done = 1; }
                        } else if (strcmp(fn, "kmeans") == 0 && mvnum == 2) {
                            value_t arg, karg; matrix_t idx_m, cent; int k;
                            if (parse_value(&arg) && current_token.type == TOK_COMMA) { next_token();
                                if (parse_value(&karg) && current_token.type == TOK_RPAREN) { next_token();
                                    if (arg.type != VAL_MATRIX) { printf("Error: kmeans requires matrix\n"); done = 1; }
                                    else { k = (karg.type == VAL_SCALAR) ? apf_to_long(&karg.v.scalar.re) : 2;
                                        if (mat_kmeans(&idx_m, &cent, &arg.v.matrix, k)) {
                                            value_t v; v.type = VAL_MATRIX;
                                            mat_copy(&v.v.matrix, &idx_m); if (mvidx[0] >= 0) set_var_value(mvidx[0], &v); else set_named_matrix(mvnames[0], &v.v.matrix);
                                            mat_copy(&v.v.matrix, &cent); if (mvidx[1] >= 0) set_var_value(mvidx[1], &v); else set_named_matrix(mvnames[1], &v.v.matrix);
                                            has_result = 0;
                                        } else { printf("Error: kmeans failed\n"); done = 1; } }
                                } else { printf("Error: kmeans(X,k) syntax\n"); done = 1; }
                            } else { printf("Error: kmeans(X,k) syntax\n"); done = 1; }
                        } else { printf("Error: unknown [a,b]=func '%s' with %d outputs\n", fn, mvnum); done = 1; }
                    } else { printf("Error: expected '('\n"); done = 1; }
                } else { printf("Error: expected function\n"); done = 1; }
                if (!done) { if (current_token.type == TOK_SEMI) { next_token(); if (current_token.type == TOK_END) done = 1; } else if (current_token.type == TOK_END) done = 1; }
                continue;
            } else { input_ptr = mvsaved; current_token = mvtok; }
        }
        /* Handle variable assignment */
        if (current_token.type == TOK_FUNC || current_token.type == TOK_IMAG) {
            char var_name[16];
            int var_idx;
            const char *saved_pos = input_ptr;
            token_t saved_tok = current_token;
            
            if (current_token.type == TOK_IMAG) {
                strcpy(var_name, "i");
            } else {
                strcpy(var_name, current_token.func_name);
            }
            var_idx = get_var_index(var_name);
            
            next_token();
            if (var_idx >= 0 && current_token.type == TOK_ASSIGN) {
                next_token();
                if (current_token.type == TOK_END || current_token.type == TOK_SEMI) {
                    printf("Error: expected expression after =\n");
                    done = 1;
                    has_result = 0;
                    continue;
                }
                if (!parse_value(&result)) {
                    done = 1;
                    continue;
                }
                set_var_value(var_idx, &result);
                has_result = 1;
            } else {
                int leading_not = 0;
                input_ptr = saved_pos;
                current_token = saved_tok;
                while (current_token.type == TOK_NOT) {
                    leading_not = !leading_not;
                    next_token();
                }
                if (!parse_value(&result)) {
                    done = 1;
                    continue;
                }
                has_result = 1;
                
                /* Handle comparison BEFORE applying not */
                if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                    current_token.type == TOK_NE || current_token.type == TOK_LT ||
                    current_token.type == TOK_LE || current_token.type == TOK_GT ||
                    current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                    value_t rhs;
                    int cmp_result = 0;
                    token_type_t op = current_token.type;
                    next_token();
                    if (!parse_value(&rhs)) {
                        printf("Error: expected expression after comparison\n");
                        done = 1;
                        has_result = 0;
                        continue;
                    }
                    if (result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                        int cmp = apf_cmp(&result.v.scalar.re, &rhs.v.scalar.re);
                        int eq = apf_eq(&result.v.scalar.re, &rhs.v.scalar.re) &&
                                 apf_eq(&result.v.scalar.im, &rhs.v.scalar.im);
                        switch (op) {
                            case TOK_ASSIGN:
                            case TOK_EQUAL: cmp_result = eq; break;
                            case TOK_NE:    cmp_result = !eq; break;
                            case TOK_LT:    cmp_result = (cmp < 0); break;
                            case TOK_LE:    cmp_result = (cmp <= 0); break;
                            case TOK_GT:    cmp_result = (cmp > 0); break;
                            case TOK_GE:    cmp_result = (cmp >= 0); break;
                            case TOK_APPROX: {
                                apf abs_lhs, abs_rhs, upper, lower, fhi, flo;
                                apf_abs(&abs_lhs, &result.v.scalar.re);
                                apf_abs(&abs_rhs, &rhs.v.scalar.re);
                                apf_from_str(&fhi, "1.0500001");
                                apf_from_str(&flo, "0.9499999");
                                apf_mul(&upper, &abs_lhs, &fhi);
                                apf_mul(&lower, &abs_lhs, &flo);
                                if (apf_is_zero(&result.v.scalar.re)) {
                                    cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                                } else if (result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                                    cmp_result = 0;
                                } else {
                                    cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                                  apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                                }
                                break;
                            }
                            default: cmp_result = 0; break;
                        }
                    }
                    result.type = VAL_SCALAR;
                    apf_from_int(&result.v.scalar.re, cmp_result);
                    apf_zero(&result.v.scalar.im);
                    is_equality_check = 1;
                }
                
                /* Apply leading not AFTER comparison */
                if (leading_not && result.type == VAL_SCALAR) {
                    int val = apf_is_zero(&result.v.scalar.re) ? 1 : 0;
                    apf_from_int(&result.v.scalar.re, val);
                    apf_zero(&result.v.scalar.im);
                    is_equality_check = 1;
                }
            }
        } else {
            int leading_not = 0;
            while (current_token.type == TOK_NOT) {
                leading_not = !leading_not;
                next_token();
            }
            if (parse_value(&result)) {
                has_result = 1;
                
                /* Handle comparison BEFORE applying not */
                if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                    current_token.type == TOK_NE || current_token.type == TOK_LT ||
                    current_token.type == TOK_LE || current_token.type == TOK_GT ||
                    current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                    value_t rhs;
                    int cmp_result = 0;
                    token_type_t op = current_token.type;
                    next_token();
                    if (!parse_value(&rhs)) {
                        printf("Error: expected expression after comparison\n");
                        done = 1;
                        has_result = 0;
                        continue;
                    }
                    if (result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                        int cmp = apf_cmp(&result.v.scalar.re, &rhs.v.scalar.re);
                        int eq = apf_eq(&result.v.scalar.re, &rhs.v.scalar.re) &&
                                 apf_eq(&result.v.scalar.im, &rhs.v.scalar.im);
                        switch (op) {
                            case TOK_ASSIGN:
                            case TOK_EQUAL: cmp_result = eq; break;
                            case TOK_NE:    cmp_result = !eq; break;
                            case TOK_LT:    cmp_result = (cmp < 0); break;
                            case TOK_LE:    cmp_result = (cmp <= 0); break;
                            case TOK_GT:    cmp_result = (cmp > 0); break;
                            case TOK_GE:    cmp_result = (cmp >= 0); break;
                            case TOK_APPROX: {
                                apf abs_lhs, abs_rhs, upper, lower, fhi, flo;
                                apf_abs(&abs_lhs, &result.v.scalar.re);
                                apf_abs(&abs_rhs, &rhs.v.scalar.re);
                                apf_from_str(&fhi, "1.0500001");
                                apf_from_str(&flo, "0.9499999");
                                apf_mul(&upper, &abs_lhs, &fhi);
                                apf_mul(&lower, &abs_lhs, &flo);
                                if (apf_is_zero(&result.v.scalar.re)) {
                                    cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                                } else if (result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                                    cmp_result = 0;
                                } else {
                                    cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                                  apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                                }
                                break;
                            }
                            default: cmp_result = 0; break;
                        }
                    }
                    result.type = VAL_SCALAR;
                    apf_from_int(&result.v.scalar.re, cmp_result);
                    apf_zero(&result.v.scalar.im);
                    is_equality_check = 1;
                }
                
                /* Apply leading not AFTER comparison */
                if (leading_not && result.type == VAL_SCALAR) {
                    int val = apf_is_zero(&result.v.scalar.re) ? 1 : 0;
                    apf_from_int(&result.v.scalar.re, val);
                    apf_zero(&result.v.scalar.im);
                    is_equality_check = 1;
                }
            } else {
                done = 1;
                continue;
            }
        }
        
        /* Check for comparison: expr op expr where op is =, ==, <>, <, <=, >, >=, ~= */
        if (has_result && (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                           current_token.type == TOK_NE || current_token.type == TOK_LT ||
                           current_token.type == TOK_LE || current_token.type == TOK_GT ||
                           current_token.type == TOK_GE || current_token.type == TOK_APPROX)) {
            value_t rhs;
            int cmp_result = 0;
            token_type_t op = current_token.type;
            next_token();
            if (!parse_value(&rhs)) {
                printf("Error: expected expression after comparison operator\n");
                done = 1;
                has_result = 0;
                continue;
            }
            /* Compare the two values (scalars only) */
            if (result.type == VAL_SCALAR && rhs.type == VAL_SCALAR) {
                int cmp = apf_cmp(&result.v.scalar.re, &rhs.v.scalar.re);
                int eq = apf_eq(&result.v.scalar.re, &rhs.v.scalar.re) &&
                         apf_eq(&result.v.scalar.im, &rhs.v.scalar.im);
                switch (op) {
                    case TOK_ASSIGN:
                    case TOK_EQUAL: cmp_result = eq; break;
                    case TOK_NE:    cmp_result = !eq; break;
                    case TOK_LT:    cmp_result = (cmp < 0); break;
                    case TOK_LE:    cmp_result = (cmp <= 0); break;
                    case TOK_GT:    cmp_result = (cmp > 0); break;
                    case TOK_GE:    cmp_result = (cmp >= 0); break;
                    case TOK_APPROX: {
                        /* ~= : approximately equal within 5% of |lhs| */
                        apf abs_lhs, abs_rhs, upper, lower, factor_hi, factor_lo;
                        apf_abs(&abs_lhs, &result.v.scalar.re);
                        apf_abs(&abs_rhs, &rhs.v.scalar.re);
                        apf_from_str(&factor_hi, "1.0500001");
                        apf_from_str(&factor_lo, "0.9499999");
                        apf_mul(&upper, &abs_lhs, &factor_hi);
                        apf_mul(&lower, &abs_lhs, &factor_lo);
                        if (apf_is_zero(&result.v.scalar.re)) {
                            cmp_result = apf_is_zero(&rhs.v.scalar.re) ? 1 : 0;
                        } else if (result.v.scalar.re.sign != rhs.v.scalar.re.sign) {
                            cmp_result = 0;
                        } else {
                            cmp_result = (apf_cmp(&abs_rhs, &lower) >= 0 && 
                                          apf_cmp(&abs_rhs, &upper) <= 0) ? 1 : 0;
                        }
                        break;
                    }
                    default:        cmp_result = 0; break;
                }
            } else {
                cmp_result = 0;
            }
            /* Store result as 1 or 0 */
            result.type = VAL_SCALAR;
            apf_from_int(&result.v.scalar.re, cmp_result);
            apf_zero(&result.v.scalar.im);
            is_equality_check = 1;
            has_result = 1;
        }
        
        /* Handle boolean operators: and, or, xor */
        while (has_result && (current_token.type == TOK_AND || current_token.type == TOK_OR || 
                              current_token.type == TOK_XOR)) {
            token_type_t bool_op = current_token.type;
            int lhs_bool, rhs_bool, bool_result;
            value_t rhs;
            
            if (result.type != VAL_SCALAR) {
                printf("Error: boolean operators require scalar values\n");
                done = 1;
                has_result = 0;
                break;
            }
            lhs_bool = !apf_is_zero(&result.v.scalar.re);
            
            next_token();
            
            {
                int negate = 0;
                while (current_token.type == TOK_NOT) {
                    negate = !negate;
                    next_token();
                }
                
                if (!parse_value(&rhs)) {
                    printf("Error: expected expression after boolean operator\n");
                    done = 1;
                    has_result = 0;
                    break;
                }
                
                /* Check for comparison on right side */
                if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                    current_token.type == TOK_NE || current_token.type == TOK_LT ||
                    current_token.type == TOK_LE || current_token.type == TOK_GT ||
                    current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                    value_t rhs2;
                    int cmp_result2 = 0;
                    token_type_t op2 = current_token.type;
                    next_token();
                    if (!parse_value(&rhs2)) {
                        printf("Error: expected expression after comparison\n");
                        done = 1;
                        has_result = 0;
                        break;
                    }
                    if (rhs.type == VAL_SCALAR && rhs2.type == VAL_SCALAR) {
                        int cmp2 = apf_cmp(&rhs.v.scalar.re, &rhs2.v.scalar.re);
                        int eq2 = apf_eq(&rhs.v.scalar.re, &rhs2.v.scalar.re) &&
                                 apf_eq(&rhs.v.scalar.im, &rhs2.v.scalar.im);
                        switch (op2) {
                            case TOK_ASSIGN:
                            case TOK_EQUAL: cmp_result2 = eq2; break;
                            case TOK_NE:    cmp_result2 = !eq2; break;
                            case TOK_LT:    cmp_result2 = (cmp2 < 0); break;
                            case TOK_LE:    cmp_result2 = (cmp2 <= 0); break;
                            case TOK_GT:    cmp_result2 = (cmp2 > 0); break;
                            case TOK_GE:    cmp_result2 = (cmp2 >= 0); break;
                            case TOK_APPROX: {
                                apf abs_l, abs_r, up, lo, fh, fl;
                                apf_abs(&abs_l, &rhs.v.scalar.re);
                                apf_abs(&abs_r, &rhs2.v.scalar.re);
                                apf_from_str(&fh, "1.0500001");
                                apf_from_str(&fl, "0.9499999");
                                apf_mul(&up, &abs_l, &fh);
                                apf_mul(&lo, &abs_l, &fl);
                                if (apf_is_zero(&rhs.v.scalar.re)) {
                                    cmp_result2 = apf_is_zero(&rhs2.v.scalar.re) ? 1 : 0;
                                } else if (rhs.v.scalar.re.sign != rhs2.v.scalar.re.sign) {
                                    cmp_result2 = 0;
                                } else {
                                    cmp_result2 = (apf_cmp(&abs_r, &lo) >= 0 && 
                                                  apf_cmp(&abs_r, &up) <= 0) ? 1 : 0;
                                }
                                break;
                            }
                            default: cmp_result2 = 0; break;
                        }
                    }
                    rhs.type = VAL_SCALAR;
                    apf_from_int(&rhs.v.scalar.re, cmp_result2);
                    apf_zero(&rhs.v.scalar.im);
                }
                
                if (rhs.type != VAL_SCALAR) {
                    printf("Error: boolean operators require scalar values\n");
                    done = 1;
                    has_result = 0;
                    break;
                }
                
                rhs_bool = !apf_is_zero(&rhs.v.scalar.re);
                if (negate) rhs_bool = !rhs_bool;
            }
            
            switch (bool_op) {
                case TOK_AND: bool_result = lhs_bool && rhs_bool; break;
                case TOK_OR:  bool_result = lhs_bool || rhs_bool; break;
                case TOK_XOR: bool_result = lhs_bool != rhs_bool; break;
                default:      bool_result = 0; break;
            }
            
            result.type = VAL_SCALAR;
            apf_from_int(&result.v.scalar.re, bool_result);
            apf_zero(&result.v.scalar.im);
            is_equality_check = 1;
        }
        
        /* Update ans immediately after each expression */
        if (has_result) {
            last_ans = result;
            last_ans_valid = 1;
        }
        
        if (current_token.type == TOK_SEMI) {
            next_token();
            if (current_token.type == TOK_END) {
                done = 1;
            }
            is_equality_check = 0;
        } else if (current_token.type == TOK_END) {
            done = 1;
        } else {
            printf("Error: unexpected characters\n");
            done = 1;
            has_result = 0;
        }
    }
    
    /* Print final result */
    if (has_result) {
        if ((is_equality_check || result_is_boolean) && result.type == VAL_SCALAR) {
            int val = apf_is_zero(&result.v.scalar.re) ? 0 : 1;
            printf("%s\n", val ? "true" : "false");
            result_is_boolean = 0;  /* Reset for next call */
            return val ? 0 : 1;  /* Exit 0 for true, 1 for false */
        } else if (result.type == VAL_SCALAR) {
            apfc_to_str(buf, sizeof(buf), &result.v.scalar, display_digits);
            if (quiet) {
                printf("%s\n", buf);
            } else {
                printf("= %s\n", buf);
            }
            return 0;  /* Non-boolean numeric: always exit 0 */
        } else {
            /* Matrix */
            if (!quiet) printf("= ");
            mat_print(&result.v.matrix);
        }
    }
    
    return 0;  /* Matrix or no result: exit 0 */
}
