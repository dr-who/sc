/* main.c - sc - Scientific Calculator entry point
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"

/* Forward declaration for evaluating a single expression */
static int eval_expr_line(const char *line, int quiet);

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
    int len, mant_digits;
    long max_bin_exp, min_bin_exp;
    int interactive;

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
        static char expr[MAX_INPUT];
        expr[0] = '\0';
        
        /* Concatenate all arguments with spaces */
        for (i = 1; i < argc; i++) {
            if (i > 1) strcat(expr, " ");
            strncat(expr, argv[i], MAX_INPUT - strlen(expr) - 2);
        }
        
        eval_expr_line(expr, 1);  /* quiet=1: just print result, no "= " prefix */
        return 0;
    }
    
    /* Check for pipe input: echo "1+2" | sc */
    interactive = !is_pipe_input();
    
    /* Interactive mode (or pipe mode with same logic, just no banner/prompts) */
    if (interactive) {
        dos_init_screen();

        /* Calculate precision info */
        mant_digits = (int)((long)AP_BITS * 301L / 1000L);
        max_bin_exp = APF_EXP_MAX - (AP_BITS - 1);
        min_bin_exp = APF_EXP_MIN;

        printf("\nsc - Scientific Calculator\n");
        printf("  Mantissa: %d bits (~%d decimal digits)\n", AP_BITS, mant_digits);
        printf("  Max value: ~2^%ld  Min positive: ~2^%ld\n", max_bin_exp, min_bin_exp);
        printf("  Matrices: up to %dx%d\n", MAT_MAX_ROWS, MAT_MAX_COLS);
        printf("Type 'help' for commands, 'vars' to list variables.\n\n");
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
        if (str_eq(input, "test")) { run_tests(); continue; }
        if (str_eq(input, "bench")) { run_bench(); continue; }
        if (str_starts(input, "mode")) { handle_mode(input); continue; }
        
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
        if (str_starts(input, "round")) {
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
            for (i = 0; i < MAX_FUNCTIONS; i++) {
                if (user_funcs[i].defined) {
                    printf("  %c(%c) = %s\n", 'a'+i, user_funcs[i].param, user_funcs[i].body);
                    count++;
                }
            }
            if (count == 0) printf("No functions defined.\n");
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
            while (*p == ' ') p++;
            if (*p) {
                const char *comma;
                apf_from_str(&x, p);
                comma = strchr(p, ',');
                if (comma) {
                    apf_from_str(&y, comma + 1);
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
            
            do_lorenz_text(sigma, rho, beta, steps);
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
            
            do_rossler_text(a, b, c, steps);
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
                do_parametric_text(xfunc, yfunc, tmin, tmax, steps);
            }
            continue;
        }
        
        /* Plot command: plot f(x) -5:0.1:5 or plot f(x) xrange(-5,5) */
        if (str_starts(input, "plot ") || str_starts(input, "textplot ")) {
            int text_mode = str_starts(input, "textplot ");
            char func_name[16];
            char *p = input + (text_mode ? 9 : 5);
            int fi = 0;
            apf xmin, xmax;
            
            /* Skip whitespace */
            while (*p == ' ') p++;
            
            /* Get function name */
            while (*p && *p != '(' && *p != ' ' && fi < 15) {
                func_name[fi++] = *p++;
            }
            func_name[fi] = '\0';
            
            /* Skip optional (x) */
            if (*p == '(') {
                while (*p && *p != ')') p++;
                if (*p == ')') p++;
            }
            
            /* Skip whitespace */
            while (*p == ' ') p++;
            
            /* Default range */
            apf_from_int(&xmin, -10);
            apf_from_int(&xmax, 10);
            
            /* Parse range: support -5:0.1:5 or xrange(-5,5) */
            if (*p && *p != '\0') {
                const char *saved = input_ptr;
                
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
                do_textplot(func_name, &xmin, &xmax);
            } else {
                do_plot(func_name, &xmin, &xmax);
            }
#else
            /* On Linux, both plot and textplot use text mode */
            do_plot(func_name, &xmin, &xmax);
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
            
            last_result.type = VAL_SCALAR;
            apf_zero(&last_result.v.scalar.re);
            apf_zero(&last_result.v.scalar.im);
            
            while (!done) {
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
                        next_token();
                        if (!parse_value(&last_result)) {
                            done = 1;
                            continue;
                        }
                        set_var_value(var_idx, &last_result);
                        has_result = 1;
                    } else {
                        /* Not an assignment, restore and parse as expression */
                        input_ptr = saved_pos;
                        current_token = saved_tok;
                        if (!parse_value(&last_result)) {
                            done = 1;
                            continue;
                        }
                        has_result = 1;
                    }
                } else if (parse_value(&last_result)) {
                    has_result = 1;
                } else {
                    done = 1;
                    continue;
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
                } else if (current_token.type == TOK_END) {
                    done = 1;
                } else {
                    printf("Error: unexpected characters\n");
                    done = 1;
                    has_result = 0;
                }
            }
            
            if (has_result) {
                print_value(&last_result);
            }
        }
    }

    if (interactive) printf("Goodbye!\n");
    return 0;
}

/* Evaluate a single expression line (for command-line/pipe mode) */
static int eval_expr_line(const char *line, int quiet)
{
    value_t result;
    static char buf[512];
    int done = 0;
    int has_result = 0;
    
    result.type = VAL_SCALAR;
    apf_zero(&result.v.scalar.re);
    apf_zero(&result.v.scalar.im);
    
    input_ptr = line;
    next_token();
    
    while (!done) {
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
                if (!parse_value(&result)) {
                    done = 1;
                    continue;
                }
                set_var_value(var_idx, &result);
                has_result = 1;
            } else {
                input_ptr = saved_pos;
                current_token = saved_tok;
                if (!parse_value(&result)) {
                    done = 1;
                    continue;
                }
                has_result = 1;
            }
        } else if (parse_value(&result)) {
            has_result = 1;
        } else {
            done = 1;
            continue;
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
        if (result.type == VAL_SCALAR) {
            apfc_to_str(buf, sizeof(buf), &result.v.scalar, display_digits);
            if (quiet) {
                printf("%s\n", buf);
            } else {
                printf("= %s\n", buf);
            }
        } else {
            /* Matrix */
            if (!quiet) printf("= ");
            mat_print(&result.v.matrix);
        }
    }
    
    return has_result;
}
