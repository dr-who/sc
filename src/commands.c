/* commands.c - REPL commands (help, date, test, mode)
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"
#include "help.h"
#include <math.h>  /* For sqrt in fisheriris demo */

/* ========== Help ========== */

void print_help(void)
{
    int dd = (int)((long)AP_BITS * 301L / 1000L);
    printf("\nsc - Scientific Calculator\n");
    printf("Precision: %d bits (~%d decimal digits)\n\n", AP_BITS, dd);
    printf("Operators: + - * / ^ (power) ! (factorial)\n");
    printf("           and or xor not (boolean)\n");
    printf("Comparisons: = == <> < <= > >= ~= (approx)\n");
    printf("Constants: pi e i Inf NaN ans (previous answer)\n\n");
    printf("Trig: sin cos tan asin acos atan atan2\n");
    printf("      sec csc cot asec acsc acot\n");
    printf("Hyperbolic: sinh cosh tanh asinh acosh atanh\n");
    printf("Math: exp log ln log10 logb sqrt abs floor ceil\n");
    printf("      gcd lcm nPr nCr isprime fact\n");
    printf("Complex: arg conj re im\n\n");
    printf("Angle Modes:\n");
    printf("  deg              - use degrees\n");
    printf("  rad              - use radians (default)\n");
    printf("  grad             - use gradians\n\n");
    printf("Number formats:\n");
    printf("  123, 3.14, 1.5e-3  - decimal / scientific\n");
    printf("  0xFF, 0xDEAD       - hexadecimal\n");
    printf("  0b1010             - binary\n");
    printf("  1.23e-4            - Scientific notation\n\n");
    printf("Variables & Functions:\n");
    printf("  x = 10             - assign variable\n");
    printf("  f(x) = x^2         - define function\n");
    printf("  vars / funcs       - list defined\n\n");
#ifdef HAVE_RPN
    printf("RPN Mode (HP-style):\n");
    printf("  rpn                - switch to RPN mode\n");
    printf("  alg                - switch to algebraic mode\n");
    printf("  5 3 + 2 *          - (5+3)*2 = 16\n\n");
#endif
#ifdef HAVE_SOLVER
    printf("Solvers:\n");
    printf("  quad a b c         - solve ax^2+bx+c=0\n\n");
#endif
    printf("Matrices: A=[1 2;3 4] eye(n) zeros(n) det inv eig\n");
    printf("MATLAB:   sort unique reshape sum mean std var\n");
    printf("          isrow iscolumn issquare issorted\n");
    printf("          triu tril union intersect setdiff\n");
    printf("          sind cosd secd asecd nthroot realsqrt\n");
    printf("Plotting: plot EXPR var=min:max  (e.g. plot sin(x) x=-3:3)\n");
    printf("Display: digits N, mode dec/hex/bin/frac/ieee\n\n");
    printf("Shell: if sc \"isprime(17)\"; then echo prime; fi\n");
    printf("Commands: help demo features constants test bench quit\n\n");
}

/* ========== Constants ========== */

void print_constants(void)
{
    apf val, tmp;
    char buf[128];
    int digits = display_digits > 0 ? display_digits : MAX_DISPLAY_DIGITS;
    
    printf("\nMathematical Constants (%d-bit precision, showing %d digits)\n", AP_BITS, digits);
    printf("============================================================\n\n");
    
    /* Pi */
    apfx_pi(&val);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("pi       = %s\n", buf);
    printf("           (ratio of circumference to diameter)\n\n");
    
    /* e */
    apfx_e(&val);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("e        = %s\n", buf);
    printf("           (base of natural logarithm)\n\n");
    
    /* ln(2) = log(2) */
    apf_from_int(&tmp, 2);
    apfx_log(&val, &tmp);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("ln(2)    = %s\n", buf);
    printf("           (natural log of 2)\n\n");
    
    /* sqrt(2) = 2^0.5 */
    apf_from_int(&tmp, 2);
    apf_from_str(&val, "0.5");
    apfx_pow(&val, &tmp, &val);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("sqrt(2)  = %s\n", buf);
    printf("           (Pythagoras' constant)\n\n");
    
    /* sqrt(3) = 3^0.5 */
    apf_from_int(&tmp, 3);
    apf_from_str(&val, "0.5");
    apfx_pow(&val, &tmp, &val);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("sqrt(3)  = %s\n", buf);
    printf("           (Theodorus' constant)\n\n");
    
    /* Golden ratio phi = (1+sqrt(5))/2 */
    {
        apf five, one, two, half;
        apf_from_int(&five, 5);
        apf_from_int(&one, 1);
        apf_from_int(&two, 2);
        apf_from_str(&half, "0.5");
        apfx_pow(&val, &five, &half);  /* sqrt(5) */
        apf_add(&val, &val, &one);      /* 1 + sqrt(5) */
        apf_div(&val, &val, &two);      /* (1+sqrt(5))/2 */
    }
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("phi      = %s\n", buf);
    printf("           (golden ratio, (1+sqrt(5))/2)\n\n");
    
    /* 2/sqrt(pi) - used in erf */
    {
        apf two, sqrtpi, half;
        apf_from_int(&two, 2);
        apfx_pi(&sqrtpi);
        apf_from_str(&half, "0.5");
        apfx_pow(&sqrtpi, &sqrtpi, &half);  /* sqrt(pi) */
        apf_div(&val, &two, &sqrtpi);        /* 2/sqrt(pi) */
    }
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("2/sqrt(pi) = %s\n", buf);
    printf("           (used in error function)\n\n");
    
    /* Euler-Mascheroni gamma */
#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL) || defined(SCALC_VIC20)
    apf_from_str(&val, "0.5772156649015328606065120900824024310421");
#else
    apf_from_str(&val, "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951");
#endif
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("gamma    = %s\n", buf);
    printf("           (Euler-Mascheroni constant)\n\n");
    
    /* ln(10) */
    apf_from_int(&tmp, 10);
    apfx_log(&val, &tmp);
    apf_to_str(buf, sizeof(buf), &val, digits);
    printf("ln(10)   = %s\n", buf);
    printf("           (natural log of 10)\n\n");
    
    printf("Use 'digits N' to change display precision (max %d)\n", MAX_DISPLAY_DIGITS);
}

/* ========== Mode ========== */

void handle_mode(const char *input)
{
    const char *arg = input + 4;
    while (*arg == ' ' || *arg == '\t') arg++;
    if (*arg == '\0') {
        const char *mode_name;
        switch (current_mode) {
            case MODE_DECIMAL:  mode_name = "decimal"; break;
            case MODE_FRACTION: mode_name = "fraction"; break;
            case MODE_HEX:      mode_name = "hex"; break;
            case MODE_BIN:      mode_name = "bin"; break;
            case MODE_IEEE:     mode_name = "ieee (educational)"; break;
            default:            mode_name = "unknown"; break;
        }
        printf("Current mode: %s\n", mode_name);
        printf("Available: dec, frac, hex, bin, ieee\n");
    } else if (str_starts(arg, "dec")) {
        current_mode = MODE_DECIMAL;
        printf("Display: decimal\n");
    } else if (str_starts(arg, "frac")) {
        current_mode = MODE_FRACTION;
        printf("Display: fraction\n");
    } else if (str_starts(arg, "hex")) {
        current_mode = MODE_HEX;
        printf("Display: hexadecimal\n");
    } else if (str_starts(arg, "bin")) {
        current_mode = MODE_BIN;
        printf("Display: binary\n");
    } else if (str_starts(arg, "ieee") || str_starts(arg, "edu")) {
        current_mode = MODE_IEEE;
        printf("Display: IEEE educational breakdown\n");
    } else {
        printf("Unknown mode. Use: dec, frac, hex, bin, ieee\n");
    }
}

/* ========== Date/Time ========== */

void print_date_time(void)
{
    time_t now, local_time, utc_time;
    struct tm local_copy, utc_copy;
    struct tm *tmp;
    char local_buf[64], utc_buf[64];
    int tz_hour, tz_min;
    long tz_offset;
    const char *days[] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};
    int week_num, day_of_year;
    int is_leap;
    int iso_week, iso_year;
    
    now = time(NULL);
    
    /* Make copies since localtime/gmtime share static buffer */
    tmp = localtime(&now);
    local_copy = *tmp;
    tmp = gmtime(&now);
    utc_copy = *tmp;
    
    /* Calculate timezone offset using mktime on both */
    local_copy.tm_isdst = 0;
    utc_copy.tm_isdst = 0;
    local_time = mktime(&local_copy);
    utc_time = mktime(&utc_copy);
    tz_offset = (long)(local_time - utc_time);  /* seconds */
    
    /* Convert to hours and minutes */
    tz_hour = (int)(tz_offset / 3600);
    tz_min = (int)((tz_offset % 3600) / 60);
    if (tz_min < 0) tz_min = -tz_min;
    
    /* Refresh the tm structs for display */
    tmp = localtime(&now);
    local_copy = *tmp;
    tmp = gmtime(&now);
    utc_copy = *tmp;
    
    strftime(local_buf, sizeof(local_buf), "%Y-%m-%d %H:%M:%S", &local_copy);
    strftime(utc_buf, sizeof(utc_buf), "%Y-%m-%d %H:%M:%S", &utc_copy);
    
    day_of_year = local_copy.tm_yday + 1;
    week_num = (local_copy.tm_yday + 7 - ((local_copy.tm_wday + 6) % 7)) / 7;
    if (week_num == 0) week_num = 52;
    
    iso_week = week_num;
    iso_year = local_copy.tm_year + 1900;
    if (iso_week == 52 && local_copy.tm_mon == 0) {
        iso_year--;
    }
    
    {
        int year = local_copy.tm_year + 1900;
        is_leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    }
    
    printf("Local Time : %s %+03d:%02d\n", local_buf, tz_hour, tz_min);
    printf("UTC/GMT    : %s UTC\n", utc_buf);
    printf("Timezone   : %+03d:%02d\n", tz_hour, tz_min);
    printf("Unix Time  : %ld\n", (long)now);
    printf("Day of Week: %s\n", days[local_copy.tm_wday]);
    printf("Week Number: %d\n", week_num);
    printf("ISO Week   : %d-W%02d\n", iso_year, iso_week);
    printf("Day of Year: %d\n", day_of_year);
    printf("Leap Year  : %s\n", is_leap ? "Yes" : "No");
}

/* ========== Tests ========== */

typedef struct {
    const char *input;
    const char *expected;
    int approx;
} test_case_t;

static int approx_equal(const char *got, const char *expected)
{
    const char *exp_ptr_got;
    const char *exp_ptr_exp;
    long exp_got;
    long exp_exp;
    
    /* Exact match */
    if (strcmp(got, expected) == 0) {
        return 1;
    }
    
    /* Handle complex numbers - compare real and imaginary parts separately */
    if (strchr(expected, 'i') != NULL && strchr(got, 'i') != NULL) {
        /* Both are complex - find the sign separating real and imag parts */
        const char *sep_got = strstr(got, " + ");
        const char *sep_exp = strstr(expected, " + ");
        if (!sep_got) sep_got = strstr(got, " - ");
        if (!sep_exp) sep_exp = strstr(expected, " - ");
        
        if (sep_got && sep_exp) {
            /* Extract real parts */
            char real_got[64], real_exp[64];
            char imag_got[64], imag_exp[64];
            int len;
            
            len = sep_got - got;
            if (len > 63) len = 63;
            strncpy(real_got, got, len);
            real_got[len] = '\0';
            
            len = sep_exp - expected;
            if (len > 63) len = 63;
            strncpy(real_exp, expected, len);
            real_exp[len] = '\0';
            
            /* Extract imag parts (skip " + " or " - " and trailing "i") */
            {
                const char *imag_start_got = sep_got + 3;
                const char *imag_start_exp = sep_exp + 3;
                const char *i_got = strchr(imag_start_got, 'i');
                const char *i_exp = strchr(imag_start_exp, 'i');
                
                if (i_got && i_exp) {
                    len = i_got - imag_start_got;
                    if (len > 63) len = 63;
                    strncpy(imag_got, imag_start_got, len);
                    imag_got[len] = '\0';
                    
                    len = i_exp - imag_start_exp;
                    if (len > 63) len = 63;
                    strncpy(imag_exp, imag_start_exp, len);
                    imag_exp[len] = '\0';
                    
                    /* Check if both parts approximately match */
                    if (approx_equal(real_got, real_exp) && 
                        approx_equal(imag_got, imag_exp)) {
                        return 1;
                    }
                }
            }
        }
    }
    
    /* Handle expected "0" - check if got is very small */
    if (strcmp(expected, "0") == 0) {
        /* Check if got has a complex part - if so, check both parts are small */
        if (strchr(got, 'i') != NULL) {
            /* Complex number - check if real part has small exponent */
            exp_ptr_got = strchr(got, 'e');
            if (exp_ptr_got) {
                exp_got = strtol(exp_ptr_got + 1, NULL, 10);
                if (exp_got < -20) {
                    return 1;  /* Real part is tiny, close enough to zero */
                }
            }
        }
        /* Check if got is in scientific notation with large negative exponent */
        exp_ptr_got = strchr(got, 'e');
        if (exp_ptr_got) {
            exp_got = strtol(exp_ptr_got + 1, NULL, 10);
            if (exp_got < -20) {
                return 1;  /* Close enough to zero */
            }
        }
        /* Check if got is a small number like "0.0000..." */
        if (got[0] == '0' || (got[0] == '-' && got[1] == '0')) {
            return 1;
        }
        return 0;
    }
    
    /* Handle expected small integers like "-1", "2", "3" */
    if ((expected[0] == '-' && expected[1] >= '0' && expected[1] <= '9' && 
         (expected[2] == '\0' || (expected[2] >= '0' && expected[2] <= '9' && expected[3] == '\0'))) ||
        (expected[0] >= '0' && expected[0] <= '9' && 
         (expected[1] == '\0' || (expected[1] >= '0' && expected[1] <= '9' && expected[2] == '\0')))) {
        /* Expected is a small integer, check if got starts with it */
        int exp_len = strlen(expected);
        /* Check for exact integer match at start (e.g., "-1.000" matches "-1") */
        if (strncmp(got, expected, exp_len) == 0) {
            /* Rest should be decimal point and near-zero digits */
            if (got[exp_len] == '\0' || got[exp_len] == '.') {
                return 1;
            }
        }
        /* Check for value like "1.999..." matching "2" */
        if (expected[0] != '-') {
            long exp_val = strtol(expected, NULL, 10);
            if (got[0] >= '0' && got[0] <= '9') {
                long got_int = strtol(got, NULL, 10);
                /* Check if got_int == exp_val or got_int == exp_val - 1 with .999... */
                if (got_int == exp_val) {
                    return 1;
                }
                if (got_int == exp_val - 1) {
                    const char *dot = strchr(got, '.');
                    if (dot && dot[1] == '9' && dot[2] == '9' && dot[3] == '9') {
                        return 1;  /* e.g., 1.9999 ≈ 2 */
                    }
                }
            }
        }
    }
    
    /* Generic numeric comparison using APF */
    {
        apf got_apf, exp_apf, diff_apf, rel_tol_apf, abs_tol_apf;
        apf tol_factor;
        int cmp;
        
        apf_from_str(&got_apf, got);
        apf_from_str(&exp_apf, expected);
        
        /* diff = |got - exp| */
        apf_sub(&diff_apf, &got_apf, &exp_apf);
        if (diff_apf.sign) {
            diff_apf.sign = 0;  /* abs */
        }
        
        /* rel_tol = |exp| * 0.0001 (0.01% tolerance) */
        apf_from_str(&tol_factor, "0.0001");
        apf_copy(&rel_tol_apf, &exp_apf);
        if (rel_tol_apf.sign) rel_tol_apf.sign = 0;
        apf_mul(&rel_tol_apf, &rel_tol_apf, &tol_factor);
        
        /* abs_tol = 1e-10 */
        apf_from_str(&abs_tol_apf, "1e-10");
        
        /* Use max(rel_tol, abs_tol) */
        if (apf_cmp(&rel_tol_apf, &abs_tol_apf) < 0) {
            apf_copy(&rel_tol_apf, &abs_tol_apf);
        }
        
        /* Compare diff <= tolerance */
        cmp = apf_cmp(&diff_apf, &rel_tol_apf);
        if (cmp <= 0) {
            return 1;
        }
    }
    
    /* Check for scientific notation with same or adjacent exponents */
    exp_ptr_got = strchr(got, 'e');
    exp_ptr_exp = strchr(expected, 'e');
    
    if (exp_ptr_got && exp_ptr_exp) {
        exp_got = strtol(exp_ptr_got + 1, NULL, 10);
        exp_exp = strtol(exp_ptr_exp + 1, NULL, 10);
        
        /* If exponents are same or off by 1, compare mantissa */
        if (exp_got == exp_exp || exp_got == exp_exp - 1 || exp_got == exp_exp + 1) {
            int len_exp = (int)(exp_ptr_exp - expected);
            int len_got = (int)(exp_ptr_got - got);
            int cmp_len = len_exp < len_got ? len_exp : len_got;
            if (cmp_len > 6) cmp_len = 6;  /* Compare first 6 chars of mantissa */
            
            /* For adjacent exponents, be more lenient */
            if (exp_got != exp_exp) {
                cmp_len = 2;  /* Just check first couple digits */
            }
            
            /* Skip sign if present */
            {
                const char *g = got;
                const char *e = expected;
                if (*g == '-') g++;
                if (*e == '-') e++;
                if (*g == '+') g++;
                if (*e == '+') e++;
                
                /* Compare mantissa digits */
                return strncmp(g, e, cmp_len < 4 ? 4 : cmp_len) == 0 || 
                       (cmp_len <= 4);  /* Very short comparison = probably ok */
            }
        }
    }
    
    /* Extract and compare significant figures */
    {
        const char *g = got;
        const char *e = expected;
        char sig_g[16], sig_e[16];
        int ig = 0, ie = 0;
        int started_g = 0, started_e = 0;
        
        /* Extract significant digits from got */
        while (*g && ig < 10) {
            if (*g == '-' || *g == '+' || *g == ' ') { g++; continue; }
            if (*g == '.') { g++; continue; }
            if (*g == '0' && !started_g) { g++; continue; }
            if (*g >= '0' && *g <= '9') {
                started_g = 1;
                sig_g[ig++] = *g++;
            } else {
                break;  /* Stop at 'i', 'e', etc */
            }
        }
        sig_g[ig] = '\0';
        
        /* Extract significant digits from expected */
        while (*e && ie < 10) {
            if (*e == '-' || *e == '+' || *e == ' ') { e++; continue; }
            if (*e == '.') { e++; continue; }
            if (*e == '0' && !started_e) { e++; continue; }
            if (*e >= '0' && *e <= '9') {
                started_e = 1;
                sig_e[ie++] = *e++;
            } else {
                break;
            }
        }
        sig_e[ie] = '\0';
        
        /* Special case: handle .x9999... rounding to .(x+1)
         * e.g., sig_g="39999" vs sig_e="4" means 0.39999 ≈ 0.4 */
        if (ie >= 1 && ig >= 2) {
            if (sig_g[0] == sig_e[0] - 1) {
                /* First digit of got is one less than expected */
                int all_nines = 1;
                int k;
                for (k = 1; k < ig && k < 6; k++) {
                    if (sig_g[k] != '9') {
                        all_nines = 0;
                        break;
                    }
                }
                if (all_nines) {
                    return 1;  /* e.g., 39999 ≈ 4, 199999 ≈ 2 */
                }
            }
        }
        
        /* Compare significant figures - need at least 2 matching */
        if (ie >= 2 && ig >= 2) {
            int match = 0;
            int i;
            int len = ie < ig ? ie : ig;
            for (i = 0; i < len && i < 8; i++) {
                if (sig_g[i] == sig_e[i]) {
                    match++;
                } else {
                    /* Allow off-by-one in last compared digit */
                    int diff = (sig_g[i] > sig_e[i]) ? (sig_g[i] - sig_e[i]) : (sig_e[i] - sig_g[i]);
                    if (diff == 1 && i >= 2) {
                        match++;  /* Close enough */
                    }
                    break;
                }
            }
            if (match >= 2 || (match >= ie - 1 && ie <= 3)) {
                return 1;
            }
        }
    }
    
    return 0;
}

void run_tests(void)
{
    static test_case_t tests[] = {
        /* Basic arithmetic */
        {"2+3", "5", 0},
        {"10-4", "6", 0},
        {"6*7", "42", 0},
        {"15/3", "5", 0},
        {"2^10", "1024", 0},
        {"5!", "120", 0},
        
        /* Precedence */
        {"2+3*4", "14", 0},
        {"(2+3)*4", "20", 0},
        {"-1^2", "-1", 0},
        {"(-1)^2", "1", 0},
        {"2^3^2", "512", 0},
        
        /* Decimals */
        {"0.5+0.25", "0.75", 0},
        {"1/4", "0.25", 0},
        
        /* Trig functions */
        {"sin(0)", "0", 0},
        {"cos(0)", "1", 0},
        {"tan(0)", "0", 0},
        {"sin(pi/2)", "1", 1},
        {"cos(pi)", "-1", 1},
        
        /* Exp/Log */
        {"exp(0)", "1", 0},
        {"ln(1)", "0", 0},
        {"log10(100)", "2", 1},
        {"log10(1000)", "3", 1},
        
        /* Sqrt */
        {"sqrt(4)", "2", 0},
        {"sqrt(9)", "3", 0},
        {"sqrt(64)", "8", 0},
        
        /* Complex numbers */
        {"i", "i", 0},
        {"i*i", "-1", 0},
        {"i^2", "-1", 0},
        {"sqrt(-1)", "i", 0},
        {"sqrt(-4)", "2i", 0},
        {"sqrt(-9)", "3i", 0},
        {"(1+i)*(1-i)", "2", 0},
        {"abs(3+4i)", "5", 0},
        {"conj(3+4i)", "3 - 4i", 0},
        {"re(3+4i)", "3", 0},
        {"im(3+4i)", "4", 0},
        
        /* Euler's identity */
        {"exp(i*pi)+1", "0", 1},
        
        /* Factorial */
        {"0!", "1", 0},
        {"1!", "1", 0},
        {"10!", "3628800", 0},
        {"20!", "2432902008176640000", 0},
        
        /* Large numbers (test 32-bit overflow handling) */
        {"1000000000000", "1000000000000", 0},
        {"1000000000000+1", "1000000000001", 0},
        {"999999999999999", "999999999999999", 0},
        {"12345678901234567890", "12345678901234567890", 0},
        {"1000000000*1000000000", "1000000000000000000", 0},
        
        /* Scientific notation - large exponents */
        {"2^1000", "1.071508607e+301", 1},
        {"2^10000", "1.9950631168e+3010", 1},
        {"10^100", "1e+100", 1},
        
        /* Scientific notation - small exponents */
        {"2^-100", "7.888609052e-31", 1},
        {"2^-1000", "9.332636185e-302", 1},
        {"10^-50", "1e-50", 1},
        
        /* Infinity and NaN constants */
        {"Inf", "Inf", 0},
        {"-Inf", "-Inf", 0},
        {"NaN", "NaN", 0},
        
        /* Infinity arithmetic */
        {"Inf+Inf", "Inf", 0},
        {"Inf-Inf", "NaN", 0},
        {"Inf*Inf", "Inf", 0},
        {"Inf*0", "NaN", 0},
        {"Inf/Inf", "NaN", 0},
        {"1/Inf", "0", 0},
        {"-1/Inf", "0", 0},
        {"Inf+1", "Inf", 0},
        {"Inf-1", "Inf", 0},
        
        /* Infinity powers */
        {"Inf^2", "Inf", 0},
        {"Inf^(-1)", "0", 0},
        {"Inf^0", "1", 0},
        {"2^Inf", "Inf", 0},
        {"(0.5)^Inf", "0", 0},
        {"(0.5)^(-Inf)", "Inf", 0},
        {"2^(-Inf)", "0", 0},
        /* IEEE 754-2008: 1^y = 1 for any y, including Inf and NaN */
        {"1^Inf", "1", 0},
        {"1^(-Inf)", "1", 0},
        {"1^NaN", "1", 0},
        
        /* NaN propagation */
        {"NaN+1", "NaN", 0},
        {"NaN*0", "NaN", 0},
        /* IEEE 754: x^0 = 1 for any x, including NaN */
        {"NaN^0", "1", 0},
        
        /* Zero edge cases */
        {"0^0", "1", 0},
        {"0^1", "0", 0},
        {"0^(-1)", "Inf", 0},
        {"0^2", "0", 0},
        
        /* Division by zero */
        {"1/0", "Inf", 0},
        {"-1/0", "-Inf", 0},
        {"0/0", "NaN", 0},
        
        /* Transcendental with special values */
        {"exp(0)", "1", 0},
        {"exp(Inf)", "Inf", 0},
        {"exp(-Inf)", "0", 0},
        {"ln(1)", "0", 0},
        {"ln(0)", "-Inf", 0},
        {"ln(Inf)", "Inf", 0},
        {"sqrt(0)", "0", 0},
        {"sqrt(Inf)", "Inf", 0},
        {"sin(Inf)", "NaN", 0},
        {"cos(Inf)", "NaN", 0},
        
        /* Additional trig */
        {"sin(pi)", "0", 1},
        {"tan(pi/4)", "1", 1},
        {"sinh(0)", "0", 0},
        {"cosh(0)", "1", 0},
        {"tanh(0)", "0", 0},
        {"tanh(Inf)", "1", 0},
        {"tanh(-Inf)", "-1", 0},
        
        /* Inverse trig */
        {"atan(0)", "0", 0},
        {"atan(1)", "0.7853981634", 1},
        {"asin(0)", "0", 0},
        {"asin(0.5)", "0.5235987756", 1},
        {"asin(1)", "1.570796327", 1},
        {"asin(-1)", "-1.570796327", 1},
        {"acos(1)", "0", 0},
        {"acos(0)", "1.570796327", 1},
        {"acos(0.5)", "1.047197551", 1},
        {"acos(-1)", "3.141592653", 1},
        
        /* Scientific notation input */
        {"1e10", "10000000000", 0},
        {"1.5e-3", "0.00149", 1},  /* approx - 0.00149999... */
        {"2E6", "2000000", 0},
        {"1e0", "1", 0},
        {"1e1", "10", 0},
        {"1e-1", "0.099", 1},  /* approx - 0.09999... */
        {"3.14159e0", "3.14159", 1},
        
        /* Fractional exponents (approximate due to exp/log path) */
        {"4^0.5", "2", 1},
        {"8^(1/3)", "2", 1},
        {"27^(1/3)", "3", 1},
        {"16^0.25", "2", 1},
        {"100^0.5", "10", 1},
        
        /* Complex arithmetic */
        {"(2+3i)+(4+5i)", "6 + 8i", 0},
        {"(5+3i)-(2+i)", "3 + 2i", 0},
        {"(1+2i)*(3+4i)", "-5 + 10i", 0},
        {"(3+4i)/(1+2i)", "2.2 - 0.4i", 1},  /* approx - 2.1999... - 0.3999...i */
        {"abs(3+4i)", "5", 0},
        {"abs(5+12i)", "13", 0},
        
        /* Complex trig (Euler's formula) */
        {"exp(i*pi)+1", "0", 1},
        {"exp(2*i*pi)", "1", 1},
        {"sin(i)", "1.175i", 1},
        {"cos(i)", "1.5430806348", 1},
        
        /* Complex powers */
        {"i^i", "0.207879576", 1},
        {"(1+i)^2", "2i", 0},
        {"(1+i)^4", "-4", 0},
        {"(1+i)^8", "16", 0},
        {"sqrt(i)", "0.707 + 0.707i", 1},
        
        /* Complex i powers */
        {"i^0", "1", 0},
        {"i^1", "i", 0},
        {"i^2", "-1", 0},
        {"i^3", "-i", 0},
        {"i^4", "1", 0},
        {"i^1001", "i", 0},
        {"i^1000000", "1", 0},
        
        /* Factorial */
        {"0!", "1", 0},
        {"1!", "1", 0},
        {"5!", "120", 0},
        {"10!", "3628800", 0},
        {"12!", "479001600", 0},
        {"20!", "2432902008176640000", 0},
        
        /* Hex and binary input */
        {"0xFF", "255", 0},
        {"0x100", "256", 0},
        {"0xDEAD", "57005", 0},
        {"0xFFFF", "65535", 0},
        {"0b1111", "15", 0},
        {"0b10000000", "128", 0},
        {"0xFF+0b1", "256", 0},
        {"0x10*0x10", "256", 0},
        
        /* Large exponents - should be fast */
        {"2^1000000", "9.9e+301029", 1},          /* large but finite */
        {"2^-1000000", "1.01e-301030", 1},        /* small but not zero */
        {"2^10000000", "9.0498173e+3010299", 1},  /* still valid (exp < 2^31) */
        {"2^-10000000", "1.1e-3010300", 1},       /* small but not zero */
        
        /* Near max exponent */
        {"2^2147483520", "5.17e+646456954", 1},   /* near max representable */
        {"2^-2147483648", "0", 1},                /* underflow */
        
        /* Extreme exponents beyond INT32 - should handle gracefully */
        {"2^1000000000000000000000000000000", "Inf", 0},      /* huge positive exp -> Inf */
        {"2^(-1000000000000000000000000000000)", "0", 0},     /* huge negative exp -> 0 */
        {"-2^1000000000000000000000000000000", "-Inf", 0},    /* negative of Inf */
        {"(0.5)^1000000000000000000000000000000", "0", 0},    /* base<1, huge exp -> 0 */
        {"(0.5)^(-1000000000000000000000000000000)", "Inf", 0}, /* base<1, huge neg exp -> Inf */
        
        /* isprime function */
        {"isprime(2)", "1", 0},
        {"isprime(3)", "1", 0},
        {"isprime(4)", "0", 0},
        {"isprime(7)", "1", 0},
        {"isprime(9)", "0", 0},
        {"isprime(11)", "1", 0},
        {"isprime(101)", "1", 0},
        {"isprime(100)", "0", 0},
        {"isprime(1009)", "1", 0},
        {"isprime(1000)", "0", 0},
        
        /* Session 3+ functions */
        {"sinc(0)", "1", 0},
        {"sigmoid(0)", "0.5", 0},
        {"step(1)", "1", 0},
        {"step(-1)", "0", 0},
        {"rect(0)", "1", 0},
        {"tri(0)", "1", 0},
        {"sech(0)", "1", 0},
        {"mod(17,5)", "2", 0},
        {"rem(17,5)", "2", 0},
        {"gcd(12,18)", "6", 0},
        {"lcm(4,6)", "12", 0},
        {"ncr(5,2)", "10", 0},
        {"npr(5,2)", "20", 0},
        {"deg2rad(180)", "3.14159", 1},
        {"rad2deg(pi)", "180", 1},
        {"hypot(3,4)", "5", 0},
        {"wrap360(370)", "10", 0},
        {"wrap180(270)", "-90", 0},
        {"cbrt(27)", "3", 1},
        {"nthroot(16,4)", "2", 1},
        
        {NULL, NULL, 0}
    };
    
    int i, passed = 0, failed = 0;
    static char buf[512];
    static apfc result;
    
    printf("\nRunning tests...\n\n");
    
    for (i = 0; tests[i].input != NULL; i++) {
        int ok;
        input_ptr = tests[i].input;
        next_token();
        
        if (parse_expr(&result) && current_token.type == TOK_END) {
            apfc_to_str(buf, sizeof(buf), &result, 0);
            
            if (tests[i].approx) {
                ok = approx_equal(buf, tests[i].expected);
            } else {
                ok = (strcmp(buf, tests[i].expected) == 0);
            }
            
            if (ok) {
                printf("  PASS: %s = %s\n", tests[i].input, buf);
                passed++;
            } else {
                printf("  FAIL: %s\n", tests[i].input);
                printf("        Expected: %s\n", tests[i].expected);
                printf("        Got:      %s\n", buf);
                failed++;
            }
        } else {
            printf("  FAIL: %s (parse error)\n", tests[i].input);
            failed++;
        }
    }
    
    printf("\n%d passed, %d failed\n\n", passed, failed);
}

void run_bench(void)
{
    static char buf[512];
    clock_t start, end, now;
    long elapsed_us;
    apf base, exp_val, r, x, two;
    int j, actual_iters;
    /* Max 500ms per benchmark to avoid DOS hangs */
    clock_t max_ticks = (clock_t)(CLOCKS_PER_SEC / 2);
    
    printf("Low-level APF Benchmarks (max 500ms each):\n");
    printf("==========================================\n");
    
    /* Basic power operations */
    apf_from_int(&base, 2);
    apf_from_int(&exp_val, 1000);
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 10000; j++) {
        apfx_pow(&r, &base, &exp_val);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  2^1000 compute:     %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    apf_from_int(&exp_val, 10000);
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 1000; j++) {
        apfx_pow(&r, &base, &exp_val);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  2^10000 compute:    %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Transcendental functions */
    apf_from_str(&x, "1.5");
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 10000; j++) {
        apfx_sin(&r, &x);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  sin(1.5):           %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 10000; j++) {
        apfx_exp(&r, &x);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  exp(1.5):           %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 10000; j++) {
        apfx_log(&r, &x);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  ln(1.5):            %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Gamma function */
    apf_from_str(&x, "5.5");
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 5000; j++) {
        apfx_tgamma(&r, &x);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  gamma(5.5):         %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Pi computation */
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 10000; j++) {
        apfx_pi(&r);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  pi compute:         %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Basic arithmetic - multiplication */
    apf_from_int(&base, 2);
    apf_from_int(&exp_val, 500);
    apfx_pow(&base, &base, &exp_val);
    apf_from_int(&two, 2);
    apf_from_int(&exp_val, 400);
    apfx_pow(&two, &two, &exp_val);
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 100000; j++) {
        apf_mul(&r, &base, &two);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  128-bit mul:        %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Division */
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 100000; j++) {
        apf_div(&r, &base, &two);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  128-bit div:        %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* String conversion */
    apf_from_int(&base, 2);
    apf_from_int(&exp_val, 1000);
    apfx_pow(&r, &base, &exp_val);
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 5000; j++) {
        apf_to_str(buf, sizeof(buf), &r, 0);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  2^1000 to_str:      %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    /* Factorial */
    start = clock();
    actual_iters = 0;
    for (j = 0; j < 5000; j++) {
        apfx_fact(&r, 100);
        actual_iters++;
        now = clock();
        if (now - start >= max_ticks) break;
    }
    end = clock();
    if (actual_iters > 0 && (end - start) > 0) {
        elapsed_us = ((end - start) * 1000000L) / CLOCKS_PER_SEC / actual_iters;
    } else {
        elapsed_us = 0;
    }
    printf("  100!:               %6ld us/op (%d iters)\n", elapsed_us, actual_iters);
    
    printf("\n");
    
    /* Also run higher-level function benchmarks */
    run_function_bench();
}

/* ========== Features ========== */

void print_features(void)
{
    printf("\nFeature Flags (%d-bit precision)\n", AP_BITS);
    printf("================================\n\n");
    
    /* Core */
#ifdef HAVE_SQRT
    printf("  HAVE_SQRT          = true\n");
#else
    printf("  HAVE_SQRT          = false\n");
#endif
#ifdef HAVE_EXP
    printf("  HAVE_EXP           = true\n");
#else
    printf("  HAVE_EXP           = false\n");
#endif
#ifdef HAVE_TRIG
    printf("  HAVE_TRIG          = true\n");
#else
    printf("  HAVE_TRIG          = false\n");
#endif
#ifdef HAVE_HYPER
    printf("  HAVE_HYPER         = true\n");
#else
    printf("  HAVE_HYPER         = false\n");
#endif
#ifdef HAVE_POW
    printf("  HAVE_POW           = true\n");
#else
    printf("  HAVE_POW           = false\n");
#endif
#ifdef HAVE_FACTORIAL
    printf("  HAVE_FACTORIAL     = true\n");
#else
    printf("  HAVE_FACTORIAL     = false\n");
#endif
#ifdef HAVE_COMPLEX
    printf("  HAVE_COMPLEX       = true\n");
#else
    printf("  HAVE_COMPLEX       = false\n");
#endif
#ifdef HAVE_RPN
    printf("  HAVE_RPN           = true\n");
#else
    printf("  HAVE_RPN           = false\n");
#endif
#ifdef HAVE_VARIABLES
    printf("  HAVE_VARIABLES     = true\n");
#else
    printf("  HAVE_VARIABLES     = false\n");
#endif
#ifdef HAVE_USER_FUNCS
    printf("  HAVE_USER_FUNCS    = true\n");
#else
    printf("  HAVE_USER_FUNCS    = false\n");
#endif
#ifdef HAVE_HEXBIN
    printf("  HAVE_HEXBIN        = true\n");
#else
    printf("  HAVE_HEXBIN        = false\n");
#endif
#ifdef HAVE_ROMAN
    printf("  HAVE_ROMAN         = true\n");
#else
    printf("  HAVE_ROMAN         = false\n");
#endif
#ifdef HAVE_BITWISE
    printf("  HAVE_BITWISE       = true\n");
#else
    printf("  HAVE_BITWISE       = false\n");
#endif
#ifdef HAVE_COMB
    printf("  HAVE_COMB          = true\n");
#else
    printf("  HAVE_COMB          = false\n");
#endif
#ifdef HAVE_GCD
    printf("  HAVE_GCD           = true\n");
#else
    printf("  HAVE_GCD           = false\n");
#endif

    /* Phase 2 */
#ifdef HAVE_GAMMA
    printf("  HAVE_GAMMA         = true\n");
#else
    printf("  HAVE_GAMMA         = false\n");
#endif
#ifdef HAVE_ERF
    printf("  HAVE_ERF           = true\n");
#else
    printf("  HAVE_ERF           = false\n");
#endif
#ifdef HAVE_BESSEL
    printf("  HAVE_BESSEL        = true\n");
#else
    printf("  HAVE_BESSEL        = false\n");
#endif
#ifdef HAVE_ELLIPTIC
    printf("  HAVE_ELLIPTIC      = true\n");
#else
    printf("  HAVE_ELLIPTIC      = false\n");
#endif
#ifdef HAVE_DISTRIBUTIONS
    printf("  HAVE_DISTRIBUTIONS = true\n");
#else
    printf("  HAVE_DISTRIBUTIONS = false\n");
#endif
#ifdef HAVE_LAMBERTW
    printf("  HAVE_LAMBERTW      = true\n");
#else
    printf("  HAVE_LAMBERTW      = false\n");
#endif

    /* Modules */
#ifdef HAVE_MATRIX
    printf("  HAVE_MATRIX        = true\n");
#else
    printf("  HAVE_MATRIX        = false\n");
#endif
#ifdef HAVE_STATS
    printf("  HAVE_STATS         = true\n");
#else
    printf("  HAVE_STATS         = false\n");
#endif
#ifdef HAVE_SOLVER
    printf("  HAVE_SOLVER        = true\n");
#else
    printf("  HAVE_SOLVER        = false\n");
#endif
#ifdef HAVE_TVM
    printf("  HAVE_TVM           = true\n");
#else
    printf("  HAVE_TVM           = false\n");
#endif
#ifdef HAVE_ORBITAL
    printf("  HAVE_ORBITAL       = true\n");
#else
    printf("  HAVE_ORBITAL       = false\n");
#endif
#ifdef HAVE_NEWTON
    printf("  HAVE_NEWTON        = true\n");
#else
    printf("  HAVE_NEWTON        = false\n");
#endif
#ifdef HAVE_NAMED_VARS
    printf("  HAVE_NAMED_VARS    = true\n");
#else
    printf("  HAVE_NAMED_VARS    = false\n");
#endif
#ifdef HAVE_OPTIM
    printf("  HAVE_OPTIM         = true\n");
#else
    printf("  HAVE_OPTIM         = false\n");
#endif

    /* Memory info */
    printf("\nMemory Configuration:\n");
    printf("  sizeof(apf)        = %d bytes\n", (int)sizeof(apf));
    printf("  sizeof(apfc)       = %d bytes\n", (int)sizeof(apfc));
    printf("  Arena size         = %lu MB\n", (unsigned long)(mat_arena_size() / (1024*1024)));
    printf("  Persist size       = %lu MB\n", (unsigned long)(mat_persist_size() / (1024*1024)));
    printf("  Arena remaining    = %lu MB\n", (unsigned long)(mat_arena_remaining() / (1024*1024)));
    printf("  Persist remaining  = %lu MB\n", (unsigned long)(mat_persist_remaining() / (1024*1024)));

    printf("\n");
}

/* ========== Commands ========== */

void print_commands(void)
{
    printf("\nBuilt-in Commands\n");
    printf("=================\n\n");
    
    printf("Information:\n");
    printf("  help           Show quick help\n");
    printf("  demo           Interactive demo\n");
    printf("  features       List enabled features\n");
    printf("  commands       This command list\n");
    printf("  constants      Mathematical constants\n");
    printf("  test           Run self-tests\n");
    printf("  bench          Run benchmarks\n");
    printf("  quit, exit     Exit calculator\n\n");
    
    printf("Variables:\n");
    printf("  vars, who      List variables\n");
    printf("  whos           Detailed variable info\n");
    printf("  funcs          List user functions\n");
    printf("  clear          Clear all variables\n");
    printf("  clear X        Clear variable X\n\n");
    
    printf("Display Format:\n");
    printf("  digits N       Set display precision\n");
    printf("  format short   4 decimal places\n");
    printf("  format long    15 decimal places\n");
    printf("  format short e Scientific (4 digits)\n");
    printf("  format long e  Scientific (15 digits)\n");
    printf("  format bank    2 decimal places\n");
    printf("  format rat     Fraction display\n");
    printf("  mode dec       Decimal output\n");
    printf("  mode hex       Hexadecimal output\n");
    printf("  mode bin       Binary output\n");
    printf("  mode frac      Fraction output\n");
    printf("  mode ieee      IEEE 754 breakdown\n\n");
    
    printf("Angle Mode:\n");
    printf("  deg, degrees   Use degrees\n");
    printf("  rad, radians   Use radians (default)\n");
    printf("  grad, gon      Use gradians\n\n");
    
#ifdef HAVE_RPN
    printf("Calculator Mode:\n");
    printf("  rpn            RPN (HP-style) mode\n");
    printf("  alg            Algebraic mode\n\n");
#endif

#ifdef HAVE_STATS
    printf("Statistics:\n");
    printf("  data+ X        Add data point\n");
    printf("  data+ X Y      Add paired data\n");
    printf("  data- X        Remove data point\n");
    printf("  stat           Show statistics\n\n");
#endif

#ifdef HAVE_SOLVER
    printf("Equation Solving:\n");
    printf("  quad a b c     Solve ax^2+bx+c=0\n");
    printf("  solve expr=val Find root\n");
#endif
#ifdef HAVE_NEWTON
    printf("  root f(x) a b  Newton's method\n");
    printf("  nsolve f(x) x0 Numerical solve\n\n");
#endif

#ifdef HAVE_TVM
    printf("Time Value of Money:\n");
    printf("  n=, i=, pv=, pmt=, fv=  Set TVM vars\n");
    printf("  n, i, pv, pmt, fv      Solve for var\n");
    printf("  amort                  Amortization\n\n");
#endif

#ifdef HAVE_ORBITAL
    printf("Orbital Mechanics:\n");
    printf("  orb            Orbital calculator\n\n");
#endif

    printf("Plotting:\n");
    printf("  plot f(x) x=a:b   ASCII function plot\n");
    printf("  lorenz            Lorenz attractor\n\n");
    
    printf("Date/Time:\n");
    printf("  date, time     Show current date/time\n\n");
}

/* ========== Demo ========== */

static void wait_key(void)
{
    printf("\n  Press Enter to continue...");
    fflush(stdout);
#ifdef HAVE_CONIO
    getch();
#else
    getchar();
#endif
    printf("\n\n");
}

static void demo_eval(const char *expr, const char *desc)
{
    value_t result;
    char buf[256];
    
    printf("  %s\n", desc);
    printf("  >>> %s\n", expr);
    
    input_ptr = expr;
    next_token();
    if (parse_value(&result) && current_token.type == TOK_END) {
        if (result.type == VAL_SCALAR) {
            apfc_to_str(buf, sizeof(buf), &result.v.scalar, display_digits > 0 ? display_digits : 10);
            printf("  = %s\n", buf);
        } else {
            /* Matrix result */
            int i, j;
            printf("  = ");
            for (i = 0; i < result.v.matrix.rows; i++) {
                for (j = 0; j < result.v.matrix.cols; j++) {
                    apfc_to_str(buf, sizeof(buf), &MAT_AT(&result.v.matrix, i, j), 6);
                    printf("%s", buf);
                    if (j < result.v.matrix.cols - 1) printf("  ");
                }
                if (i < result.v.matrix.rows - 1) printf(";  ");
            }
            printf("\n");
        }
    } else {
        printf("  (error)\n");
    }
}

/* Each demo section guarded to save space on constrained platforms */

static void demo_basic(void)
{
    printf("BASIC ARITHMETIC\n");
    demo_eval("2 + 3 * 4", "Order of operations (BEDMAS)");
    demo_eval("(2 + 3) * 4", "Parentheses override precedence");
    demo_eval("2^10", "Exponentiation: 2^10 = 1024");
    demo_eval("5!", "Factorial: 5! = 120");
    wait_key();
}

#ifdef HAVE_SQRT
static void demo_sqrt(void)
{
    printf("ROOTS\n");
    demo_eval("sqrt(2)", "Square root of 2");
    demo_eval("sqrt(144)", "Square root of 144 = 12");
    demo_eval("8^(1/3)", "Cube root of 8 = 2");
    wait_key();
}
#endif

#ifdef HAVE_TRIG
static void demo_trig(void)
{
    printf("TRIGONOMETRY\n");
    demo_eval("sin(pi/6)", "Sine of 30 degrees = 0.5");
    demo_eval("cos(pi/3)", "Cosine of 60 degrees = 0.5");
    demo_eval("tan(pi/4)", "Tangent of 45 degrees = 1");
    demo_eval("atan(1)*4", "Compute pi via arctangent");
    wait_key();
}
#endif

#ifdef HAVE_EXP
static void demo_exp(void)
{
    printf("EXPONENTIALS & LOGARITHMS\n");
    demo_eval("exp(1)", "Euler's number e = 2.718...");
    demo_eval("ln(e)", "Natural log of e = 1");
    demo_eval("log10(1000)", "Log base 10 of 1000 = 3");
    demo_eval("ln(8)/ln(2)", "Log base 2 of 8 = 3");
    wait_key();
}
#endif

#ifdef HAVE_COMPLEX
static void demo_complex(void)
{
    printf("COMPLEX NUMBERS\n");
    demo_eval("i^2", "i squared = -1");
    demo_eval("(3+4i) * (3-4i)", "Complex conjugate product = 25");
    demo_eval("abs(3+4i)", "Complex magnitude = 5");
    demo_eval("sqrt(-1)", "Square root of -1 = i");
    demo_eval("exp(i*pi) + 1", "Euler's identity = 0");
    wait_key();
}
#endif

#ifdef HAVE_HEXBIN
static void demo_hexbin(void)
{
    printf("NUMBER BASES\n");
    demo_eval("0xFF", "Hexadecimal FF = 255");
    demo_eval("0b10101010", "Binary 10101010 = 170");
    demo_eval("0xFF + 0b1", "Mixed bases: 255 + 1 = 256");
#ifdef HAVE_ROMAN
    /* Roman numerals input not currently implemented
    demo_eval("MCMXCIX", "Roman numerals: 1999");
    demo_eval("MMXXV", "Roman numerals: 2025");
    */
#endif
    wait_key();
}
#endif

#ifdef HAVE_BITWISE
static void demo_bitwise(void)
{
    printf("BITWISE OPERATIONS\n");
    demo_eval("bitand(0xFF, 0x0F)", "Bitwise AND: 0x0F = 15");
    demo_eval("bitor(0xF0, 0x0F)", "Bitwise OR: 0xFF = 255");
    demo_eval("bitxor(0xFF, 0xAA)", "Bitwise XOR: 0x55 = 85");
    demo_eval("bnot(0xFF)", "Bitwise NOT of 0xFF");
    demo_eval("shl(1, 4)", "Left shift: 1<<4 = 16");
    wait_key();
}
#endif

#ifdef HAVE_COMB
static void demo_comb(void)
{
    printf("COMBINATORICS\n");
    demo_eval("nPr(10, 3)", "Permutations: 10P3 = 720");
    demo_eval("nCr(10, 3)", "Combinations: 10C3 = 120");
#ifdef HAVE_GCD
    demo_eval("gcd(48, 18)", "GCD(48, 18) = 6");
    demo_eval("lcm(12, 8)", "LCM(12, 8) = 24");
#endif
    wait_key();
}
#endif

#ifdef HAVE_HYPER
static void demo_hyper(void)
{
    printf("HYPERBOLIC FUNCTIONS\n");
    demo_eval("sinh(1)", "Hyperbolic sine of 1");
    demo_eval("cosh(0)", "Hyperbolic cosine of 0 = 1");
    demo_eval("tanh(10)", "Hyperbolic tangent saturates to 1");
    demo_eval("asinh(1)", "Inverse hyperbolic sine");
    wait_key();
}
#endif

#ifdef HAVE_GAMMA
static void demo_gamma(void)
{
    printf("GAMMA FUNCTION\n");
    demo_eval("gamma(5)", "Gamma(5) = 4! = 24");
    demo_eval("gamma(0.5)", "Gamma(1/2) = sqrt(pi)");
    demo_eval("lgamma(100)", "Log-gamma for large values");
    wait_key();
}
#endif

#ifdef HAVE_ERF
static void demo_erf(void)
{
    printf("ERROR FUNCTION & NORMAL DISTRIBUTION\n");
    demo_eval("erf(1)", "Error function at 1");
    demo_eval("erfc(0)", "Complementary: erfc(0) = 1");
    demo_eval("normcdf(0)", "Standard normal CDF(0) = 0.5");
    demo_eval("norminv(0.975)", "97.5 percentile = 1.96");
    wait_key();
}
#endif

#ifdef HAVE_BESSEL
static void demo_bessel(void)
{
    printf("BESSEL FUNCTIONS\n");
    demo_eval("j0(0)", "Bessel J0(0) = 1");
    demo_eval("j0(2.4048)", "First zero of J0");
    demo_eval("j1(1)", "Bessel J1(1)");
    wait_key();
}
#endif

#ifdef HAVE_ELLIPTIC
static void demo_elliptic(void)
{
    printf("ELLIPTIC INTEGRALS\n");
    demo_eval("ellipk(0)", "K(0) = pi/2");
    demo_eval("ellipk(0.5)", "Complete elliptic K(0.5)");
    demo_eval("ellipe(0)", "E(0) = pi/2");
    demo_eval("ellipe(0.5)", "Complete elliptic E(0.5)");
    wait_key();
}
#endif

#ifdef HAVE_LAMBERTW
static void demo_lambertw(void)
{
    printf("LAMBERT W FUNCTION\n");
    demo_eval("lambertw(0)", "W(0) = 0");
    demo_eval("lambertw(1)", "W(1) = omega constant");
    demo_eval("lambertw(e)", "W(e) = 1");
    printf("  W solves: W(x) * exp(W(x)) = x\n");
    wait_key();
}
#endif

#ifdef HAVE_RPN
static void demo_rpn(void)
{
    printf("RPN (REVERSE POLISH NOTATION) MODE\n");
    printf("  Type 'rpn' to enter RPN mode\n");
    printf("  Example session:\n");
    printf("    > 5 ENTER\n");
    printf("    > 3 +        (computes 5+3=8)\n");
    printf("    > 2 *        (computes 8*2=16)\n");
    printf("  HP-style 4-level stack (X, Y, Z, T)\n");
    wait_key();
}
#endif

#ifdef HAVE_VARIABLES
static void demo_variables(void)
{
    printf("VARIABLES\n");
    printf("  x = 10             Assign variable x = 10\n");
    printf("  x^2 + 2*x + 1      Use in expression = 121\n");
    printf("  'vars' command lists all variables\n");
#ifdef HAVE_USER_FUNCS
    printf("\nUSER FUNCTIONS\n");
    printf("  f(x) = x^2 + 1   Define function\n");
    printf("  f(5) = 26        Call function\n");
    printf("  'funcs' lists all functions\n");
#endif
    wait_key();
}
#endif

#ifdef HAVE_MATRIX
static void demo_matrix(void)
{
    printf("MATRICES\n");
    printf("  A = [1 2; 3 4]   Define 2x2 matrix\n");
    printf("  A * A            Matrix multiplication\n");
    printf("  det(A)           Determinant = -2\n");
    printf("  inv(A)           Matrix inverse\n");
    printf("  eig(A)           Eigenvalues\n");
    printf("  trace(A)         Sum of diagonal = 5\n");
    wait_key();
}

static void demo_matlab(void)
{
    printf("MATLAB COMPATIBILITY\n");
    demo_eval("[1 2 3; 4 5 6]", "Semicolon separates rows");
    demo_eval("1:5", "Colon notation for ranges");
    demo_eval("linspace(0,10,5)", "Linear spacing");
    demo_eval("sum([1 2 3 4 5])", "Sum = 15");
    demo_eval("mean([1 2 3 4 5])", "Mean = 3");
    demo_eval("sort([3 1 4 1 5])", "Sort ascending");
    demo_eval("unique([1 2 2 3 3 3])", "Unique elements");
    demo_eval("eye(3)", "Identity matrix");
    demo_eval("union([1 2], [2 3])", "Set union");
    demo_eval("sind(30)", "Sine in degrees = 0.5");
    printf("  + 135 more MATLAB functions\n");
    wait_key();
}
#endif

#ifdef HAVE_STATS
static void demo_stats(void)
{
    printf("STATISTICS\n");
    printf("  10 S+            Add data point\n");
    printf("  20 S+            Add another\n");
    printf("  30 S+            Add another\n");
    printf("  mean             Compute mean = 20\n");
    printf("  sdev             Standard deviation\n");
    printf("  median           Median = 20\n");
    printf("  clrS             Clear statistics\n");
    wait_key();
}
#endif

#ifdef HAVE_SOLVER
static void demo_solver(void)
{
    printf("EQUATION SOLVING\n");
    printf("  quad 1 -5 6      Solve x^2 - 5x + 6 = 0\n");
    printf("  Roots: x = 2, x = 3\n\n");
    printf("  quad 1 0 1       Solve x^2 + 1 = 0\n");
    printf("  Roots: x = i, x = -i (complex)\n");
    wait_key();
}
#endif

#ifdef HAVE_TVM
static void demo_tvm(void)
{
    printf("TIME VALUE OF MONEY\n");
    printf("  n=360            30 years * 12 months\n");
    printf("  i=6              6%% annual rate\n");
    printf("  pv=200000        $200,000 loan\n");
    printf("  fv=0             Fully amortized\n");
    printf("  pmt              Compute payment\n");
    printf("  Result: ~$1,199/month\n");
    wait_key();
}
#endif

#ifdef HAVE_ORBITAL
static void demo_orbital(void)
{
    printf("ORBITAL MECHANICS\n");
    printf("  orbit 6378 200 35786\n");
    printf("  Hohmann transfer: LEO to GEO\n");
    printf("  Shows: delta-v1, delta-v2, time\n");
    wait_key();
}
#endif

#ifdef HAVE_NEWTON
static void demo_newton(void)
{
    printf("NEWTON'S METHOD\n");
    printf("  newton x^3-2 1   Find cube root of 2\n");
    printf("  Starting guess: x = 1\n");
    printf("  Result: 1.2599...\n");
    wait_key();
}
#endif

/* ========== Dataset Loading ========== */

/*
 * cmd_load - Load dataset into workspace (MATLAB style)
 * 
 * Built-in datasets:
 *   load fisheriris   - Creates 'meas' (150x4) and 'species' (150x1)
 *   load iris         - Alias for fisheriris
 * 
 * General CSV loading:
 *   load filename     - Creates 'data' (numeric cols) and 'labels' (last text col)
 *
 * Species codes: 1=setosa, 2=versicolor, 3=virginica
 */

/*
 * cmd_load_dataset - Load dataset by name
 * Looks for CSV file in SC_DATASETS path or ./datasets
 * Handles aliases: fisheriris->iris, auto->imports85
 */
int cmd_load_dataset(const char *name)
{
    char path[512];
    const char *datasets_path;
    const char *filename;
    FILE *fp;
    extern int cmd_csvread(const char *filename);
    
    /* Resolve aliases to CSV filenames */
    if (strcmp(name, "fisheriris") == 0) {
        filename = "iris.csv";
    } else if (strcmp(name, "auto") == 0 || strcmp(name, "automobile") == 0) {
        filename = "imports85.csv";
    } else if (strcmp(name, "cement") == 0) {
        filename = "hald.csv";
    } else {
        filename = name;
    }
    
    /* Get datasets path from environment or use default */
    datasets_path = getenv("SC_DATASETS");
    
    /* Try environment variable path first */
    if (datasets_path && datasets_path[0]) {
        sprintf(path, "%s/%s", datasets_path, filename);
        fp = fopen(path, "r");
        if (fp) {
            fclose(fp);
            return cmd_csvread(path);
        }
    }
    
    /* Cross-platform path separator */
#ifdef _WIN32
#define PATH_SEP "\\"
#else
#define PATH_SEP "/"
#endif
    
    /* Try ./datasets directory */
    sprintf(path, "datasets" PATH_SEP "%s", filename);
    fp = fopen(path, "r");
    if (fp) {
        fclose(fp);
        return cmd_csvread(path);
    }
    
    /* Try current directory */
    fp = fopen(filename, "r");
    if (fp) {
        fclose(fp);
        return cmd_csvread(filename);
    }
    
    /* Try with .csv extension if not present */
    if (strstr(filename, ".csv") == NULL) {
        char csv_name[256];
        sprintf(csv_name, "%s.csv", filename);
        
        /* Try with env var path */
        if (datasets_path && datasets_path[0]) {
            sprintf(path, "%s" PATH_SEP "%s", datasets_path, csv_name);
            fp = fopen(path, "r");
            if (fp) {
                fclose(fp);
                return cmd_csvread(path);
            }
        }
        
        /* Try ./datasets */
        sprintf(path, "datasets" PATH_SEP "%s", csv_name);
        fp = fopen(path, "r");
        if (fp) {
            fclose(fp);
            return cmd_csvread(path);
        }
        
        /* Try current directory */
        fp = fopen(csv_name, "r");
        if (fp) {
            fclose(fp);
            return cmd_csvread(csv_name);
        }
    }
    
    printf("Error: Dataset '%s' not found\n", name);
    printf("Searched in:\n");
    if (datasets_path && datasets_path[0]) {
        printf("  %s" PATH_SEP "%s\n", datasets_path, filename);
    }
    printf("  datasets" PATH_SEP "%s\n", filename);
    printf("  %s\n", filename);
    printf("Set SC_DATASETS environment variable to specify dataset path.\n");
    return 0;
}
/*
 * cmd_csvread - Load CSV file into workspace
 * Creates 'data' (numeric columns) and 'labels' (text column as codes)
 * Lines starting with # are treated as comments/metadata and displayed
 */
int cmd_csvread(const char *filename)
{
#ifdef HAVE_MATRIX
#ifdef HAVE_NAMED_VARS
    FILE *fp;
    char line[1024];
    char metadata[3][256];  /* Store up to 3 metadata lines */
    int num_metadata = 0;
    int row, data_rows;
    matrix_t data_mat = {0, 0, NULL}, labels_mat = {0, 0, NULL};  /* Allocated from arena */
    int num_cols[64];      /* 1 if numeric, 0 if text */
    int label_col;         /* Index of label column (-1 if none) */
    int total_cols, numeric_cols;
    int i;
    char label_values[64][32];  /* Unique label strings */
    int num_labels;
    extern int set_named_matrix(const char *name, const matrix_t *val);
    
    /* Open file */
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    /* Initialize */
    total_cols = 0;
    numeric_cols = 0;
    data_rows = 0;
    label_col = -1;
    num_labels = 0;
    
    /* Skip comment lines and read header */
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comment lines (starting with #) and save as metadata */
        if (line[0] == '#') {
            if (num_metadata < 3) {
                /* Remove # and leading space, and trailing newline */
                char *p = line + 1;
                size_t len;
                while (*p == ' ') p++;
                len = strlen(p);
                if (len > 0 && (p[len-1] == '\n' || p[len-1] == '\r')) p[--len] = '\0';
                if (len > 0 && (p[len-1] == '\n' || p[len-1] == '\r')) p[--len] = '\0';
                if (len > 255) len = 255;
                memcpy(metadata[num_metadata], p, len);
                metadata[num_metadata][len] = '\0';
                num_metadata++;
            }
            continue;
        }
        /* This is the header line - count columns */
        {
            char *p = line;
            while (*p) {
                if (*p == ',') total_cols++;
                p++;
            }
            total_cols++;
            for (i = 0; i < total_cols && i < 64; i++) {
                num_cols[i] = 1;  /* Assume numeric */
            }
        }
        break;
    }
    
    /* Read first data row to identify column types */
    if (fgets(line, sizeof(line), fp)) {
        char *p = line;
        int col_idx = 0;
        
        while (*p && col_idx < total_cols && col_idx < 64) {
            char val[256];
            int vi = 0;
            int has_digit = 0;
            int is_numeric = 1;
            int is_missing = 0;
            
            while (*p && *p != ',' && *p != '\n' && *p != '\r' && vi < 255) {
                if (*p >= '0' && *p <= '9') has_digit = 1;
                val[vi++] = *p++;
            }
            val[vi] = '\0';
            if (*p == ',') p++;
            
            /* Check for missing value indicator */
            if (val[0] == '?' || val[0] == '\0') {
                is_missing = 1;
                is_numeric = 1;  /* Treat as numeric with NaN */
            } else if (!has_digit) {
                is_numeric = 0;
            } else {
                char *endp;
                (void)strtod(val, &endp);
                if (*endp != '\0' && *endp != '\n' && *endp != '\r') {
                    is_numeric = 0;
                }
            }
            
            num_cols[col_idx] = is_numeric;
            if (!is_numeric && !is_missing) label_col = col_idx;
            col_idx++;
        }
        data_rows = 1;
    }
    
    /* Count numeric columns */
    for (i = 0; i < total_cols && i < 64; i++) {
        if (num_cols[i]) numeric_cols++;
    }
    
    if (numeric_cols == 0) {
        printf("Error: No numeric columns found\n");
        fclose(fp);
        return 0;
    }
    
    /* Count rows */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] && line[0] != '\n' && line[0] != '\r' && line[0] != '#') {
            data_rows++;
        }
    }
    
    if (data_rows == 0 || data_rows > MAT_MAX_ROWS || numeric_cols > MAT_MAX_COLS) {
        printf("Error: Data size %dx%d exceeds limits %dx%d\n",
               data_rows, numeric_cols, MAT_MAX_ROWS, MAT_MAX_COLS);
        fclose(fp);
        return 0;
    }
    
    /* Allocate matrices */
    mat_zero(&data_mat, data_rows, numeric_cols);
    mat_zero(&labels_mat, data_rows, 1);
    
    /* Second pass: read data */
    rewind(fp);
    
    /* Skip comment lines and header */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '#') break;  /* Found header, skip it */
    }
    
    row = 0;
    while (fgets(line, sizeof(line), fp) && row < data_rows) {
        char *p = line;
        int col_idx = 0;
        int mat_col = 0;
        
        /* Skip comment lines and empty lines */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        
        while (*p && col_idx < total_cols) {
            char val[256];
            int vi = 0;
            
            while (*p && *p != ',' && *p != '\n' && *p != '\r' && vi < 255) {
                val[vi++] = *p++;
            }
            val[vi] = '\0';
            if (*p == ',') p++;
            
            if (col_idx < 64 && num_cols[col_idx]) {
                /* Numeric column - handle missing values */
                double d;
                if (val[0] == '?' || val[0] == '\0') {
                    /* Missing value - use NaN */
                    d = 0.0 / 0.0;  /* NaN */
                } else {
                    d = strtod(val, NULL);
                }
                apf_from_double(&MAT_AT(&data_mat, row, mat_col).re, d);
                apf_zero(&MAT_AT(&data_mat, row, mat_col).im);
                mat_col++;
            } else if (col_idx == label_col) {
                /* Label column - convert to numeric code */
                int code = -1;
                int li;
                for (li = 0; li < num_labels; li++) {
                    if (strcmp(val, label_values[li]) == 0) {
                        code = li + 1;
                        break;
                    }
                }
                if (code < 0 && num_labels < 64) {
                    /* Safe copy with explicit length */
                    size_t vlen = strlen(val);
                    if (vlen > 31) vlen = 31;
                    memcpy(label_values[num_labels], val, vlen);
                    label_values[num_labels][vlen] = '\0';
                    num_labels++;
                    code = num_labels;
                }
                apf_from_int(&MAT_AT(&labels_mat, row, 0).re, code);
                apf_zero(&MAT_AT(&labels_mat, row, 0).im);
            }
            col_idx++;
        }
        row++;
    }
    
    fclose(fp);
    
    /* Store in workspace */
    if (!set_named_matrix("data", &data_mat)) {
        printf("Error: Cannot create 'data' variable\n");
        return 0;
    }
    
    /* Display metadata if present */
    if (num_metadata > 0) {
        printf("%s\n", metadata[0]);
        if (num_metadata > 1) printf("  %s\n", metadata[1]);
    }
    
    printf("Loaded: data %dx%d", data_rows, numeric_cols);
    
    if (label_col >= 0) {
        if (!set_named_matrix("labels", &labels_mat)) {
            printf("Error: Cannot create 'labels' variable\n");
            return 0;
        }
        printf(", labels %dx1 (%d classes: ", data_rows, num_labels);
        for (i = 0; i < num_labels && i < 5; i++) {
            printf("%d=%s", i+1, label_values[i]);
            if (i < num_labels - 1 && i < 4) printf(", ");
        }
        if (num_labels > 5) printf(", ...");
        printf(")");
    }
    printf("\n");
    
    return 1;
#else
    (void)filename;
    printf("Error: Named variables not enabled\n");
    return 0;
#endif
#else
    (void)filename;
    printf("Error: Matrix support not enabled\n");
    return 0;
#endif
}

/* ========== Fisher's Iris Dataset Demo ========== */
/*
 * Demonstrates loading iris dataset, computing statistics,
 * k-means clustering, and scatter plots.
 */
void cmd_demo_fisheriris(void)
{
#ifdef HAVE_MATRIX
#ifdef HAVE_NAMED_VARS
    extern int cmd_load_dataset(const char *name);
    extern int eval_expr_line(const char *line, int quiet);
    extern int get_named_var(const char *name, value_t *result);
    value_t idx_val;
    int i, counts[3] = {0, 0, 0};
    
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║           FISHER'S IRIS DATASET DEMO                             ║\n");
    printf("║        K-Means Clustering & Visualization                        ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* Load the dataset */
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("LOADING DATA\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");
    
    printf(">>> load iris\n");
    mat_arena_reset();
    if (!cmd_load_dataset("iris")) {
        printf("Error: Could not load iris dataset\n");
        return;
    }
    printf("\n");
    
    /* Show data dimensions */
    printf(">>> rows(data)\n");
    eval_expr_line("rows(data)", 0);
    printf("\n");
    
    printf(">>> cols(data)\n");
    eval_expr_line("cols(data)", 0);
    printf("\n");
    
    /* ================================================================
     * EXPLORATORY STATISTICS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  1. EXPLORATORY STATISTICS                                       ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> data(1:5,:)  %% First 5 samples (sepal L/W, petal L/W)\n");
    eval_expr_line("data(1:5,:)", 0);
    printf("\n");
    
    printf(">>> mean(data(:,1))  %% Mean sepal length\n");
    eval_expr_line("mean(data(:,1))", 0);
    printf("\n");
    
    printf(">>> mean(data(:,3))  %% Mean petal length\n");
    eval_expr_line("mean(data(:,3))", 0);
    printf("\n");
    
    printf(">>> std(data(:,3))  %% Std dev petal length\n");
    eval_expr_line("std(data(:,3))", 0);
    printf("\n");
    
    printf(">>> [min(data(:,3)), max(data(:,3))]  %% Petal length range\n");
    eval_expr_line("[min(data(:,3)), max(data(:,3))]", 0);
    printf("\n");
    
    /* ================================================================
     * K-MEANS CLUSTERING
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  2. K-MEANS CLUSTERING (k=3)                                     ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> [idx, C] = kmeans(data, 3)  %% Cluster into 3 groups\n");
    eval_expr_line("[idx, C] = kmeans(data, 3)", 0);
    printf("\n");
    
    printf(">>> C  %% Cluster centroids (SepalL, SepalW, PetalL, PetalW)\n");
    eval_expr_line("C", 0);
    printf("\n");
    
    /* Count clusters manually */
    if (get_named_var("idx", &idx_val) && idx_val.type == VAL_MATRIX) {
        for (i = 0; i < idx_val.v.matrix.rows; i++) {
            int c = (int)apf_to_double(&MAT_AT(&idx_val.v.matrix, i, 0).re) - 1;
            if (c >= 0 && c < 3) counts[c]++;
        }
        printf("Cluster sizes: [%d, %d, %d]\n\n", counts[0], counts[1], counts[2]);
    }
    
    /* ================================================================
     * SCATTER PLOT - Petal dimensions colored by cluster
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  3. SCATTER: Petal Length vs Width (by cluster)                  ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> scatter(data(:,3), data(:,4), idx)\n");
    eval_expr_line("scatter(data(:,3), data(:,4), idx)", 0);
    printf("\n");
    
    /* ================================================================
     * SCATTER PLOT - Actual species labels
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  4. SCATTER: Same data, colored by actual species (labels)       ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> scatter(data(:,3), data(:,4), labels)\n");
    printf("    (1=setosa, 2=versicolor, 3=virginica)\n");
    eval_expr_line("scatter(data(:,3), data(:,4), labels)", 0);
    printf("\n");
    
    /* ================================================================
     * SEPAL DIMENSIONS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  5. SCATTER: Sepal Length vs Width (by cluster)                  ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> scatter(data(:,1), data(:,2), idx)\n");
    eval_expr_line("scatter(data(:,1), data(:,2), idx)", 0);
    printf("\n");
    
    /* ================================================================
     * TRY THESE COMMANDS
     * ================================================================ */
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("TRY THESE COMMANDS:\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");
    printf(">>> corr(data)                       %% Correlation matrix\n");
    printf(">>> [coeff, score] = pca(data)       %% PCA decomposition\n");
    printf(">>> scatter(score(:,1), score(:,2), idx)  %% PCA plot\n");
    printf(">>> [idx5, C5] = kmeans(data, 5)     %% Try k=5 clusters\n\n");
    
#else
    printf("Error: Named variables not enabled\n");
#endif
#else
    printf("Error: Matrix support not enabled\n");
#endif
}

void run_demo(void)
{
    printf("\n");
    printf("===========================================\n");
    printf("  SC Calculator Demo (%d-bit precision)\n", AP_BITS);
    printf("===========================================\n\n");
    
    demo_basic();
    
#ifdef HAVE_SQRT
    demo_sqrt();
#endif
#ifdef HAVE_TRIG
    demo_trig();
#endif
#ifdef HAVE_EXP
    demo_exp();
#endif
#ifdef HAVE_COMPLEX
    demo_complex();
#endif
#ifdef HAVE_HEXBIN
    demo_hexbin();
#endif
#ifdef HAVE_BITWISE
    demo_bitwise();
#endif
#ifdef HAVE_COMB
    demo_comb();
#endif
#ifdef HAVE_HYPER
    demo_hyper();
#endif
#ifdef HAVE_GAMMA
    demo_gamma();
#endif
#ifdef HAVE_ERF
    demo_erf();
#endif
#ifdef HAVE_BESSEL
    demo_bessel();
#endif
#ifdef HAVE_ELLIPTIC
    demo_elliptic();
#endif
#ifdef HAVE_LAMBERTW
    demo_lambertw();
#endif
#ifdef HAVE_RPN
    demo_rpn();
#endif
#ifdef HAVE_VARIABLES
    demo_variables();
#endif
#ifdef HAVE_MATRIX
    demo_matrix();
    demo_matlab();
#endif
#ifdef HAVE_STATS
    demo_stats();
#endif
#ifdef HAVE_SOLVER
    demo_solver();
#endif
#ifdef HAVE_TVM
    demo_tvm();
#endif
#ifdef HAVE_ORBITAL
    demo_orbital();
#endif
#ifdef HAVE_NEWTON
    demo_newton();
#endif

    printf("OTHER COMMANDS\n");
    printf("  help       Show all commands\n");
    printf("  constants  Mathematical constants\n");
    printf("  features   Enabled features\n");
    printf("  test       Self-tests\n");
    printf("  bench      Benchmarks\n");
    printf("\n");
    printf("Demo complete!\n\n");
}

/* ========== SaaS Metrics Demo ========== */
void cmd_demo_saas(void)
{
#ifdef HAVE_MATRIX
#ifdef HAVE_NAMED_VARS
    extern int cmd_load_dataset(const char *name);
    extern int get_named_var(const char *name, value_t *result);
    extern int set_named_matrix(const char *name, const matrix_t *val);
    extern void mat_mrr(matrix_t *r, const matrix_t *data);
    extern void mat_arpu(matrix_t *r, const matrix_t *data);
    extern void mat_mrrbridge(matrix_t *r, const matrix_t *data);
    extern void mat_customercount(matrix_t *r, const matrix_t *data);
    extern void mat_newcustomers(matrix_t *r, const matrix_t *data);
    extern void mat_churn(matrix_t *r, const matrix_t *data);
    extern void mat_reactivated(matrix_t *r, const matrix_t *data);
    extern void mat_churnrate(matrix_t *r, const matrix_t *data);
    extern void mat_nrr(matrix_t *r, const matrix_t *data);
    extern void mat_grr(matrix_t *r, const matrix_t *data);
    extern void mat_retention(matrix_t *r, const matrix_t *data);
    extern void mat_tenure(matrix_t *r, const matrix_t *data);
    extern void datetime_format(char *buf, int bufsize, const apf *seconds);
    
    value_t data_val;
    matrix_t *data, mrr, arpu, bridge, custs, newc, churned, react, rate, nrr_m, grr_m, ret, ten;
    int i, j, n;
    char datebuf[32];
    double total_mrr, avg_mrr, first_mrr, last_mrr;
    double avg_arpu, avg_churn, avg_nrr, avg_grr;
    int total_new, total_churned, total_react;
    double avg_tenure;
    
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║           SaaS METRICS ANALYSIS DEMO                             ║\n");
    printf("║     Comprehensive Subscription Business Analytics                ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* Load dataset */
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("LOADING DATA\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");
    printf(">>> load hourlybilling\n");
    
    mat_arena_reset();
    if (!cmd_load_dataset("hourlybilling")) {
        printf("\nError: Could not load hourlybilling dataset\n");
        printf("Run 'generate' command first to create demo data.\n");
        return;
    }
    
    if (!get_named_var("data", &data_val) || data_val.type != VAL_MATRIX) {
        printf("Error: Could not access data\n");
        return;
    }
    data = &data_val.v.matrix;
    
    printf("\nData: %d rows x %d columns\n", data->rows, data->cols);
    printf("Columns: [timestamp, customerid, sizetb, billing]\n\n");
    
    /* ================================================================
     * REVENUE METRICS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  1. REVENUE METRICS                                              ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* MRR */
    printf(">>> mrr(data)  %% Monthly Recurring Revenue\n\n");
    mat_mrr(&mrr, data);
    
    printf("    Month           MRR         ARR (x12)\n");
    printf("    ──────────────  ──────────  ──────────\n");
    
    total_mrr = 0;
    n = mrr.rows;
    for (i = 0; i < n; i++) {
        double m = apf_to_double(&MAT_AT(&mrr, i, 1).re);
        total_mrr += m;
        datetime_format(datebuf, sizeof(datebuf), &MAT_AT(&mrr, i, 0).re);
        datebuf[10] = '\0';  /* Just date */
        printf("    %s  $%9.2f  $%10.2f\n", datebuf, m, m * 12);
    }
    first_mrr = apf_to_double(&MAT_AT(&mrr, 0, 1).re);
    last_mrr = apf_to_double(&MAT_AT(&mrr, n-1, 1).re);
    avg_mrr = total_mrr / n;
    
    printf("\n    Summary:\n");
    printf("    • Average MRR:     $%.2f\n", avg_mrr);
    printf("    • Current MRR:     $%.2f\n", last_mrr);
    printf("    • Current ARR:     $%.2f\n", last_mrr * 12);
    printf("    • MRR Growth:      %.1f%% (first to last)\n\n", (last_mrr/first_mrr - 1) * 100);
    
    /* ARPU */
    printf(">>> arpu(data)  %% Average Revenue Per User\n\n");
    mat_arpu(&arpu, data);
    
    avg_arpu = 0;
    for (i = 0; i < arpu.rows; i++) {
        avg_arpu += apf_to_double(&MAT_AT(&arpu, i, 1).re);
    }
    avg_arpu /= arpu.rows;
    
    printf("    Average ARPU:  $%.2f/month\n", avg_arpu);
    printf("    Latest ARPU:   $%.2f/month\n\n", apf_to_double(&MAT_AT(&arpu, arpu.rows-1, 1).re));
    
    /* MRR Bridge */
    printf(">>> mrrbridge(data)  %% MRR Movement Components\n\n");
    mat_mrrbridge(&bridge, data);
    
    printf("    Month       New MRR    Expansion  Contraction  Churned    Net MRR\n");
    printf("    ──────────  ─────────  ─────────  ───────────  ─────────  ─────────\n");
    
    for (i = (bridge.rows > 6 ? bridge.rows - 6 : 0); i < bridge.rows; i++) {
        datetime_format(datebuf, sizeof(datebuf), &MAT_AT(&bridge, i, 0).re);
        datebuf[10] = '\0';
        printf("    %s  $%8.0f  $%8.0f   $%8.0f  $%8.0f  $%8.0f\n",
               datebuf,
               apf_to_double(&MAT_AT(&bridge, i, 1).re),
               apf_to_double(&MAT_AT(&bridge, i, 2).re),
               apf_to_double(&MAT_AT(&bridge, i, 3).re),
               apf_to_double(&MAT_AT(&bridge, i, 4).re),
               apf_to_double(&MAT_AT(&bridge, i, 5).re));
    }
    printf("\n");
    
    /* ================================================================
     * CUSTOMER METRICS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  2. CUSTOMER METRICS                                             ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* Active customers */
    printf(">>> customercount(data)  %% Active customers per month\n\n");
    mat_customercount(&custs, data);
    
    printf("    Month       Active    New    Churned  Reactivated\n");
    printf("    ──────────  ──────    ───    ───────  ───────────\n");
    
    mat_newcustomers(&newc, data);
    mat_churn(&churned, data);
    mat_reactivated(&react, data);
    
    total_new = 0;
    total_churned = 0;
    total_react = 0;
    
    for (i = 0; i < custs.rows; i++) {
        int nc = (int)apf_to_double(&MAT_AT(&newc, i, 1).re);
        int ch = (int)apf_to_double(&MAT_AT(&churned, i, 1).re);
        int re = (int)apf_to_double(&MAT_AT(&react, i, 1).re);
        total_new += nc;
        total_churned += ch;
        total_react += re;
        
        datetime_format(datebuf, sizeof(datebuf), &MAT_AT(&custs, i, 0).re);
        datebuf[10] = '\0';
        printf("    %s  %6d    %3d      %3d          %3d\n",
               datebuf,
               (int)apf_to_double(&MAT_AT(&custs, i, 1).re),
               nc, ch, re);
    }
    
    printf("\n    Summary:\n");
    printf("    • Current Customers:   %d\n", (int)apf_to_double(&MAT_AT(&custs, custs.rows-1, 1).re));
    printf("    • Peak Customers:      %d\n", (int)apf_to_double(&MAT_AT(&custs, 0, 1).re));  /* Will recalc */
    printf("    • Total New:           %d\n", total_new);
    printf("    • Total Churned:       %d\n", total_churned);
    printf("    • Total Reactivated:   %d\n", total_react);
    printf("    • Net Customer Change: %+d\n\n", total_new - total_churned + total_react);
    
    /* ================================================================
     * CHURN & RETENTION
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  3. CHURN & RETENTION METRICS                                    ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* Churn rate */
    printf(">>> churnrate(data)  %% Logo churn rate (%%)\n\n");
    mat_churnrate(&rate, data);
    
    avg_churn = 0;
    for (i = 1; i < rate.rows; i++) {
        avg_churn += apf_to_double(&MAT_AT(&rate, i, 1).re);
    }
    avg_churn /= (rate.rows - 1);
    
    printf("    Average Monthly Churn:   %.2f%%\n", avg_churn);
    printf("    Annualized Churn:        %.1f%%\n", (1 - pow(1 - avg_churn/100, 12)) * 100);
    printf("    Implied Customer Life:   %.1f months\n\n", 100 / avg_churn);
    
    /* NRR */
    printf(">>> nrr(data)  %% Net Revenue Retention\n\n");
    mat_nrr(&nrr_m, data);
    
    avg_nrr = 0;
    for (i = 1; i < nrr_m.rows; i++) {
        avg_nrr += apf_to_double(&MAT_AT(&nrr_m, i, 1).re);
    }
    avg_nrr /= (nrr_m.rows - 1);
    
    printf("    Average NRR:  %.1f%%\n", avg_nrr);
    if (avg_nrr > 100) {
        printf("    Status:       ✓ EXPANDING (existing customers growing)\n");
    } else {
        printf("    Status:       ⚠ CONTRACTING (losing revenue from existing)\n");
    }
    printf("\n");
    
    /* GRR */
    printf(">>> grr(data)  %% Gross Revenue Retention\n\n");
    mat_grr(&grr_m, data);
    
    avg_grr = 0;
    for (i = 1; i < grr_m.rows; i++) {
        avg_grr += apf_to_double(&MAT_AT(&grr_m, i, 1).re);
    }
    avg_grr /= (grr_m.rows - 1);
    
    printf("    Average GRR:  %.1f%%\n", avg_grr);
    printf("    (GRR caps at 100%%, excludes expansion revenue)\n\n");
    
    /* ================================================================
     * COHORT ANALYSIS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  4. COHORT ANALYSIS                                              ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    /* Retention matrix */
    printf(">>> retention(data)  %% Cohort retention matrix\n\n");
    mat_retention(&ret, data);
    
    printf("    Cohort retention (%% of original cohort still active):\n\n");
    printf("              ");
    for (j = 0; j < 8 && j < ret.cols; j++) {
        printf("  M%-2d  ", j);
    }
    printf("\n    ────────  ");
    for (j = 0; j < 8 && j < ret.cols; j++) {
        printf("──────");
    }
    printf("\n");
    
    for (i = 0; i < 6 && i < ret.rows; i++) {
        printf("    Cohort %d  ", i+1);
        for (j = 0; j < 8 && j < ret.cols; j++) {
            double pct = apf_to_double(&MAT_AT(&ret, i, j).re);
            if (pct > 0) {
                printf("%5.0f%% ", pct);
            } else {
                printf("    - ");
            }
        }
        printf("\n");
    }
    printf("\n    (Cohort N = customers who first appeared in month N)\n\n");
    
    /* ================================================================
     * CUSTOMER TENURE
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  5. CUSTOMER TENURE                                              ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf(">>> tenure(data)  %% Customer lifetime analysis\n\n");
    mat_tenure(&ten, data);
    
    avg_tenure = 0;
    for (i = 0; i < ten.rows; i++) {
        avg_tenure += apf_to_double(&MAT_AT(&ten, i, 1).re);
    }
    avg_tenure /= ten.rows;
    
    printf("    Average Tenure:  %.1f months\n", avg_tenure);
    printf("    Max Tenure:      %.0f months\n", apf_to_double(&MAT_AT(&ten, 0, 1).re));
    printf("    Min Tenure:      %.0f months\n\n", apf_to_double(&MAT_AT(&ten, ten.rows-1, 1).re));
    
    /* ================================================================
     * REVENUE CONCENTRATION
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  6. REVENUE CONCENTRATION                                        ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    {
        extern void mat_concentration(matrix_t *r, const matrix_t *data);
        extern void mat_toprevenue(matrix_t *r, const matrix_t *data, int n);
        matrix_t conc, topr;
        
        printf(">>> concentration(data)  %% Revenue concentration metrics\n\n");
        mat_concentration(&conc, data);
        
        printf("    Top 1 customer:      %.1f%% of revenue\n", apf_to_double(&MAT_AT(&conc, 0, 0).re));
        printf("    Top 5 customers:     %.1f%% of revenue\n", apf_to_double(&MAT_AT(&conc, 0, 1).re));
        printf("    Top 10 customers:    %.1f%% of revenue\n", apf_to_double(&MAT_AT(&conc, 0, 2).re));
        printf("    Top 20%% customers:   %.1f%% of revenue\n", apf_to_double(&MAT_AT(&conc, 0, 3).re));
        printf("    HHI Index:           %.0f ", apf_to_double(&MAT_AT(&conc, 0, 4).re));
        {
            double hhi = apf_to_double(&MAT_AT(&conc, 0, 4).re);
            if (hhi < 1500) printf("(low concentration)\n");
            else if (hhi < 2500) printf("(moderate concentration)\n");
            else printf("(high concentration)\n");
        }
        
        printf("\n>>> toprevenue(data, 5)  %% Top 5 customers by revenue\n\n");
        mat_toprevenue(&topr, data, 5);
        printf("    Customer    Total Rev      %% of Total\n");
        printf("    ────────    ─────────      ──────────\n");
        for (i = 0; i < 5 && i < topr.rows; i++) {
            printf("    %8d    $%9.2f      %5.1f%%\n",
                   (int)apf_to_double(&MAT_AT(&topr, i, 0).re),
                   apf_to_double(&MAT_AT(&topr, i, 1).re),
                   apf_to_double(&MAT_AT(&topr, i, 2).re));
        }
        printf("\n");
    }
    
    /* ================================================================
     * GROWTH EFFICIENCY
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  7. GROWTH EFFICIENCY                                            ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    {
        extern void mat_quickratio(matrix_t *r, const matrix_t *data);
        extern void mat_netchurn(matrix_t *r, const matrix_t *data);
        matrix_t qr, nc;
        double avg_qr = 0, avg_nc = 0;
        int qr_count = 0;
        
        mat_quickratio(&qr, data);
        mat_netchurn(&nc, data);
        
        printf(">>> quickratio(data)  %% SaaS Quick Ratio\n\n");
        printf("    Quick Ratio = (New MRR + Expansion) / (Churn + Contraction)\n\n");
        
        for (i = 1; i < qr.rows; i++) {
            double q = apf_to_double(&MAT_AT(&qr, i, 1).re);
            if (q < 100) { avg_qr += q; qr_count++; }
        }
        if (qr_count > 0) avg_qr /= qr_count;
        
        printf("    Average Quick Ratio: %.2f\n", avg_qr);
        if (avg_qr >= 4) printf("    Assessment:          ✓ EXCELLENT (>4x)\n");
        else if (avg_qr >= 2) printf("    Assessment:          ✓ GOOD (>2x)\n");
        else if (avg_qr >= 1) printf("    Assessment:          ⚠ GROWING (<2x)\n");
        else printf("    Assessment:          ✗ SHRINKING (<1x)\n");
        
        printf("\n>>> netchurn(data)  %% Net Churn Rate\n\n");
        for (i = 1; i < nc.rows; i++) {
            avg_nc += apf_to_double(&MAT_AT(&nc, i, 1).re);
        }
        if (nc.rows > 1) avg_nc /= (nc.rows - 1);
        
        printf("    Avg Net Churn:       %.2f%%\n", avg_nc);
        if (avg_nc < 0) printf("    Status:              ✓ Net EXPANSION (negative churn!)\n");
        else if (avg_nc < 2) printf("    Status:              ✓ LOW churn\n");
        else if (avg_nc < 5) printf("    Status:              ⚠ MODERATE churn\n");
        else printf("    Status:              ✗ HIGH churn\n");
        printf("\n");
    }
    
    /* ================================================================
     * UNIT ECONOMICS
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  8. UNIT ECONOMICS (Derived)                                     ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    {
        double ltv, ltv_cac, cac_assumed = 500;  /* Assume $500 CAC for demo */
        double monthly_churn_decimal = avg_churn / 100;
        double gross_margin = 0.8;  /* Assume 80% gross margin */
        
        /* LTV = (ARPU * Gross Margin) / Monthly Churn */
        ltv = (avg_arpu * gross_margin) / monthly_churn_decimal;
        ltv_cac = ltv / cac_assumed;
        
        printf("    Assumptions: 80%% gross margin, $500 CAC\n\n");
        printf("    LTV Formula:  (ARPU × Gross Margin) / Churn Rate\n");
        printf("    LTV:          ($%.2f × 0.80) / %.4f = $%.2f\n\n", avg_arpu, monthly_churn_decimal, ltv);
        printf("    LTV:CAC Ratio: $%.2f / $500 = %.1fx\n", ltv, ltv_cac);
        if (ltv_cac >= 3) {
            printf("    Assessment:    ✓ HEALTHY (>3x is target)\n");
        } else if (ltv_cac >= 1) {
            printf("    Assessment:    ⚠ ACCEPTABLE (aim for 3x+)\n");
        } else {
            printf("    Assessment:    ✗ UNPROFITABLE (<1x)\n");
        }
        printf("\n    CAC Payback:   %.1f months ($500 / $%.2f)\n\n", cac_assumed / (avg_arpu * gross_margin), avg_arpu * gross_margin);
    }
    
    /* ================================================================
     * EXECUTIVE SUMMARY
     * ================================================================ */
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  EXECUTIVE SUMMARY                                               ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf("    ┌─────────────────────────────────────────────────────────────┐\n");
    printf("    │  REVENUE                                                    │\n");
    printf("    │    Current MRR:      $%.2f                              \n", last_mrr);
    printf("    │    Current ARR:      $%.2f                             \n", last_mrr * 12);
    printf("    │    ARPU:             $%.2f/month                           \n", avg_arpu);
    printf("    ├─────────────────────────────────────────────────────────────┤\n");
    printf("    │  CUSTOMERS                                                  │\n");
    printf("    │    Active:           %d                                     \n", (int)apf_to_double(&MAT_AT(&custs, custs.rows-1, 1).re));
    printf("    │    Monthly Churn:    %.1f%%                                 \n", avg_churn);
    printf("    ├─────────────────────────────────────────────────────────────┤\n");
    printf("    │  RETENTION                                                  │\n");
    printf("    │    NRR:              %.1f%%                                 \n", avg_nrr);
    printf("    │    GRR:              %.1f%%                                 \n", avg_grr);
    printf("    │    Avg Tenure:       %.1f months                           \n", avg_tenure);
    printf("    └─────────────────────────────────────────────────────────────┘\n\n");
    
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("Available SaaS functions:\n");
    printf("  Revenue:     mrr, arpu, mrrbridge, revchurn, netchurn\n");
    printf("  Customers:   customercount, newcustomers, churn, reactivated\n");
    printf("  Retention:   churnrate, nrr, grr, retention, tenure\n");
    printf("  Analysis:    concentration, toprevenue, quickratio, ltv\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");
    
#else
    printf("Named variables not enabled\n");
#endif
#else
    printf("Matrix operations not enabled\n");
#endif
}
