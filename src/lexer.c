/* lexer.c - Tokenizer for scalc
 * Handles decimal, hex (0x), binary (0b), and imaginary numbers
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"

/* Global tokenizer state */
const char *input_ptr = NULL;
token_t current_token;

void skip_whitespace(void)
{
    while (*input_ptr == ' ' || *input_ptr == '\t') input_ptr++;
}

static int is_alpha(char c)
{
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

static int isxdigit_local(char c)
{
    return (c >= '0' && c <= '9') || 
           (c >= 'a' && c <= 'f') || 
           (c >= 'A' && c <= 'F');
}

void next_token(void)
{
    char func_buf[16];
    int func_len;

    skip_whitespace();

    if (*input_ptr == '\0' || *input_ptr == '\n') {
        current_token.type = TOK_END;
        return;
    }

    /* Hex number (0x...) */
    if (*input_ptr == '0' && (*(input_ptr+1) == 'x' || *(input_ptr+1) == 'X')) {
        apf int_part, sixteen, digit_apf, temp;
        input_ptr += 2;  /* skip 0x */
        apf_from_int(&int_part, 0);
        apf_from_int(&sixteen, 16);
        
        while (isxdigit_local(*input_ptr)) {
            int d;
            char c = *input_ptr;
            if (c >= '0' && c <= '9') d = c - '0';
            else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
            else d = c - 'A' + 10;
            
            apf_mul(&temp, &int_part, &sixteen);
            apf_from_int(&digit_apf, d);
            apf_add(&int_part, &temp, &digit_apf);
            input_ptr++;
        }
        
        current_token.type = TOK_NUM;
        apf_copy(&current_token.value, &int_part);
        return;
    }
    
    /* Binary number (0b...) */
    if (*input_ptr == '0' && (*(input_ptr+1) == 'b' || *(input_ptr+1) == 'B')) {
        apf int_part, two, digit_apf, temp;
        input_ptr += 2;  /* skip 0b */
        apf_from_int(&int_part, 0);
        apf_from_int(&two, 2);
        
        while (*input_ptr == '0' || *input_ptr == '1') {
            int d = *input_ptr - '0';
            apf_mul(&temp, &int_part, &two);
            apf_from_int(&digit_apf, d);
            apf_add(&int_part, &temp, &digit_apf);
            input_ptr++;
        }
        
        current_token.type = TOK_NUM;
        apf_copy(&current_token.value, &int_part);
        return;
    }

    /* Decimal number */
    if (isdigit(*input_ptr) || (*input_ptr == '.' && isdigit(*(input_ptr + 1)))) {
        apf int_part, ten, digit_apf, temp;

        /* Parse integer part using arbitrary precision */
        apf_from_int(&int_part, 0);
        apf_from_int(&ten, 10);
        
        while (isdigit(*input_ptr)) {
            int d = *input_ptr - '0';
            apf_mul(&temp, &int_part, &ten);
            apf_from_int(&digit_apf, d);
            apf_add(&int_part, &temp, &digit_apf);
            input_ptr++;
        }

        /* Parse fractional part using arbitrary precision */
        if (*input_ptr == '.') {
            apf frac, divisor;
            apf_from_int(&frac, 0);
            apf_from_int(&divisor, 1);
            input_ptr++;
            
            while (isdigit(*input_ptr)) {
                int d = *input_ptr - '0';
                /* frac = frac * 10 + d */
                apf_mul(&temp, &frac, &ten);
                apf_from_int(&digit_apf, d);
                apf_add(&frac, &temp, &digit_apf);
                /* divisor *= 10 */
                apf_mul(&temp, &divisor, &ten);
                apf_copy(&divisor, &temp);
                input_ptr++;
            }
            
            /* Add frac/divisor to int_part */
            if (!apf_is_zero(&frac)) {
                apf_div(&temp, &frac, &divisor);
                apf_add(&frac, &int_part, &temp);
                apf_copy(&int_part, &frac);
            }
        }

        current_token.type = TOK_NUM;
        apf_copy(&current_token.value, &int_part);
        
        /* Parse scientific notation (e.g., 1e10, 1.5e-3, 2E6) */
        if (*input_ptr == 'e' || *input_ptr == 'E') {
            const char *exp_start = input_ptr;
            int exp_sign = 0;
            long exp_val = 0;
            apf multiplier;
            
            input_ptr++;
            
            /* Optional sign */
            if (*input_ptr == '-') {
                exp_sign = 1;
                input_ptr++;
            } else if (*input_ptr == '+') {
                input_ptr++;
            }
            
            /* Exponent digits */
            if (isdigit(*input_ptr)) {
                while (isdigit(*input_ptr)) {
                    exp_val = exp_val * 10 + (*input_ptr - '0');
                    if (exp_val > 1000000000L) exp_val = 1000000000L; /* Cap exponent */
                    input_ptr++;
                }
                
                /* Apply 10^exp to the number */
                apf_from_int(&multiplier, 10);
                if (exp_sign) {
                    /* Divide by 10^exp_val */
                    while (exp_val > 0) {
                        apf_div(&temp, &int_part, &multiplier);
                        apf_copy(&int_part, &temp);
                        exp_val--;
                    }
                } else {
                    /* Multiply by 10^exp_val */
                    while (exp_val > 0) {
                        apf_mul(&temp, &int_part, &multiplier);
                        apf_copy(&int_part, &temp);
                        exp_val--;
                    }
                }
                apf_copy(&current_token.value, &int_part);
            } else {
                /* Not a valid exponent, backtrack */
                input_ptr = exp_start;
            }
        }

        /* Check for trailing 'i' for imaginary */
        if (*input_ptr == 'i' && !is_alpha(*(input_ptr + 1))) {
            input_ptr++;
            current_token.type = TOK_IMAG;
        }
        return;
    }

    /* Standalone 'i' */
    if (*input_ptr == 'i' && !is_alpha(*(input_ptr + 1))) {
        input_ptr++;
        current_token.type = TOK_IMAG;
        apf_from_int(&current_token.value, 1);
        return;
    }

    /* Function or constant name */
    if (is_alpha(*input_ptr)) {
        func_len = 0;
        while (is_alpha(*input_ptr) && func_len < 15) {
            func_buf[func_len++] = *input_ptr++;
        }
        while (isdigit(*input_ptr) && func_len < 15) {
            func_buf[func_len++] = *input_ptr++;
        }
        func_buf[func_len] = '\0';
        
        /* Check for keywords */
        if (strcmp(func_buf, "for") == 0) {
            current_token.type = TOK_FOR;
            return;
        }
        if (strcmp(func_buf, "end") == 0) {
            current_token.type = TOK_ENDFOR;
            return;
        }
        if (strcmp(func_buf, "in") == 0) {
            current_token.type = TOK_IN;
            return;
        }
        
        /* Check for ans (previous answer) - case insensitive */
        if ((func_buf[0] == 'a' || func_buf[0] == 'A') &&
            (func_buf[1] == 'n' || func_buf[1] == 'N') &&
            (func_buf[2] == 's' || func_buf[2] == 'S') &&
            func_buf[3] == '\0') {
            current_token.type = TOK_ANS;
            return;
        }
        
        /* Check for Roman numeral (all uppercase M, D, C, L, X, V, I) */
        {
            int is_roman = 1;
            int ri;
            long roman_val = 0;
            
            for (ri = 0; func_buf[ri]; ri++) {
                char rc = func_buf[ri];
                if (rc != 'M' && rc != 'D' && rc != 'C' && rc != 'L' && 
                    rc != 'X' && rc != 'V' && rc != 'I') {
                    is_roman = 0;
                    break;
                }
            }
            
            if (is_roman && func_len > 0) {
                /* Parse Roman numeral with subtractive notation */
                for (ri = 0; func_buf[ri]; ri++) {
                    long cur_val = 0;
                    long next_val = 0;
                    char rc = func_buf[ri];
                    char rn = func_buf[ri + 1];
                    
                    switch (rc) {
                        case 'M': cur_val = 1000; break;
                        case 'D': cur_val = 500; break;
                        case 'C': cur_val = 100; break;
                        case 'L': cur_val = 50; break;
                        case 'X': cur_val = 10; break;
                        case 'V': cur_val = 5; break;
                        case 'I': cur_val = 1; break;
                    }
                    
                    if (rn) {
                        switch (rn) {
                            case 'M': next_val = 1000; break;
                            case 'D': next_val = 500; break;
                            case 'C': next_val = 100; break;
                            case 'L': next_val = 50; break;
                            case 'X': next_val = 10; break;
                            case 'V': next_val = 5; break;
                            case 'I': next_val = 1; break;
                        }
                    }
                    
                    if (next_val > cur_val) {
                        roman_val -= cur_val;
                    } else {
                        roman_val += cur_val;
                    }
                }
                
                current_token.type = TOK_NUM;
                apf_from_int(&current_token.value, roman_val);
                return;
            }
        }
        
        current_token.type = TOK_FUNC;
        strcpy(current_token.func_name, func_buf);
        return;
    }

    /* Operators */
    switch (*input_ptr) {
        case '+': current_token.type = TOK_PLUS; input_ptr++; return;
        case '-': current_token.type = TOK_MINUS; input_ptr++; return;
        case '*': current_token.type = TOK_MUL; input_ptr++; return;
        case '/': current_token.type = TOK_DIV; input_ptr++; return;
        case '\\': current_token.type = TOK_BACKSLASH; input_ptr++; return;
        case '^': current_token.type = TOK_POW; input_ptr++; return;
        case '!': current_token.type = TOK_FACT; input_ptr++; return;
        case '(': current_token.type = TOK_LPAREN; input_ptr++; return;
        case ')': current_token.type = TOK_RPAREN; input_ptr++; return;
        case '=': current_token.type = TOK_ASSIGN; input_ptr++; return;
        case ';': current_token.type = TOK_SEMI; input_ptr++; return;
        case '[': current_token.type = TOK_LBRACKET; input_ptr++; return;
        case ']': current_token.type = TOK_RBRACKET; input_ptr++; return;
        case ':': current_token.type = TOK_COLON; input_ptr++; return;
        case ',': current_token.type = TOK_COMMA; input_ptr++; return;
        case '\'': current_token.type = TOK_TRANSPOSE; input_ptr++; return;
        case '.':
            /* Check for element-wise operators: .* ./ .^ .' */
            if (*(input_ptr + 1) == '*') {
                current_token.type = TOK_DOT_MUL;
                input_ptr += 2;
                return;
            }
            if (*(input_ptr + 1) == '/') {
                current_token.type = TOK_DOT_DIV;
                input_ptr += 2;
                return;
            }
            if (*(input_ptr + 1) == '^') {
                current_token.type = TOK_DOT_POW;
                input_ptr += 2;
                return;
            }
            if (*(input_ptr + 1) == '\'') {
                current_token.type = TOK_DOT_TRANSPOSE;
                input_ptr += 2;
                return;
            }
            /* Otherwise fall through to number parsing (for .5 etc) */
            if (isdigit(*(input_ptr + 1))) {
                /* Will be handled by number parsing above */
            }
            current_token.type = TOK_ERROR;
            return;
        default: current_token.type = TOK_ERROR; return;
    }
}
