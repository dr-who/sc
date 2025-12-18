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
    char func_buf[32];
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
            apf pow10, base10, temp2;
            
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
                
                /* Compute 10^exp_val using fast exponentiation by squaring: O(log N) */
                apf_from_int(&pow10, 1);
                apf_from_int(&base10, 10);
                while (exp_val > 0) {
                    if (exp_val & 1) {
                        apf_mul(&temp2, &pow10, &base10);
                        apf_copy(&pow10, &temp2);
                    }
                    apf_mul(&temp2, &base10, &base10);
                    apf_copy(&base10, &temp2);
                    exp_val >>= 1;
                }
                
                /* Apply 10^exp to the number */
                if (exp_sign) {
                    apf_div(&temp, &int_part, &pow10);
                } else {
                    apf_mul(&temp, &int_part, &pow10);
                }
                apf_copy(&current_token.value, &temp);
            } else {
                /* Not a valid exponent, backtrack */
                input_ptr = exp_start;
            }
        }

        /* Check for trailing 'i' or 'j' for imaginary */
        if ((*input_ptr == 'i' || *input_ptr == 'j') && !is_alpha(*(input_ptr + 1))) {
            input_ptr++;
            current_token.type = TOK_IMAG;
        }
        return;
    }

    /* Standalone 'i' or 'j' */
    if ((*input_ptr == 'i' || *input_ptr == 'j') && !is_alpha(*(input_ptr + 1))) {
        input_ptr++;
        current_token.type = TOK_IMAG;
        apf_from_int(&current_token.value, 1);
        return;
    }

    /* Function or constant name */
    if (is_alpha(*input_ptr) || *input_ptr == '_') {
        func_len = 0;
        while ((is_alpha(*input_ptr) || isdigit(*input_ptr) || *input_ptr == '_') && func_len < 31) {
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
        /* Boolean operators */
        if (strcmp(func_buf, "and") == 0) {
            current_token.type = TOK_AND;
            return;
        }
        if (strcmp(func_buf, "or") == 0) {
            current_token.type = TOK_OR;
            return;
        }
        if (strcmp(func_buf, "not") == 0) {
            current_token.type = TOK_NOT;
            return;
        }
        if (strcmp(func_buf, "xor") == 0) {
            current_token.type = TOK_XOR;
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
        
        current_token.type = TOK_FUNC;
        strcpy(current_token.func_name, func_buf);
        return;
    }

    /* Operators */
    /* Check for Unicode multiplication and division symbols */
    /* × U+00D7 (UTF-8: 0xC3 0x97) */
    if ((unsigned char)*input_ptr == 0xC3 && (unsigned char)*(input_ptr+1) == 0x97) {
        current_token.type = TOK_MUL;
        input_ptr += 2;
        return;
    }
    /* ÷ U+00F7 (UTF-8: 0xC3 0xB7) */
    if ((unsigned char)*input_ptr == 0xC3 && (unsigned char)*(input_ptr+1) == 0xB7) {
        current_token.type = TOK_DIV;
        input_ptr += 2;
        return;
    }
    /* ✕ U+2715 (UTF-8: 0xE2 0x9C 0x95) - ballot X */
    if ((unsigned char)*input_ptr == 0xE2 && (unsigned char)*(input_ptr+1) == 0x9C && (unsigned char)*(input_ptr+2) == 0x95) {
        current_token.type = TOK_MUL;
        input_ptr += 3;
        return;
    }
    /* ✖ U+2716 (UTF-8: 0xE2 0x9C 0x96) - heavy ballot X */
    if ((unsigned char)*input_ptr == 0xE2 && (unsigned char)*(input_ptr+1) == 0x9C && (unsigned char)*(input_ptr+2) == 0x96) {
        current_token.type = TOK_MUL;
        input_ptr += 3;
        return;
    }
    /* ⋅ U+22C5 (UTF-8: 0xE2 0x8B 0x85) - dot operator */
    if ((unsigned char)*input_ptr == 0xE2 && (unsigned char)*(input_ptr+1) == 0x8B && (unsigned char)*(input_ptr+2) == 0x85) {
        current_token.type = TOK_MUL;
        input_ptr += 3;
        return;
    }
    
    switch (*input_ptr) {
        case '+': current_token.type = TOK_PLUS; input_ptr++; return;
        case '-': current_token.type = TOK_MINUS; input_ptr++; return;
        case '*': current_token.type = TOK_MUL; input_ptr++; return;
        case '/': current_token.type = TOK_DIV; input_ptr++; return;
        case '%':
            /* MATLAB-style comment: skip to end of line */
            while (*input_ptr && *input_ptr != '\n') input_ptr++;
            /* Return TOK_END since this line is done */
            current_token.type = TOK_END;
            return;
        case '\\': current_token.type = TOK_BACKSLASH; input_ptr++; return;
        case '^': current_token.type = TOK_POW; input_ptr++; return;
        case '!':
            if (*(input_ptr+1) == '=') {
                current_token.type = TOK_NE;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_FACT;
            input_ptr++;
            return;
        case '(': current_token.type = TOK_LPAREN; input_ptr++; return;
        case ')': current_token.type = TOK_RPAREN; input_ptr++; return;
        case '=':
            /* Check for == (equality) */
            if (*(input_ptr+1) == '=') {
                current_token.type = TOK_EQUAL;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_ASSIGN;
            input_ptr++;
            return;
        case '>':
            if (*(input_ptr+1) == '=') {
                current_token.type = TOK_GE;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_GT;
            input_ptr++;
            return;
        case '<':
            if (*(input_ptr+1) == '=') {
                current_token.type = TOK_LE;
                input_ptr += 2;
                return;
            }
            if (*(input_ptr+1) == '>') {
                current_token.type = TOK_NE;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_LT;
            input_ptr++;
            return;
        case ';': current_token.type = TOK_SEMI; input_ptr++; return;
        case '&':
            if (*(input_ptr+1) == '&') {
                current_token.type = TOK_AND;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_ERROR;
            input_ptr++;
            return;
        case '|':
            if (*(input_ptr+1) == '|') {
                current_token.type = TOK_OR;
                input_ptr += 2;
                return;
            }
            current_token.type = TOK_ERROR;
            input_ptr++;
            return;
        case '~':
            if (*(input_ptr+1) == '=') {
                /* MATLAB: ~= is not-equal */
                current_token.type = TOK_NE;
                input_ptr += 2;
                return;
            }
            /* ~ alone is NOT */
            current_token.type = TOK_NOT;
            input_ptr++;
            return;
        case '[': current_token.type = TOK_LBRACKET; input_ptr++; return;
        case ']': current_token.type = TOK_RBRACKET; input_ptr++; return;
        case ':': current_token.type = TOK_COLON; input_ptr++; return;
        case ',': current_token.type = TOK_COMMA; input_ptr++; return;
        case '\'':
            /* Check if this is a string literal or transpose
             * String literal: after ( , ; [ { or at start
             * Transpose: after ) ] } number or identifier
             */
            {
                /* Look back to see context - simplified: if next char is alphanumeric or digit, it's a string */
                const char *p = input_ptr + 1;
                int is_string = 0;
                
                /* Check if looks like datetime string: starts with digit */
                if (*p >= '0' && *p <= '9') {
                    is_string = 1;
                }
                /* Also check for alphabetic start (for other string uses) */
                if (*p >= 'a' && *p <= 'z') {
                    is_string = 1;
                }
                if (*p >= 'A' && *p <= 'Z') {
                    is_string = 1;
                }
                
                if (is_string) {
                    /* Parse string literal */
                    int len = 0;
                    input_ptr++;  /* skip opening quote */
                    while (*input_ptr && *input_ptr != '\'' && len < 63) {
                        current_token.str_value[len++] = *input_ptr++;
                    }
                    current_token.str_value[len] = '\0';
                    if (*input_ptr == '\'') input_ptr++;  /* skip closing quote */
                    current_token.type = TOK_STRING;
                    return;
                }
            }
            current_token.type = TOK_TRANSPOSE;
            input_ptr++;
            return;
        case '.':
            /* Check for element-wise operators: .* ./ .^ .' .\ */
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
            if (*(input_ptr + 1) == '\\') {
                current_token.type = TOK_DOT_BACKSLASH;
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
            /* Line continuation: ... */
            if (*(input_ptr + 1) == '.' && *(input_ptr + 2) == '.') {
                /* Skip to end of line and continue */
                input_ptr += 3;
                while (*input_ptr && *input_ptr != '\n') input_ptr++;
                if (*input_ptr == '\n') input_ptr++;
                next_token();  /* Get next real token */
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
