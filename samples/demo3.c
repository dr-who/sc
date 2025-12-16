/* demo3.c - Expression parsing
 * 
 * Demonstrates parsing and evaluating mathematical expressions
 * using the scalc expression parser.
 * 
 * This is a more complex example that uses the full parser.
 * 
 * Compile with full calculator (simpler approach):
 *   Use scalc interactively or with: echo "1+(3*4)" | ./sc
 * 
 * For embedding, see the eval_expr_with_var() function in runtime.c
 */

#include <stdio.h>
#include <string.h>
#include "config.h"
#include "apf.h"
#include "apfc.h"

/* External declarations from scalc */
extern const char *input_ptr;
extern void next_token(void);
extern int parse_expr(apfc *result);
extern void init_user_funcs(void);
extern void init_variables(void);
extern void apfx_init(void);

int main(void)
{
    const char *expressions[] = {
        "1 + (3 * 4)",
        "2^10",
        "sqrt(2)",
        "sin(3.14159/2)",
        "(1+2)*(3+4)",
        "10!",
        NULL
    };
    int i;
    
    printf("Expression Parser Demo\n");
    printf("======================\n\n");
    
    /* Initialize the calculator */
    apfx_init();
    init_user_funcs();
    init_variables();
    
    /* Evaluate each expression */
    for (i = 0; expressions[i] != NULL; i++) {
        apfc result;
        char buf[64];
        
        input_ptr = expressions[i];
        next_token();
        
        if (parse_expr(&result)) {
            apf_to_str(buf, sizeof(buf), &result.re, 10);
            printf("%s = %s\n", expressions[i], buf);
        } else {
            printf("%s = (error)\n", expressions[i]);
        }
    }
    
    return 0;
}
