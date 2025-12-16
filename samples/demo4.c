/* demo4.c - Boolean comparisons
 * 
 * Demonstrates comparing arbitrary precision numbers
 * and getting boolean results.
 * 
 * Compile: gcc -I../src -DAP_LIMBS=4 demo4.c ../src/apf.c -o demo4
 */

#include <stdio.h>
#include "apf.h"

int main(void)
{
    apf a, b;
    apf zero, inf, nan_val;
    apf tolerance, diff;
    
    printf("APF Boolean Comparison Demo\n");
    printf("===========================\n\n");
    
    /* Test 1: Compare 1 and 1.03 */
    apf_from_str(&a, "1");
    apf_from_str(&b, "1.03");
    
    printf("a = 1, b = 1.03\n");
    printf("a == b: %s\n", apf_eq(&a, &b) ? "true" : "false");
    printf("a != b: %s\n", !apf_eq(&a, &b) ? "true" : "false");
    printf("a < b:  %s\n", apf_cmp(&a, &b) < 0 ? "true" : "false");
    printf("a > b:  %s\n", apf_cmp(&a, &b) > 0 ? "true" : "false");
    printf("a <= b: %s\n", apf_cmp(&a, &b) <= 0 ? "true" : "false");
    printf("a >= b: %s\n", apf_cmp(&a, &b) >= 0 ? "true" : "false");
    printf("\n");
    
    /* Test 2: Compare equal values */
    apf_from_str(&a, "3.14159");
    apf_from_str(&b, "3.14159");
    
    printf("a = 3.14159, b = 3.14159\n");
    printf("a == b: %s\n", apf_eq(&a, &b) ? "true" : "false");
    printf("a != b: %s\n", !apf_eq(&a, &b) ? "true" : "false");
    printf("\n");
    
    /* Test 3: Check special values */
    apf_zero(&zero);
    apf_set_inf(&inf, 0);  /* positive infinity */
    apf_set_nan(&nan_val);
    
    printf("Special values:\n");
    printf("0 == 0:    %s\n", apf_iszero(&zero) ? "true" : "false");
    printf("inf check: %s\n", apf_isinf(&inf) ? "true" : "false");
    printf("nan check: %s\n", apf_isnan(&nan_val) ? "true" : "false");
    printf("\n");
    
    /* Test 4: Approximate equality (within tolerance) */
    apf_from_str(&a, "1.0000001");
    apf_from_str(&b, "1.0000002");
    apf_from_str(&tolerance, "0.000001");
    
    apf_sub(&diff, &a, &b);
    apf_abs(&diff, &diff);
    
    printf("a = 1.0000001, b = 1.0000002\n");
    printf("tolerance = 0.000001\n");
    printf("|a - b| <= tolerance: %s\n", 
           apf_cmp(&diff, &tolerance) <= 0 ? "true" : "false");
    
    return 0;
}
