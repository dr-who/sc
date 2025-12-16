/* demo1.c - Basic soft-float arithmetic
 * 
 * Demonstrates using APF (Arbitrary Precision Float) library
 * for basic calculations without hardware floating point.
 * 
 * Compile: gcc -I../src -DAP_LIMBS=4 demo1.c ../src/apf.c -o demo1
 */

#include <stdio.h>
#include "apf.h"

int main(void)
{
    apf a, b, result;
    char buf[64];
    
    /* Initialize two numbers */
    apf_from_str(&a, "3.14159");
    apf_from_str(&b, "2.71828");
    
    printf("APF Soft-Float Demo\n");
    printf("===================\n\n");
    
    /* Display the values */
    apf_to_str(buf, sizeof(buf), &a, 10);
    printf("a = %s\n", buf);
    
    apf_to_str(buf, sizeof(buf), &b, 10);
    printf("b = %s\n", buf);
    
    /* Multiply */
    apf_mul(&result, &a, &b);
    apf_to_str(buf, sizeof(buf), &result, 10);
    printf("\na * b = %s\n", buf);
    
    /* Add */
    apf_add(&result, &a, &b);
    apf_to_str(buf, sizeof(buf), &result, 10);
    printf("a + b = %s\n", buf);
    
    /* Subtract */
    apf_sub(&result, &a, &b);
    apf_to_str(buf, sizeof(buf), &result, 10);
    printf("a - b = %s\n", buf);
    
    /* Divide */
    apf_div(&result, &a, &b);
    apf_to_str(buf, sizeof(buf), &result, 10);
    printf("a / b = %s\n", buf);
    
    /* Square root */
    apf_sqrt(&result, &a);
    apf_to_str(buf, sizeof(buf), &result, 10);
    printf("sqrt(a) = %s\n", buf);
    
    return 0;
}
