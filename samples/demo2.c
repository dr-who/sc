/* demo2.c - Complex number arithmetic
 * 
 * Demonstrates using APFC (Arbitrary Precision Float Complex) library
 * for complex number calculations without hardware floating point.
 * 
 * Compile: gcc -I../src -DAP_LIMBS=4 demo2.c ../src/apf.c ../src/apfc.c -o demo2
 */

#include <stdio.h>
#include "apf.h"
#include "apfc.h"

static void print_complex(const char *name, const apfc *z)
{
    char re_buf[64], im_buf[64];
    apf im_abs;
    
    apf_to_str(re_buf, sizeof(re_buf), &z->re, 10);
    apf_to_str(im_buf, sizeof(im_buf), &z->im, 10);
    
    if (apf_iszero(&z->im)) {
        printf("%s = %s\n", name, re_buf);
    } else if (apf_iszero(&z->re)) {
        printf("%s = %si\n", name, im_buf);
    } else if (z->im.sign) {
        /* Get absolute value for display */
        apf_abs(&im_abs, &z->im);
        apf_to_str(im_buf, sizeof(im_buf), &im_abs, 10);
        printf("%s = %s - %si\n", name, re_buf, im_buf);
    } else {
        printf("%s = %s + %si\n", name, re_buf, im_buf);
    }
}

int main(void)
{
    apfc a, b, result;
    
    printf("APFC Complex Number Demo\n");
    printf("========================\n\n");
    
    /* Create complex numbers: a = 3 + 4i, b = 1 - 2i */
    apf_from_int(&a.re, 3);
    apf_from_int(&a.im, 4);
    
    apf_from_int(&b.re, 1);
    apf_from_int(&b.im, -2);
    
    print_complex("a", &a);
    print_complex("b", &b);
    printf("\n");
    
    /* Addition */
    apfc_add(&result, &a, &b);
    print_complex("a + b", &result);
    
    /* Subtraction */
    apfc_sub(&result, &a, &b);
    print_complex("a - b", &result);
    
    /* Multiplication: (3+4i)(1-2i) = 3-6i+4i-8i² = 3-2i+8 = 11-2i */
    apfc_mul(&result, &a, &b);
    print_complex("a * b", &result);
    
    /* Division */
    apfc_div(&result, &a, &b);
    print_complex("a / b", &result);
    
    /* Conjugate */
    apfc_conj(&result, &a);
    print_complex("conj(a)", &result);
    
    /* Create i directly: 0 + 1i */
    printf("\n");
    apf_zero(&a.re);
    apf_from_int(&a.im, 1);
    print_complex("i", &a);
    
    /* i^2 should be -1 */
    apfc_mul(&result, &a, &a);
    print_complex("i * i", &result);
    
    return 0;
}
