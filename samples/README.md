# SCALC Sample Programs

These sample programs demonstrate how to use the APF (Arbitrary Precision Float) 
and APFC (Arbitrary Precision Float Complex) libraries standalone.

## Demo Programs

### demo1.c - Basic Soft-Float Arithmetic
Demonstrates basic arithmetic operations (add, subtract, multiply, divide, sqrt)
using the APF library without any hardware floating point.

```bash
gcc -I../src -DAP_LIMBS=4 demo1.c ../src/apf.c -o demo1
./demo1
```

### demo2.c - Complex Number Operations
Demonstrates complex number arithmetic using the APFC library.

```bash
gcc -I../src -DAP_LIMBS=4 demo2.c ../src/apf.c ../src/apfc.c ../src/apfx.c -o demo2
./demo2
```

### demo3.c - Expression Parsing
Shows how to parse and evaluate mathematical expressions. This example
requires linking with the full scalc library.

The easiest way to evaluate expressions is through the calculator:
```bash
echo "1+(3*4)" | ../bin/sc
```

### demo4.c - Boolean Comparisons
Demonstrates comparing arbitrary precision numbers and getting boolean results.

```bash
gcc -I../src -DAP_LIMBS=4 demo4.c ../src/apf.c -o demo4
./demo4
```

## Building All Demos

```bash
make
```

## Precision Settings

The `-DAP_LIMBS=N` flag controls precision:
- `AP_LIMBS=4` - 64-bit precision (~19 decimal digits)
- `AP_LIMBS=8` - 128-bit precision (~38 decimal digits, default)
- `AP_LIMBS=16` - 256-bit precision (~77 decimal digits)

## Using in Your Own Projects

To use the APF library in your own project:

1. Copy `apf.c` and `apf.h` to your project
2. Include `apf.h` in your source files
3. Compile with `-DAP_LIMBS=N` to set precision
4. No external dependencies required (pure C89)

For complex numbers, also include `apfc.c` and `apfc.h`.
