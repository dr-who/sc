# Makefile for scalc sample programs

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall -Wextra -O2 -I../src -DAP_LIMBS=4
SRCDIR = ../src

# Simple demos that only need apf.c
DEMOS_SIMPLE = demo1 demo4

# Demos that need apfc.c
DEMOS_COMPLEX = demo2

# Full parser demos (need more sources)
# demo3 requires the full calculator

all: $(DEMOS_SIMPLE) $(DEMOS_COMPLEX)

demo1: demo1.c $(SRCDIR)/apf.c
	$(CC) $(CFLAGS) -o $@ demo1.c $(SRCDIR)/apf.c

demo2: demo2.c $(SRCDIR)/apf.c $(SRCDIR)/apfc.c $(SRCDIR)/apfx.c
	$(CC) $(CFLAGS) -o $@ demo2.c $(SRCDIR)/apf.c $(SRCDIR)/apfc.c $(SRCDIR)/apfx.c

demo4: demo4.c $(SRCDIR)/apf.c
	$(CC) $(CFLAGS) -o $@ demo4.c $(SRCDIR)/apf.c

clean:
	rm -f demo1 demo2 demo4

test: all
	@echo "=== Demo 1: Basic Arithmetic ==="
	./demo1
	@echo ""
	@echo "=== Demo 2: Complex Numbers ==="
	./demo2
	@echo ""
	@echo "=== Demo 4: Boolean Comparisons ==="
	./demo4

.PHONY: all clean test
