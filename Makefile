# Makefile for sc - Scientific Calculator
# C89 compliant for maximum portability
#
# Builds three precision versions:
#   sc      - 128-bit precision (8 x 16-bit limbs) ~38 decimal digits
#   sc-256  - 256-bit precision (16 x 16-bit limbs) ~77 decimal digits
#   sc-64   - 64-bit precision (4 x 16-bit limbs) ~19 decimal digits

CC = gcc
CFLAGS_BASE = -std=c89 -pedantic -Wall -Wextra -Werror -O2
LDFLAGS = -lm

# Directories
SRCDIR = src
BINDIR = bin
OBJDIR = obj

# Version
VERSION = 1.0.0

# Source files
SOURCES = apf.c apfx.c apfc.c apf_native.c \
          matrix.c matrix_ops.c matrix_linalg.c ml.c \
          matrix_rand.c lexer.c parser.c runtime.c format.c commands.c \
          repl.c plot.c rpn.c solver.c mathx.c stats.c tvm.c newton.c \
          orbital.c optim.c help.c func_registry.c datetime.c table.c \
          scalar_funcs.c matrix_funcs.c decomp_funcs.c saas.c main.c

# Function documentation is now all in func_registry.c

# Object directories for each precision
OBJDIR_128 = $(OBJDIR)/128
OBJDIR_256 = $(OBJDIR)/256
OBJDIR_64 = $(OBJDIR)/64

OBJECTS_128 = $(patsubst %.c,$(OBJDIR_128)/%.o,$(SOURCES))
OBJECTS_256 = $(patsubst %.c,$(OBJDIR_256)/%.o,$(SOURCES))
OBJECTS_64 = $(patsubst %.c,$(OBJDIR_64)/%.o,$(SOURCES))

# Targets
TARGET_128 = $(BINDIR)/sc
TARGET_256 = $(BINDIR)/sc-256
TARGET_64 = $(BINDIR)/sc-64

# Feature flags
FEATURES = -DHAVE_TERMIOS -DHAVE_RPN -DHAVE_SOLVER -DHAVE_STATS -DHAVE_TVM \
           -DHAVE_NEWTON -DHAVE_ORBITAL

CFLAGS_128 = $(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=8
CFLAGS_256 = $(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=16
CFLAGS_64 = $(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=4

# Installation directories
PREFIX ?= /usr
BINPREFIX ?= $(PREFIX)/bin
MANPREFIX ?= $(PREFIX)/share/man

.PHONY: all clean install uninstall test rpm srpm dist dos dos-tiny vic20 minimal sc sc-256 sc-64

# Default: build all three versions
all: $(TARGET_128) $(TARGET_256) $(TARGET_64)
	@echo ""
	@echo "Built:"
	@ls -la $(BINDIR)/sc $(BINDIR)/sc-256 $(BINDIR)/sc-64

# Individual targets
sc: $(TARGET_128)
sc-256: $(TARGET_256)
sc-64: $(TARGET_64)

# 128-bit version (default, 8 limbs)
$(TARGET_128): $(OBJECTS_128) | $(BINDIR)
	$(CC) $(OBJECTS_128) $(LDFLAGS) -o $@

$(OBJDIR_128)/%.o: $(SRCDIR)/%.c | $(OBJDIR_128)
	$(CC) $(CFLAGS_128) -c $< -o $@

$(OBJDIR_128):
	mkdir -p $(OBJDIR_128)

# 256-bit version (16 limbs)
$(TARGET_256): $(OBJECTS_256) | $(BINDIR)
	$(CC) $(OBJECTS_256) $(LDFLAGS) -o $@

$(OBJDIR_256)/%.o: $(SRCDIR)/%.c | $(OBJDIR_256)
	$(CC) $(CFLAGS_256) -c $< -o $@

$(OBJDIR_256):
	mkdir -p $(OBJDIR_256)

# 64-bit version (4 limbs)
$(TARGET_64): $(OBJECTS_64) | $(BINDIR)
	$(CC) $(OBJECTS_64) $(LDFLAGS) -o $@

$(OBJDIR_64)/%.o: $(SRCDIR)/%.c | $(OBJDIR_64)
	$(CC) $(CFLAGS_64) -c $< -o $@

$(OBJDIR_64):
	mkdir -p $(OBJDIR_64)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(OBJDIR) $(BINDIR)
	rm -rf sc-$(VERSION) sc-$(VERSION).tar.gz

install: $(TARGET_128) $(TARGET_256) $(TARGET_64)
	install -d $(DESTDIR)$(BINPREFIX)
	install -m 755 $(TARGET_128) $(DESTDIR)$(BINPREFIX)/sc
	install -m 755 $(TARGET_256) $(DESTDIR)$(BINPREFIX)/sc-256
	install -m 755 $(TARGET_64) $(DESTDIR)$(BINPREFIX)/sc-64
	install -d $(DESTDIR)$(MANPREFIX)/man1
	install -m 644 packaging/sc.1 $(DESTDIR)$(MANPREFIX)/man1/sc.1
	gzip -f $(DESTDIR)$(MANPREFIX)/man1/sc.1 2>/dev/null || true

uninstall:
	rm -f $(DESTDIR)$(BINPREFIX)/sc
	rm -f $(DESTDIR)$(BINPREFIX)/sc-256
	rm -f $(DESTDIR)$(BINPREFIX)/sc-64
	rm -f $(DESTDIR)$(MANPREFIX)/man1/sc.1.gz
	rm -f $(DESTDIR)$(MANPREFIX)/man1/sc.1

# Create distribution tarball
dist: clean
	mkdir -p sc-$(VERSION)
	cp -r src packaging tests Makefile Makefile.dos README.md COMPARISON.md LICENSE sc-$(VERSION)/
	tar czf sc-$(VERSION).tar.gz sc-$(VERSION)
	rm -rf sc-$(VERSION)

# Build RPM package
rpm: dist
	@echo "Building RPM package..."
	mkdir -p ~/rpmbuild/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
	cp sc-$(VERSION).tar.gz ~/rpmbuild/SOURCES/
	cp packaging/sc.spec ~/rpmbuild/SPECS/
	sed 's/^Version:.*/Version:        $(VERSION)/' packaging/sc.spec > ~/rpmbuild/SPECS/sc.spec
	rpmbuild -ba ~/rpmbuild/SPECS/sc.spec
	@echo ""
	@echo "RPM packages built:"
	@ls -la ~/rpmbuild/RPMS/*/*.rpm 2>/dev/null || true
	@ls -la ~/rpmbuild/SRPMS/*.rpm 2>/dev/null || true

# Build source RPM only
srpm: dist
	mkdir -p ~/rpmbuild/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
	cp sc-$(VERSION).tar.gz ~/rpmbuild/SOURCES/
	cp packaging/sc.spec ~/rpmbuild/SPECS/
	sed 's/^Version:.*/Version:        $(VERSION)/' packaging/sc.spec > ~/rpmbuild/SPECS/sc.spec
	rpmbuild -bs ~/rpmbuild/SPECS/sc.spec

# Run tests (uses 128-bit version)
test: $(TARGET_128)
	@echo "Running comprehensive test suite..."
	@./tests/test_all_functions.sh $(TARGET_128)
	@echo ""
	@echo "Running BEDMAS tests..."
	@./tests/test_bedmas.sh $(TARGET_128)
	@echo ""
	@echo "Running core tests..."
	@./tests/test_main.sh $(TARGET_128)
	@echo ""
	@echo "Running help coverage test..."
	@./tests/test_help_coverage.sh $(TARGET_128)

# Test all versions
test-all: $(TARGET_128) $(TARGET_256) $(TARGET_64)
	@echo "=== Testing sc (128-bit) ==="
	@./tests/test_bedmas.sh $(TARGET_128)
	@./tests/test_main.sh $(TARGET_128)
	@echo ""
	@echo "=== Testing sc-256 (256-bit) ==="
	@./tests/test_bedmas.sh $(TARGET_256)
	@./tests/test_main.sh $(TARGET_256)
	@echo ""
	@echo "=== Testing sc-64 (64-bit) ==="
	@./tests/test_bedmas.sh $(TARGET_64)
	@./tests/test_main.sh $(TARGET_64)

# Platform-specific builds (128-bit only for simplicity)
dos:
	@echo "Building for DOS 16-bit (SCALC_MEDIUM)..."
	$(MAKE) clean
	$(MAKE) sc CFLAGS_128="$(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=8 -DSCALC_MEDIUM"

dos-tiny:
	@echo "Building for DOS 16-bit minimal (SCALC_TINY)..."
	$(MAKE) clean
	$(MAKE) sc CFLAGS_128="$(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=8 -DSCALC_TINY"

vic20:
	@echo "Building for VIC-20 (SCALC_VIC20, 64-bit precision)..."
	@echo "Note: Requires cc65 cross-compiler for actual VIC-20 target"
	$(MAKE) clean
	$(MAKE) sc-64 CFLAGS_64="$(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=4 -DSCALC_VIC20"

minimal:
	@echo "Building minimal embedded version (SCALC_MINIMAL)..."
	$(MAKE) clean
	$(MAKE) sc-64 CFLAGS_64="$(CFLAGS_BASE) $(FEATURES) -DAP_LIMBS=4 -DSCALC_MINIMAL"

# Cross-compile for DOS using Open Watcom (if available)
dos-watcom:
	@echo "Cross-compiling for DOS with Open Watcom..."
	@command -v wcl >/dev/null 2>&1 || { echo "Error: Open Watcom not found."; exit 1; }
	mkdir -p $(BINDIR)
	wcl -0 -bt=dos -DSCALC_MEDIUM -DAP_LIMBS=8 $(addprefix $(SRCDIR)/,$(SOURCES)) -fe=$(BINDIR)/sc.exe
	wcl -0 -bt=dos -DSCALC_MEDIUM -DAP_LIMBS=16 $(addprefix $(SRCDIR)/,$(SOURCES)) -fe=$(BINDIR)/sc-256.exe
	wcl -0 -bt=dos -DSCALC_MEDIUM -DAP_LIMBS=4 $(addprefix $(SRCDIR)/,$(SOURCES)) -fe=$(BINDIR)/sc-64.exe

# Download IEEE 754 test suite (requires internet)
download-ieee754-tests:
	@echo "Downloading IEEE 754 test suite..."
	@cd tests && ./download_ieee754_suite.sh

# Dependencies (simplified)
$(OBJDIR_128)/main.o $(OBJDIR_256)/main.o $(OBJDIR_64)/main.o: $(SRCDIR)/main.c $(SRCDIR)/sc.h
$(OBJDIR_128)/parser.o $(OBJDIR_256)/parser.o $(OBJDIR_64)/parser.o: $(SRCDIR)/parser.c $(SRCDIR)/sc.h
$(OBJDIR_128)/lexer.o $(OBJDIR_256)/lexer.o $(OBJDIR_64)/lexer.o: $(SRCDIR)/lexer.c $(SRCDIR)/sc.h
$(OBJDIR_128)/apf.o $(OBJDIR_256)/apf.o $(OBJDIR_64)/apf.o: $(SRCDIR)/apf.c $(SRCDIR)/apf.h
$(OBJDIR_128)/apfx.o $(OBJDIR_256)/apfx.o $(OBJDIR_64)/apfx.o: $(SRCDIR)/apfx.c $(SRCDIR)/apfx.h $(SRCDIR)/apf.h
$(OBJDIR_128)/apfc.o $(OBJDIR_256)/apfc.o $(OBJDIR_64)/apfc.o: $(SRCDIR)/apfc.c $(SRCDIR)/apfc.h $(SRCDIR)/apf.h

sample_cluster:
	@bash samples/pca_cluster_demo.sh bin/sc


# Add func_dispatch to objects
obj/128/func_dispatch.o: src/func_dispatch.c src/func_dispatch.h
	$(CC) $(CFLAGS) $(CFLAGS_128) -c src/func_dispatch.c -o obj/128/func_dispatch.o
