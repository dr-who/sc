# Makefile for sc - Scientific Calculator
# C89 compliant for maximum portability

CC = gcc
CFLAGS = -std=c89 -pedantic -Wall -Wextra -O2
LDFLAGS = -lm

# Directories
SRCDIR = src
BINDIR = bin
OBJDIR = obj

# Version
VERSION = 1.0.0

# Source files
SOURCES = apf.c apfx.c apfc.c matrix.c matrix_ops.c matrix_linalg.c \
          matrix_rand.c lexer.c parser.c runtime.c format.c commands.c \
          repl.c plot.c rpn.c solver.c mathx.c stats.c tvm.c newton.c \
          orbital.c main.c

OBJECTS = $(patsubst %.c,$(OBJDIR)/%.o,$(SOURCES))
TARGET = $(BINDIR)/sc

# Feature flags
CFLAGS += -DHAVE_TERMIOS
CFLAGS += -DHAVE_RPN -DHAVE_SOLVER -DHAVE_STATS -DHAVE_TVM
CFLAGS += -DHAVE_NEWTON -DHAVE_ORBITAL

# Installation directories
PREFIX ?= /usr
BINPREFIX ?= $(PREFIX)/bin
MANPREFIX ?= $(PREFIX)/share/man

.PHONY: all clean install uninstall test rpm srpm dist

all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(BINDIR)
	rm -rf sc-$(VERSION) sc-$(VERSION).tar.gz

install: $(TARGET)
	install -d $(DESTDIR)$(BINPREFIX)
	install -m 755 $(TARGET) $(DESTDIR)$(BINPREFIX)/sc
	install -d $(DESTDIR)$(MANPREFIX)/man1
	install -m 644 packaging/sc.1 $(DESTDIR)$(MANPREFIX)/man1/sc.1
	gzip -f $(DESTDIR)$(MANPREFIX)/man1/sc.1 2>/dev/null || true

uninstall:
	rm -f $(DESTDIR)$(BINPREFIX)/sc
	rm -f $(DESTDIR)$(MANPREFIX)/man1/sc.1.gz
	rm -f $(DESTDIR)$(MANPREFIX)/man1/sc.1

# Create distribution tarball
dist: clean
	mkdir -p sc-$(VERSION)
	cp -r src packaging tests Makefile README.md LICENSE sc-$(VERSION)/
	tar czf sc-$(VERSION).tar.gz sc-$(VERSION)
	rm -rf sc-$(VERSION)

# Build RPM package
rpm: dist
	@echo "Building RPM package..."
	mkdir -p ~/rpmbuild/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
	cp sc-$(VERSION).tar.gz ~/rpmbuild/SOURCES/
	cp packaging/sc.spec ~/rpmbuild/SPECS/
	sed -i 's/^Version:.*/Version:        $(VERSION)/' ~/rpmbuild/SPECS/sc.spec
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
	sed -i 's/^Version:.*/Version:        $(VERSION)/' ~/rpmbuild/SPECS/sc.spec
	rpmbuild -bs ~/rpmbuild/SPECS/sc.spec

# Run tests
test: $(TARGET)
	@echo "Running test suite..."
	@./tests/test_bedmas.sh $(TARGET)
	@./tests/test_main.sh $(TARGET)

# Download IEEE 754 test suite (requires internet)
download-ieee754-tests:
	@echo "Downloading IEEE 754 test suite..."
	@cd tests && ./download_ieee754_suite.sh

# Dependencies (simplified)
$(OBJDIR)/main.o: $(SRCDIR)/main.c $(SRCDIR)/sc.h
$(OBJDIR)/parser.o: $(SRCDIR)/parser.c $(SRCDIR)/sc.h
$(OBJDIR)/lexer.o: $(SRCDIR)/lexer.c $(SRCDIR)/sc.h
$(OBJDIR)/runtime.o: $(SRCDIR)/runtime.c $(SRCDIR)/sc.h
$(OBJDIR)/format.o: $(SRCDIR)/format.c $(SRCDIR)/sc.h
$(OBJDIR)/commands.o: $(SRCDIR)/commands.c $(SRCDIR)/sc.h
$(OBJDIR)/repl.o: $(SRCDIR)/repl.c $(SRCDIR)/sc.h
$(OBJDIR)/apf.o: $(SRCDIR)/apf.c $(SRCDIR)/apf.h
$(OBJDIR)/apfx.o: $(SRCDIR)/apfx.c $(SRCDIR)/apfx.h $(SRCDIR)/apf.h
$(OBJDIR)/apfc.o: $(SRCDIR)/apfc.c $(SRCDIR)/apfc.h $(SRCDIR)/apf.h
$(OBJDIR)/matrix.o: $(SRCDIR)/matrix.c $(SRCDIR)/matrix.h
