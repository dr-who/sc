Name:           sc
Version:        1.0.0
Release:        1%{?dist}
Summary:        Arbitrary precision scientific calculator

License:        MIT
URL:            https://github.com/dr-who/sc
Source0:        %{name}-%{version}.tar.gz

BuildRequires:  gcc
BuildRequires:  make

# Disable debugsource package (not needed for simple C project)
%global debug_package %{nil}

%description
sc is an interactive arbitrary precision scientific calculator supporting
complex numbers, matrices, and a wide range of mathematical functions.
It uses 128-bit quad precision floating point internally for high accuracy.

Features:
- 128-bit quad precision (approximately 34 significant digits)
- Complex number arithmetic
- Matrix operations and linear algebra
- Scientific functions (trig, hyperbolic, exponential, logarithmic)
- Multiple number formats (decimal, hex, binary, octal, Roman numerals)
- IEEE 754-2008 compliant arithmetic (proper Inf/NaN handling)
- Comparison operators (<, <=, >, >=, ==, <>, ~=)
- Boolean operators (and, or, xor, not)
- Approximately equal operator (~=) for tolerance comparisons
- Primality testing with isprime()
- Shell-friendly exit codes (0=true/non-zero, 1=false/zero)
- Interactive REPL with command history
- Command-line and pipe modes
- HP-style RPN mode
- Time Value of Money (TVM) calculations
- Statistical functions (mean, stddev, variance, etc.)
- Equation solving (Newton-Raphson, quadratic formula)
- Combinatorics (nPr, nCr, factorial)

%prep
%setup -q

%build
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT PREFIX=/usr

%files
%license LICENSE
%doc README.md
%{_bindir}/sc
%{_bindir}/sc-256
%{_bindir}/sc-64
%{_mandir}/man1/sc.1*

%changelog
* Sun Dec 15 2024 Stu <stu@example.com> - 1.0.0-1
- Initial RPM release
- 128-bit quad precision arithmetic
- Complex number support
- Matrix operations
- IEEE 754-2008 compliant (Inf/NaN)
- RPN mode
- TVM calculations
- Statistical functions
- Boolean operators (and, or, xor, not)
- Comparison operators (<, <=, >, >=, ==, <>)
- Approximately equal (~=) for tolerance comparisons
- Primality testing with isprime()
- Shell-friendly exit codes for scripting
