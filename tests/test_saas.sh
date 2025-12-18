#!/bin/bash
# Test SaaS metrics functions

SC="${SC:-./bin/sc}"
PASS=0
FAIL=0

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

pass() { echo -e "${GREEN}✓${NC} $1"; ((PASS++)); }
fail() { echo -e "${RED}✗${NC} $1"; ((FAIL++)); }

echo "Running SaaS Metrics Tests..."
echo "=============================="
echo ""

# Generate test data if not present
if [ ! -f datasets/hourlybilling.csv ]; then
    echo "Generating test data..."
    ./gen_hourly_billing.sh datasets/hourlybilling.csv 50000 > /dev/null 2>&1
fi

# Test: Load data
echo "Testing: Data Loading"
result=$($SC << 'INPUT'
load hourlybilling
rows(data)
quit
INPUT
)
if echo "$result" | grep -q "Loaded:"; then
    pass "load hourlybilling"
else
    fail "load hourlybilling"
    echo "  Got: $result"
fi

# Test: mrr function
echo ""
echo "Testing: mrr(data)"
result=$($SC << 'INPUT'
load hourlybilling
m = mrr(data)
rows(m)
cols(m)
quit
INPUT
)
if echo "$result" | grep -qE "^[0-9]+$" && echo "$result" | grep -q "^2$"; then
    pass "mrr returns Nx2 matrix"
else
    fail "mrr returns Nx2 matrix"
    echo "  Got: $result"
fi

# Test: arpu function
echo ""
echo "Testing: arpu(data)"
result=$($SC << 'INPUT'
load hourlybilling
a = arpu(data)
rows(a) > 0
cols(a) == 2
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "arpu returns valid matrix"
else
    fail "arpu returns valid matrix"
fi

# Test: mrrbridge function
echo ""
echo "Testing: mrrbridge(data)"
result=$($SC << 'INPUT'
load hourlybilling
b = mrrbridge(data)
cols(b)
quit
INPUT
)
if echo "$result" | grep -q "^6$"; then
    pass "mrrbridge returns 6 columns"
else
    fail "mrrbridge returns 6 columns"
    echo "  Got: $result"
fi

# Test: customercount function
echo ""
echo "Testing: customercount(data)"
result=$($SC << 'INPUT'
load hourlybilling
c = customercount(data)
min(c(:,2)) >= 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "customercount returns valid counts"
else
    fail "customercount returns valid counts"
fi

# Test: newcustomers function
echo ""
echo "Testing: newcustomers(data)"
result=$($SC << 'INPUT'
load hourlybilling
n = newcustomers(data)
sum(n(:,2)) > 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "newcustomers finds new customers"
else
    fail "newcustomers finds new customers"
fi

# Test: churn function
echo ""
echo "Testing: churn(data)"
result=$($SC << 'INPUT'
load hourlybilling
ch = churn(data)
min(ch(:,2)) >= 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "churn returns non-negative values"
else
    fail "churn returns non-negative values"
fi

# Test: reactivated function
echo ""
echo "Testing: reactivated(data)"
result=$($SC << 'INPUT'
load hourlybilling
re = reactivated(data)
min(re(:,2)) >= 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "reactivated returns non-negative values"
else
    fail "reactivated returns non-negative values"
fi

# Test: churnrate function
echo ""
echo "Testing: churnrate(data)"
result=$($SC << 'INPUT'
load hourlybilling
cr = churnrate(data)
min(cr(:,2)) >= 0
max(cr(:,2)) <= 100
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "churnrate returns valid percentages"
else
    fail "churnrate returns valid percentages"
fi

# Test: nrr function
echo ""
echo "Testing: nrr(data)"
result=$($SC << 'INPUT'
load hourlybilling
nr = nrr(data)
min(nr(:,2)) > 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "nrr returns positive values"
else
    fail "nrr returns positive values"
fi

# Test: grr function
echo ""
echo "Testing: grr(data)"
result=$($SC << 'INPUT'
load hourlybilling
gr = grr(data)
max(gr(:,2)) <= 100
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "grr capped at 100%"
else
    fail "grr capped at 100%"
fi

# Test: retention function
echo ""
echo "Testing: retention(data)"
result=$($SC << 'INPUT'
load hourlybilling
R = retention(data)
R(1,1)
quit
INPUT
)
if echo "$result" | grep -q "100"; then
    pass "retention shows 100% in month 0"
else
    fail "retention shows 100% in month 0"
    echo "  Got: $result"
fi

# Test: tenure function
echo ""
echo "Testing: tenure(data)"
result=$($SC << 'INPUT'
load hourlybilling
t = tenure(data)
cols(t)
min(t(:,2)) >= 1
quit
INPUT
)
if echo "$result" | grep -q "^4$"; then
    pass "tenure returns 4 columns"
else
    fail "tenure returns 4 columns"
    echo "  Got: $result"
fi

# Test: Consistency checks
echo ""
echo "Testing: Consistency Checks"

# MRR should equal sum of customer revenues
result=$($SC << 'INPUT'
load hourlybilling
m = mrr(data)
m(1,2) > 0
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "MRR is positive for first month"
else
    fail "MRR is positive for first month"
fi

# ARPU = MRR / customers
result=$($SC << 'INPUT'
load hourlybilling
m = mrr(data)
a = arpu(data)
c = customercount(data)
abs(a(1,2) - m(1,2)/c(1,2)) < 0.01
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "ARPU = MRR / customer_count"
else
    fail "ARPU = MRR / customer_count"
fi

# New customers in month 1 = total active in month 1
result=$($SC << 'INPUT'
load hourlybilling
n = newcustomers(data)
c = customercount(data)
n(1,2) == c(1,2)
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "Month 1: new_customers = active_customers"
else
    fail "Month 1: new_customers = active_customers"
fi

# Test new SaaS functions

echo ""
echo "Testing: concentration(data)"
result=$($SC << 'INPUT'
load hourlybilling
c = concentration(data)
cols(c) == 5
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "concentration returns 5 columns"
else
    fail "concentration returns 5 columns"
fi

echo ""
echo "Testing: toprevenue(data, 5)"
result=$($SC << 'INPUT'
load hourlybilling
t = toprevenue(data, 5)
rows(t) == 5
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "toprevenue returns correct rows"
else
    fail "toprevenue returns correct rows"
fi

echo ""
echo "Testing: quickratio(data)"
result=$($SC << 'INPUT'
load hourlybilling
q = quickratio(data)
cols(q) == 2
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "quickratio returns 2 columns"
else
    fail "quickratio returns 2 columns"
fi

echo ""
echo "Testing: ltv(data)"
result=$($SC << 'INPUT'
load hourlybilling
l = ltv(data)
cols(l) == 5
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "ltv returns 5 columns"
else
    fail "ltv returns 5 columns"
fi

echo ""
echo "Testing: revchurn(data)"
result=$($SC << 'INPUT'
load hourlybilling
r = revchurn(data)
cols(r) == 4
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "revchurn returns 4 columns"
else
    fail "revchurn returns 4 columns"
fi

echo ""
echo "Testing: netchurn(data)"
result=$($SC << 'INPUT'
load hourlybilling
n = netchurn(data)
cols(n) == 2
quit
INPUT
)
if echo "$result" | grep -q "true"; then
    pass "netchurn returns 2 columns"
else
    fail "netchurn returns 2 columns"
fi

echo ""
echo "=============================="
echo "SaaS Tests: $PASS passed, $FAIL failed"

if [ $FAIL -eq 0 ]; then
    exit 0
else
    exit 1
fi
