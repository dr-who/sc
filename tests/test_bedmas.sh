#!/bin/bash
# test_bedmas.sh - Comprehensive BEDMAS/PEMDAS order of operations tests
# Tests Brackets, Exponents, Division, Multiplication, Addition, Subtraction

SC="${1:-./bin/sc}"

PASS=0
FAIL=0

# Colors
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    CYAN='\033[0;36m'
    NC='\033[0m'
else
    GREEN=''
    RED=''
    CYAN=''
    NC=''
fi

# Test a single expression
test_eq() {
    local expr="$1"
    local expected="$2"
    local result
    
    result=$($SC "$expr" 2>/dev/null | tr -d '\n\r')
    
    if [ "$result" = "$expected" ]; then
        echo -e "  ${GREEN}PASS${NC}: $expr = $expected"
        ((PASS++))
    else
        echo -e "  ${RED}FAIL${NC}: $expr"
        echo "        got:      '$result'"
        echo "        expected: '$expected'"
        ((FAIL++))
    fi
}

echo "=============================================="
echo "  BEDMAS/PEMDAS Order of Operations Tests"
echo "=============================================="
echo ""

# =============================================
echo -e "${CYAN}=== Basic Operations ===${NC}"
# =============================================

test_eq "2+3" "5"
test_eq "5-2" "3"
test_eq "3*4" "12"
test_eq "12/4" "3"
test_eq "2^3" "8"

# =============================================
echo ""
echo -e "${CYAN}=== Unary Minus (Critical!) ===${NC}"
# =============================================

test_eq "-5" "-5"
test_eq "-5+3" "-2"
test_eq "-5-3" "-8"
test_eq "3+-5" "-2"
test_eq "3--5" "8"
test_eq "-5*2" "-10"
test_eq "2*-5" "-10"
test_eq "-5/2" "-2.5"
test_eq "-2^2" "-4"
test_eq "(-2)^2" "4"
test_eq "-(2+3)" "-5"
test_eq "--5" "5"
test_eq "---5" "-5"

# =============================================
echo ""
echo -e "${CYAN}=== Addition and Subtraction (Left to Right) ===${NC}"
# =============================================

test_eq "1+2+3" "6"
test_eq "10-5-2" "3"
test_eq "1+2-3" "0"
test_eq "10-5+2" "7"
test_eq "1-2+3-4+5" "3"
test_eq "100-50-25-10-5" "10"

# =============================================
echo ""
echo -e "${CYAN}=== Multiplication and Division (Left to Right) ===${NC}"
# =============================================

test_eq "2*3*4" "24"
test_eq "24/4/2" "3"
test_eq "24/4*2" "12"
test_eq "2*24/4" "12"
test_eq "100/10/2/5" "1"
test_eq "2*3/6*4" "4"

# =============================================
echo ""
echo -e "${CYAN}=== Multiplication/Division Before Addition/Subtraction ===${NC}"
# =============================================

test_eq "2+3*4" "14"
test_eq "3*4+2" "14"
test_eq "10-2*3" "4"
test_eq "2*3-10" "-4"
test_eq "10/2+3" "8"
test_eq "3+10/2" "8"
test_eq "10-6/2" "7"
test_eq "6/2-10" "-7"
test_eq "1+2*3+4" "11"
test_eq "1*2+3*4" "14"
test_eq "10-2*3+4/2" "6"

# =============================================
echo ""
echo -e "${CYAN}=== Exponentiation (Right Associative, Highest Precedence) ===${NC}"
# =============================================

test_eq "2^3^2" "512"        # 2^(3^2) = 2^9 = 512, NOT (2^3)^2 = 64
test_eq "2^2^3" "256"        # 2^(2^3) = 2^8 = 256
test_eq "3^2^1" "9"          # 3^(2^1) = 3^2 = 9
test_eq "2^3*4" "32"         # (2^3)*4 = 8*4 = 32
test_eq "4*2^3" "32"         # 4*(2^3) = 4*8 = 32
test_eq "2+3^2" "11"         # 2+(3^2) = 2+9 = 11
test_eq "3^2+2" "11"         # (3^2)+2 = 9+2 = 11
test_eq "2*3^2" "18"         # 2*(3^2) = 2*9 = 18
test_eq "3^2*2" "18"         # (3^2)*2 = 9*2 = 18
test_eq "10-2^3" "2"         # 10-(2^3) = 10-8 = 2
test_eq "2^3-10" "-2"        # (2^3)-10 = 8-10 = -2
test_eq "16/2^2" "4"         # 16/(2^2) = 16/4 = 4
test_eq "2^2/16" "0.25"      # (2^2)/16 = 4/16 = 0.25

# =============================================
echo ""
echo -e "${CYAN}=== Brackets Override Precedence ===${NC}"
# =============================================

test_eq "(2+3)*4" "20"
test_eq "2*(3+4)" "14"
test_eq "(10-2)*3" "24"
test_eq "10-(2*3)" "4"
test_eq "(10+2)/3" "4"
test_eq "10/(2+3)" "2"
test_eq "(2+3)^2" "25"
test_eq "2^(3+1)" "16"
test_eq "(2^3)^2" "64"
test_eq "((2+3)*4)/2" "10"
test_eq "2*((3+4)*5)" "70"
test_eq "(1+2)*(3+4)" "21"
test_eq "(10-5)*(4+2)" "30"
test_eq "((1+2)*(3+4))^2" "441"

# =============================================
echo ""
echo -e "${CYAN}=== Nested Brackets ===${NC}"
# =============================================

test_eq "((2+3))" "5"
test_eq "(((2+3)))" "5"
test_eq "((2+3)*4)" "20"
test_eq "(2*(3+4))" "14"
test_eq "((2+3)*(4+5))" "45"
test_eq "(((2+3)*(4+5))+10)" "55"
test_eq "2*((3+(4*5))+6)" "58"
test_eq "((10-2)*3)+((4+5)*2)" "42"
test_eq "(2^(3^2))" "512"
test_eq "((2^3)^2)" "64"

# =============================================
echo ""
echo -e "${CYAN}=== Complex Mixed Operations ===${NC}"
# =============================================

test_eq "1+2*3-4/2" "5"
test_eq "2^2+3*4-5" "11"
test_eq "10/2*3+4-1" "18"
test_eq "2*3+4*5" "26"
test_eq "10-2+3*4/2" "14"
test_eq "2^3+3^2" "17"
test_eq "2^3*3^2" "72"
test_eq "2^3-3^2" "-1"
test_eq "100-50+25*2-10/2" "95"
test_eq "(2+3)*(4-1)/(5-2)" "5"
test_eq "2^(1+2)*3" "24"
test_eq "(2^3+4)*(5-2)" "36"

# =============================================
echo ""
echo -e "${CYAN}=== Division Edge Cases ===${NC}"
# =============================================

test_eq "6/3/2" "1"          # (6/3)/2 = 1, NOT 6/(3/2) = 4
test_eq "24/2/3/4" "1"       # Left to right
test_eq "100/10/5/2" "1"
test_eq "1/2/2" "0.25"
test_eq "8/4/2/1" "1"
test_eq "2/1/2/1/2" "0.5"

# =============================================
echo ""
echo -e "${CYAN}=== Modulo Operator ===${NC}"
# =============================================

test_eq "10%3" "1"
test_eq "17%5" "2"
test_eq "100%7" "2"
test_eq "15%5" "0"
test_eq "23%10" "3"
test_eq "-10%3" "2"          # Floor mod: -10 - 3*floor(-10/3) = -10 - 3*(-4) = 2
test_eq "10%-3" "-2"
test_eq "mod(10, 3)" "1"
test_eq "mod(17, 5)" "2"
test_eq "7.5%2" "1.5"

# =============================================
echo ""
echo -e "${CYAN}=== Subtraction Edge Cases ===${NC}"
# =============================================

test_eq "10-5-3" "2"         # (10-5)-3 = 2, NOT 10-(5-3) = 8
test_eq "100-50-25-10" "15"
test_eq "1-2-3-4" "-8"
test_eq "10-1-1-1-1-1" "5"

# =============================================
echo ""
echo -e "${CYAN}=== Large Expression Stress Tests ===${NC}"
# =============================================

test_eq "1+2+3+4+5+6+7+8+9+10" "55"
test_eq "1*2*3*4*5" "120"
test_eq "2^2^2" "16"         # 2^(2^2) = 2^4 = 16
test_eq "1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1" "20"
test_eq "2*2*2*2*2*2*2*2*2*2" "1024"
test_eq "(1+2)*(3+4)*(5+6)" "231"
test_eq "((1+2)*3+4)*5" "65"
test_eq "1+((2+3)*4)+5" "26"
test_eq "10-9+8-7+6-5+4-3+2-1" "5"
test_eq "1*2+3*4+5*6+7*8+9*10" "190"

# =============================================
echo ""
echo -e "${CYAN}=== Fraction and Decimal Operations ===${NC}"
# =============================================

test_eq "1/2+1/2" "1"
test_eq "1/3+1/3+1/3" "1"
test_eq "0.5+0.5" "1"
test_eq "0.1+0.2" "0.3"
test_eq "1.5*2" "3"
test_eq "3/1.5" "2"
test_eq "2.5^2" "6.25"
test_eq "(1/2)^2" "0.25"

# =============================================
echo ""
echo -e "${CYAN}=== Ans Variable with Semicolons ===${NC}"
# =============================================

# These test the ans bug fix
echo "Testing ans with semicolons..."
result=$($SC "1/2;ans*2" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "1" ]; then
    echo -e "  ${GREEN}PASS${NC}: 1/2;ans*2 = 1"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 1/2;ans*2"
    echo "        got:      '$result'"
    echo "        expected: '1'"
    ((FAIL++))
fi

result=$($SC "2+3;ans*2" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "10" ]; then
    echo -e "  ${GREEN}PASS${NC}: 2+3;ans*2 = 10"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 2+3;ans*2"
    echo "        got:      '$result'"
    echo "        expected: '10'"
    ((FAIL++))
fi

result=$($SC "3^2;ans+1" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "10" ]; then
    echo -e "  ${GREEN}PASS${NC}: 3^2;ans+1 = 10"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 3^2;ans+1"
    echo "        got:      '$result'"
    echo "        expected: '10'"
    ((FAIL++))
fi

# =============================================
echo ""
echo -e "${CYAN}=== Parenthesized Division and Multiplication ===${NC}"
# =============================================

test_eq "(1/0)*(1/0)" "Inf"
test_eq "(2/1)*(3/1)" "6"
test_eq "(6/2)*(4/2)" "6"
test_eq "(10/5)*(20/4)" "10"
test_eq "(1+1)*(2+2)" "8"
test_eq "(3-1)*(5-2)" "6"

# =============================================
echo ""
echo -e "${CYAN}=== Unicode Operators × and ÷ ===${NC}"
# =============================================

test_eq "3 × 4" "12"
test_eq "12 ÷ 4" "3"
test_eq "3 + 4 × (8 - 2) - 10 ÷ 5" "25"
test_eq "2 × 3 + 4 ÷ 2" "8"
test_eq "10 ÷ 2 × 5" "25"

# =============================================
echo ""
echo -e "${CYAN}=== Equality Checking ===${NC}"
# =============================================

# Using = for equality
result=$($SC "1+3=4" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 1+3=4 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 1+3=4"
    echo "        got:      '$result'"
    echo "        expected: 'true'"
    ((FAIL++))
fi

result=$($SC "1+3=5" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "false" ]; then
    echo -e "  ${GREEN}PASS${NC}: 1+3=5 = false"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 1+3=5"
    echo "        got:      '$result'"
    echo "        expected: 'false'"
    ((FAIL++))
fi

# Using == for equality
result=$($SC "2*3==6" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 2*3==6 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 2*3==6"
    echo "        got:      '$result'"
    echo "        expected: 'true'"
    ((FAIL++))
fi

result=$($SC "2*3==7" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "false" ]; then
    echo -e "  ${GREEN}PASS${NC}: 2*3==7 = false"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 2*3==7"
    echo "        got:      '$result'"
    echo "        expected: 'false'"
    ((FAIL++))
fi

# More equality tests
result=$($SC "sqrt(4)==2" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: sqrt(4)==2 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: sqrt(4)==2"
    echo "        got:      '$result'"
    echo "        expected: 'true'"
    ((FAIL++))
fi

result=$($SC "2^10=1024" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 2^10=1024 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 2^10=1024"
    echo "        got:      '$result'"
    echo "        expected: 'true'"
    ((FAIL++))
fi

# =============================================
echo ""
echo -e "${CYAN}=== Comparison Operators ===${NC}"
# =============================================

result=$($SC "3<4" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 3<4 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 3<4"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "4<3" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "false" ]; then
    echo -e "  ${GREEN}PASS${NC}: 4<3 = false"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 4<3"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "4<=4" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 4<=4 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 4<=4"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "5>3" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 5>3 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 5>3"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "5>=5" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 5>=5 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 5>=5"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "3<>4" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "true" ]; then
    echo -e "  ${GREEN}PASS${NC}: 3<>4 = true"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 3<>4"
    echo "        got:      '$result'"
    ((FAIL++))
fi

result=$($SC "4<>4" 2>/dev/null | tr -d '\n\r')
if [ "$result" = "false" ]; then
    echo -e "  ${GREEN}PASS${NC}: 4<>4 = false"
    ((PASS++))
else
    echo -e "  ${RED}FAIL${NC}: 4<>4"
    echo "        got:      '$result'"
    ((FAIL++))
fi

# =============================================
echo ""
echo "=============================================="
echo -e "  RESULTS: ${GREEN}$PASS passed${NC}, ${RED}$FAIL failed${NC}"
TOTAL=$((PASS + FAIL))
echo "  Total: $TOTAL tests"
echo "=============================================="

# Exit with failure if any tests failed
[ $FAIL -eq 0 ]
