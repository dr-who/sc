#!/bin/bash
# Generate realistic SaaS hourly billing data with complex customer behaviors
# Usage: ./gen_hourly_billing.sh [output.csv] [num_rows]

OUTPUT="${1:-datasets/hourlybilling.csv}"
TARGET_ROWS="${2:-100000}"

# Ensure datasets directory exists
mkdir -p "$(dirname "$OUTPUT")"

awk -v target="$TARGET_ROWS" -v seed="$RANDOM" '
BEGIN {
    srand(seed)
    
    # Configuration
    START_TS = 1704067200  # Jan 1, 2024 00:00:00 UTC
    HOURS_PER_MONTH = 730
    NUM_MONTHS = 24
    
    # Customer segments with different behaviors
    # Segment 1: Enterprise (high value, low churn)
    # Segment 2: Mid-market (medium value, medium churn)  
    # Segment 3: SMB (low value, high churn)
    # Segment 4: Startup (variable, high growth or churn)
    
    NUM_ENTERPRISE = 8
    NUM_MIDMARKET = 25
    NUM_SMB = 60
    NUM_STARTUP = 20
    TOTAL_CUSTOMERS = NUM_ENTERPRISE + NUM_MIDMARKET + NUM_SMB + NUM_STARTUP
    
    # Initialize customer properties
    cust_id = 0
    
    # Enterprise customers (IDs 1-8)
    for (i = 1; i <= NUM_ENTERPRISE; i++) {
        cust_id++
        segment[cust_id] = "enterprise"
        base_size[cust_id] = 50 + rand() * 150  # 50-200 TB
        base_rate[cust_id] = 0.08 + rand() * 0.04  # $0.08-0.12/TB/hr
        start_month[cust_id] = 1 + int(rand() * 3)  # Start in months 1-3
        churn_prob[cust_id] = 0.01  # 1% monthly churn
        expand_prob[cust_id] = 0.15  # 15% chance to expand
        contract_prob[cust_id] = 0.05  # 5% chance to contract
        reactivate_prob[cust_id] = 0.3  # 30% reactivation if churned
    }
    
    # Mid-market customers (IDs 9-33)
    for (i = 1; i <= NUM_MIDMARKET; i++) {
        cust_id++
        segment[cust_id] = "midmarket"
        base_size[cust_id] = 10 + rand() * 40  # 10-50 TB
        base_rate[cust_id] = 0.10 + rand() * 0.05  # $0.10-0.15/TB/hr
        start_month[cust_id] = 1 + int(rand() * 6)  # Start in months 1-6
        churn_prob[cust_id] = 0.03  # 3% monthly churn
        expand_prob[cust_id] = 0.12
        contract_prob[cust_id] = 0.08
        reactivate_prob[cust_id] = 0.2
    }
    
    # SMB customers (IDs 34-93)
    for (i = 1; i <= NUM_SMB; i++) {
        cust_id++
        segment[cust_id] = "smb"
        base_size[cust_id] = 1 + rand() * 9  # 1-10 TB
        base_rate[cust_id] = 0.12 + rand() * 0.08  # $0.12-0.20/TB/hr
        start_month[cust_id] = 1 + int(rand() * 12)  # Start throughout year 1
        churn_prob[cust_id] = 0.06  # 6% monthly churn
        expand_prob[cust_id] = 0.08
        contract_prob[cust_id] = 0.10
        reactivate_prob[cust_id] = 0.15
    }
    
    # Startup customers (IDs 94-113) - high variance
    for (i = 1; i <= NUM_STARTUP; i++) {
        cust_id++
        segment[cust_id] = "startup"
        base_size[cust_id] = 0.5 + rand() * 5  # 0.5-5.5 TB
        base_rate[cust_id] = 0.15 + rand() * 0.10  # $0.15-0.25/TB/hr
        start_month[cust_id] = 1 + int(rand() * 18)  # Start anytime
        churn_prob[cust_id] = 0.08  # 8% monthly churn
        expand_prob[cust_id] = 0.25  # High growth potential
        contract_prob[cust_id] = 0.12
        reactivate_prob[cust_id] = 0.1
    }
    
    # Track customer state per month
    for (c = 1; c <= TOTAL_CUSTOMERS; c++) {
        current_size[c] = base_size[c]
        is_active[c] = 0
        churned_month[c] = 0
        months_since_churn[c] = 0
    }
    
    # Print header
    print "timestamp,customerid,sizetb,billing"
    
    rows_written = 0
    hours_per_row = 1
    
    # Calculate how many hours to skip to hit target rows
    total_possible_hours = NUM_MONTHS * HOURS_PER_MONTH
    # Estimate active customers on average
    avg_active = TOTAL_CUSTOMERS * 0.6
    total_possible_rows = total_possible_hours * avg_active
    if (total_possible_rows > target) {
        hours_per_row = int(total_possible_rows / target) + 1
    }
    
    # Generate data month by month
    for (month = 1; month <= NUM_MONTHS; month++) {
        month_start_ts = START_TS + (month - 1) * HOURS_PER_MONTH * 3600
        
        # Update customer states for this month
        for (c = 1; c <= TOTAL_CUSTOMERS; c++) {
            # Check if customer should start
            if (!is_active[c] && churned_month[c] == 0 && month >= start_month[c]) {
                is_active[c] = 1
                current_size[c] = base_size[c] * (0.8 + rand() * 0.4)  # 80-120% of base
            }
            
            # Check for reactivation
            if (!is_active[c] && churned_month[c] > 0) {
                months_since_churn[c]++
                if (rand() < reactivate_prob[c] && months_since_churn[c] >= 2) {
                    is_active[c] = 1
                    churned_month[c] = 0
                    months_since_churn[c] = 0
                    current_size[c] = base_size[c] * (0.5 + rand() * 0.5)  # Come back smaller
                }
            }
            
            # Check for churn (only if active)
            if (is_active[c] && rand() < churn_prob[c]) {
                is_active[c] = 0
                churned_month[c] = month
                months_since_churn[c] = 0
                continue
            }
            
            # Check for expansion/contraction (only if active)
            if (is_active[c]) {
                if (rand() < expand_prob[c]) {
                    growth = 1.1 + rand() * 0.4  # 10-50% growth
                    current_size[c] = current_size[c] * growth
                } else if (rand() < contract_prob[c]) {
                    shrink = 0.6 + rand() * 0.3  # 10-40% shrink
                    current_size[c] = current_size[c] * shrink
                    if (current_size[c] < 0.1) current_size[c] = 0.1
                }
                
                # Add some monthly noise
                noise = 0.95 + rand() * 0.1  # +/- 5%
                current_size[c] = current_size[c] * noise
            }
        }
        
        # Generate hourly records for this month
        for (hour = 0; hour < HOURS_PER_MONTH; hour += hours_per_row) {
            ts = month_start_ts + hour * 3600
            
            for (c = 1; c <= TOTAL_CUSTOMERS; c++) {
                if (is_active[c]) {
                    # Add hourly variation
                    hourly_size = current_size[c] * (0.9 + rand() * 0.2)
                    billing = hourly_size * base_rate[c]
                    
                    printf "%d,%d,%.4f,%.2f\n", ts, c, hourly_size, billing
                    rows_written++
                    
                    if (rows_written >= target) {
                        exit
                    }
                }
            }
        }
        
        if (rows_written >= target) {
            exit
        }
    }
}
' > "$OUTPUT"

# Show summary
LINES=$(wc -l < "$OUTPUT")
echo "Generated $OUTPUT with $((LINES - 1)) data rows"
echo ""
head -5 "$OUTPUT"
