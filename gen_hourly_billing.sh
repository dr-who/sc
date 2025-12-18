#!/bin/bash
# Generate hourly billing test data for SaaS metrics
# Creates CSV with columns: timestamp (unix), customerid (numeric), sizetb, billing

OUTPUT="${1:-datasets/hourlybilling.csv}"
NUM_ROWS="${2:-100000}"

echo "Generating $NUM_ROWS rows of hourly billing data..."

awk -v rows="$NUM_ROWS" -v seed="$RANDOM" '
BEGIN {
    srand(seed)
    
    # 25 customers
    num_custs = 25
    for (i = 1; i <= num_custs; i++) {
        if (i <= 10) {
            cust_start[i] = 0
        } else {
            cust_start[i] = int(rand() * 2000)
        }
        cust_size[i] = 0.5 + rand() * 5
        cust_price[i] = (80 + rand() * 40) / 720
        cust_churn[i] = (0.02 + rand() * 0.06) / 720
        cust_rejoin[i] = (0.3 + rand() * 0.3) / 720
        cust_growth[i] = (-0.01 + rand() * 0.03) / 720
        cust_active[i] = 1
    }
    
    # Header
    print "timestamp,customerid,sizetb,billing"
    
    row_count = 0
    hour = 0
    base_ts = 1704067200  # 2024-01-01 00:00:00 UTC
    
    while (row_count < rows && hour < 17520) {  # ~2 years max
        ts = base_ts + hour * 3600
        
        for (i = 1; i <= num_custs; i++) {
            if (hour < cust_start[i]) continue
            
            if (cust_active[i]) {
                if (rand() < cust_churn[i]) {
                    cust_active[i] = 0
                    continue
                }
            } else {
                if (rand() < cust_rejoin[i]) {
                    cust_active[i] = 1
                    cust_size[i] = cust_size[i] * (0.6 + rand() * 0.8)
                }
                continue
            }
            
            cust_size[i] = cust_size[i] * (1 + cust_growth[i])
            if (cust_size[i] < 0.1) cust_size[i] = 0.1
            if (cust_size[i] > 20) cust_size[i] = 20
            
            size = cust_size[i] * (0.85 + rand() * 0.3)
            billing = size * cust_price[i]
            
            printf "%d,%d,%.4f,%.2f\n", ts, i, size, billing
            row_count++
            
            if (row_count >= rows) exit
        }
        
        hour++
    }
}' > "$OUTPUT"

echo "Generated: $OUTPUT"
wc -l "$OUTPUT"
head -5 "$OUTPUT"
