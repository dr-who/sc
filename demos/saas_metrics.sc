% SaaS Metrics Demo
% Run commands one at a time to see results

% Step 1: Load data
load hourlybilling

% Step 2: Check data size
rows(data)
cols(data)

% Step 3: MRR - Monthly Recurring Revenue
m = mrr(data)

% Step 4: ARPU - Average Revenue Per User
a = arpu(data)
mean(a(:,2))

% Step 5: MRR Bridge (New, Expansion, Contraction, Churn, Net)
b = mrrbridge(data)

% Step 6: Customer counts
c = customercount(data)

% Step 7: New customers
n = newcustomers(data)
sum(n(:,2))

% Step 8: Churned customers
ch = churn(data)
sum(ch(:,2))

% Step 9: Churn rate (%)
cr = churnrate(data)
mean(cr(2:rows(cr),2))

% Step 10: Net Revenue Retention (%)
nr = nrr(data)
mean(nr(2:rows(nr),2))

% Step 11: Gross Revenue Retention (%)
gr = grr(data)
mean(gr(2:rows(gr),2))

% Step 12: Customer tenure
t = tenure(data)
mean(t(:,2))
max(t(:,2))

% Step 13: Cohort retention matrix (first 6x6)
R = retention(data)
R(1:6, 1:6)
