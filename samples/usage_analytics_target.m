%% Usage analytics with a timetable (time, customerid, usage)
% This script:
%  1) creates sample usage events
%  2) stores them in a timetable
%  3) shows the table
%  4) computes:
%       - number of unique customers
%       - new customers per month
%       - leaving customers per month (churn by last-seen month)
%
% NOTE: This is the TARGET sample. Full implementation requires:
%   - timetable() / table() data structures
%   - groupsummary() for aggregation
%   - outerjoin() for table joins
%   - String arrays for customer IDs
%
% See usage_analytics_simple.m for a working subset

clear; clc;

%% 1) Create sample event-level data
rng(7);

nEvents = 2500;

% Random timestamps across ~18 months
t0 = datetime(2024,1,1,0,0,0,'TimeZone','UTC');
t1 = datetime(2025,6,30,23,59,59,'TimeZone','UTC');
time = t0 + seconds(randi(seconds(t1-t0), nEvents, 1));

% Customer IDs
nCustPool = 400;
customerid = "C" + string(randi([1 nCustPool], nEvents, 1));

% Usage amount (e.g. GB)
usage = max(0, 10 + 50*randn(nEvents,1));          % continuous-ish
usage = round(usage, 2);

%% 2) Build a timetable
TT = timetable(time, customerid, usage);

% Sort by time for clean display & predictable grouping
TT = sortrows(TT);

%% 3) Show the table (head + basic summary)
disp("=== Usage timetable (first 15 rows) ===");
disp(TT(1:min(15,height(TT)), :));

disp("=== Summary ===");
fprintf("Total events: %d\n", height(TT));
fprintf("Time span:   %s  ..  %s\n", string(TT.time(1)), string(TT.time(end)));

%% 4) Core metrics

% 4.1 Number of customers (unique)
uniqueCustomers = unique(TT.customerid);
nCustomers = numel(uniqueCustomers);
fprintf("\nUnique customers: %d\n", nCustomers);

% 4.2 Monthly buckets for each event
monthOfEvent = dateshift(TT.time, 'start', 'month');

% 4.3 First-seen month (cohort month) per customer
custFirstTime = groupsummary(TT, "customerid", "min", "time");
custFirstMonth = dateshift(custFirstTime.min_time, 'start', 'month');

newCustomersPerMonth = table(custFirstMonth, custFirstTime.customerid, ...
    'VariableNames', ["month","customerid"]);
newCustomersPerMonth = groupsummary(newCustomersPerMonth, "month", "numel", "customerid");
newCustomersPerMonth.Properties.VariableNames(end) = "new_customers";

% 4.4 Last-seen month per customer (used for "leaving customers per month")
custLastTime = groupsummary(TT, "customerid", "max", "time");
custLastMonth = dateshift(custLastTime.max_time, 'start', 'month');

leavingCustomersPerMonth = table(custLastMonth, custLastTime.customerid, ...
    'VariableNames', ["month","customerid"]);
leavingCustomersPerMonth = groupsummary(leavingCustomersPerMonth, "month", "numel", "customerid");
leavingCustomersPerMonth.Properties.VariableNames(end) = "leaving_customers";

%% 5) Optional: Monthly active customers + usage totals (handy dashboard columns)
monthlyUsage = table(monthOfEvent, TT.customerid, TT.usage, ...
    'VariableNames', ["month","customerid","usage"]);

% Total usage per month
usagePerMonth = groupsummary(monthlyUsage, "month", "sum", "usage");
usagePerMonth.Properties.VariableNames(end) = "usage_sum";

% Active customers per month (unique customers with any usage that month)
activeCustPerMonth = groupsummary(monthlyUsage, "month", @(x) numel(unique(x)), "customerid");
activeCustPerMonth.Properties.VariableNames(end) = "active_customers";

%% 6) Combine into one monthly report (outer join so missing months don't vanish)
report = outerjoin(usagePerMonth, activeCustPerMonth, "Keys","month", "MergeKeys",true);
report = outerjoin(report, newCustomersPerMonth, "Keys","month", "MergeKeys",true);
report = outerjoin(report, leavingCustomersPerMonth, "Keys","month", "MergeKeys",true);

% Fill missing with zeros
for v = ["usage_sum","active_customers","new_customers","leaving_customers"]
    if ismember(v, report.Properties.VariableNames)
        report.(v) = fillmissing(report.(v), "constant", 0);
    end
end

report = sortrows(report, "month");

%% 7) Show the monthly report
disp(" ");
disp("=== Monthly report ===");
disp(report);

%% 8) Helper: "current customers" at end of dataset (last-seen month == last month means not left *yet*)
lastMonthInData = dateshift(TT.time(end), 'start', 'month');
stillActiveAtEnd = sum(custLastMonth == lastMonthInData);
fprintf("\nCustomers active in last month (%s): %d\n", string(lastMonthInData), stillActiveAtEnd);

%% Notes:
% - "Leaving customers per month" here means: customers whose last recorded event falls in that month.
% - If you want churn based on inactivity (e.g., no events for 60 days), say so and I'll adjust the logic.
