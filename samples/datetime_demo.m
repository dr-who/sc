t = datetime({'2024-01-15 10:30:00'; '2024-01-20 14:45:00'; '2024-02-05 09:15:00'; '2024-02-15 16:30:00'; '2024-03-01 11:00:00'})
year(t(1))
month(t(1))
day(t(1))
hour(t(1))
minute(t(1))
second(t(1))
startofmonth(t(1))
startofyear(t(1))
startofday(t(1))
datenum(2024, 6, 15, 12, 0, 0)
m = [startofmonth(t(1)); startofmonth(t(2)); startofmonth(t(3)); startofmonth(t(4)); startofmonth(t(5))]
unique(m)
y = [10.5; 25.3; 15.8; 42.1; 8.7]
sum(y)
mean(y)
max(y)
min(y)
hours(2)
minutes(30)
days(1)
t2 = datetime({'2024-01-15 10:00:00'; '2024-01-15 11:30:00'; '2024-01-15 13:15:00'; '2024-01-15 14:45:00'})
v = [100; 150; 120; 180]
retime t2 v hourly linear
rt_times
rt_data
quit
