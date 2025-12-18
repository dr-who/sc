t = datetime({'2024-01-15 10:30:00'; '2024-01-20 14:45:00'; '2024-02-05 09:15:00'; '2024-02-15 16:30:00'; '2024-03-01 11:00:00'})
year(t(1))
year(t(5))
month(t(1))
month(t(5))
m = [startofmonth(t(1)); startofmonth(t(2)); startofmonth(t(3)); startofmonth(t(4)); startofmonth(t(5))]
u = unique(m)
y = [10.5; 25.3; 15.8; 42.1; 8.7]
sum(y)
mean(y)
numel(u)
quit
