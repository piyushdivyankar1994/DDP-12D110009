A = load("log.txt);
avg = mean(A(:,2));
var = mean(A(:,3));
var = var/100;
disp(avg)
disp(var)
plot(A(:,2));
