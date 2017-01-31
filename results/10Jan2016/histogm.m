j = 1;
c = 0;
for i = 1:3200000
if(a(i) < 5e-4)
c(j) = a(i);
j++;
endif
endfor

plot(c,100)
