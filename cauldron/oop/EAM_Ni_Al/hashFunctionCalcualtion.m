A = load("pAl_r(A)_E(eV).txt");
for i = 1:5000
B(i,1) = A(i,1) - (i)*(max(A(:,1)) - min(A(:,1)))/5000 - 1;
endfor

plot(B)