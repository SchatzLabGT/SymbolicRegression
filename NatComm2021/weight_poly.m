function p = weight_poly(x,m,k)
%{
Polynomial piece of weighting function used to satisfy BC

A = d^k/dx^k[ (x^2 - 1)^m ]

x: independent variable
m: power of base function
k: order of derivative
%}
a = zeros(m*2 + 1,1); % initial coefficent vector
for l = 0:m
    a(2*l+1) = (-1)^(m-l)*nchoosek(m,l); % set polynomial coefficients
end 

c = zeros(2*m+1,1); % final coefficient vector
for n = 0:(2*m - k)
    c(n+1) = a(n+1+k)*factorial(n+k)/factorial(n);
end

p = 0;
for n = 0:(2*m-k)
    p = p + c(n+1)*x.^n; % final windowing function
end

end