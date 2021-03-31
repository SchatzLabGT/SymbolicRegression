function W = weight_full(k,p,x,t)
%{
Assemble the 1D weight functions into the full weight

k = [kx,ky,kt]: order of derivative(s)
p = [px,py,pt]: exponents of weight polynomials
%}

if length(k) == 3
    wx = weight_poly(x,p(1),k(1));
    wy = weight_poly(x,p(2),k(2));   
    wt = weight_poly(t,p(3),k(3));
    [wX,wY,wT] = meshgrid(wx,wy,wt);
    W = wX.*wY.*wT;
elseif length(k) == 2
    wx = weight_poly(x,p(1),k(1));   
    wt = weight_poly(t,p(2),k(2));
    [wT,wX] = meshgrid(wt,wx);
    W = wX.*wT;
end

end