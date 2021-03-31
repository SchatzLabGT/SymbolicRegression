%% -----------------------------------------------------------------------
%                       FULL 3D WINDOWING FUNCTION
% ------------------------------------------------------------------------
function W = weight_full_harm(k,p,q,x,y,t)
%{
Combine base windows and harmonics in 1D, then assemble into 3D
k = [kx,ky,kt]: order of derivative(s)
p = [px,py,pt] : weight function exponents
q = [q,r,s] : harmonic numbers
%}

% General Liebniz Rule 
wx = 0;
for n = 0:k(1)
    f = weight_poly(x,p(1),k(1)-n); 
    g = weight_harm(x,n,q(1));
    wx = wx + nchoosek(k(1),n)*f.*g;  
end
wy = 0;
for n = 0:k(2)
    f = weight_poly(y,p(2),k(2)-n);
    g = weight_harm(y,n,q(2));
    wy = wy + nchoosek(k(2),n)*f.*g;
end
wt = 0;
for n = 0:k(3)
    f = weight_poly(t,p(3),k(3)-n);
    g = weight_harm(t,n,q(3));
    wt = wt + nchoosek(k(3),n)*f.*g;
end

[wX,wY,wT] = meshgrid(wx,wy,wt);        % Make 1D components 3D

W = wX.*wY.*wT;                         % combine all components

end