%% -----------------------------------------------------------------------
%                     WEIGHT FUNCTION HARMONICS
% ------------------------------------------------------------------------
function g = weight_harm(x,k,q)
% g_j = sin(q*pi*x) 
if mod(k,2)                                 % odd k
    n = (1/2)*(k - 1);
    g = (-1)^n*(pi*q)^k*cos(pi*q*x);
else                                        % even k
    n = k/2; 
    g = (-1)^n*(pi*q)^k*sin(pi*q*x);
end
if q == 0 && k == 0
    g = ones(size(x));
elseif q == 0 && k > 0
    g = zeros(size(x));
end
end