function W = weight_full(k,F,o)
global var
%{
Assemble 1D weight functions into (a derivative of) the full weight
Third entries in arguments optional (only used for 3D PDEs)

k = [kx,ky,kt]: order of derivative(s)
F = [Fx,Fy,Ft]: power of (x^2-1) term
w = [ox,oy,ot]: frequency of sin/cos
%}
s = pi;
ifsinx = mod(o(1),2);
o1 = s*floor((o(1)+1)/2);
for i=0:k(1)
    if k(1)==0
        coeff = 1;
    else
        coeff = nchoosek(k(1),i);
    end
    wxi(:,i+1) = coeff*weight_poly(var.x,F(1),k(1)-i);
    sign = (-1)^(floor((i+1-ifsinx)/2));
    if mod(i+ifsinx,2)==0
       trxi(:,i+1)=sign*o1^i*cos(o1*var.x);
    else
       trxi(:,i+1)=sign*o1^i*sin(o1*var.x); 
    end
end
if length(k)==3
    ifsiny = mod(o(2),2);
    o2 = s*floor((o(2)+1)/2);
    for i=0:k(2)
        if k(2)==0
            coeff = 1;
        else
            coeff = nchoosek(k(2),i);
        end
        wyi(:,i+1) = coeff*weight_poly(var.x,F(2),k(2)-i);
        sign = (-1)^(floor((i+1-ifsiny)/2));
        if mod(i+ifsiny,2)==0
           tryi(:,i+1)=sign*o2^i*cos(o2*var.x);
        else
           tryi(:,i+1)=sign*o2^i*sin(o2*var.x); 
        end
    end
    ifsint = mod(o(3),2);
    o3 = s*floor((o(3)+1)/2);
    for i=0:k(3)
        if k(3)==0
            coeff = 1;
        else
            coeff = nchoosek(k(3),i);
        end
        wti(:,i+1) = coeff*weight_poly(var.t,F(3),k(3)-i);
        sign = (-1)^(floor((i+1-ifsint)/2));
        if mod(i+ifsint,2)==0
           trti(:,i+1)=sign*o3^i*cos(o3*var.t);
        else
           trti(:,i+1)=sign*o3^i*sin(o3*var.t); 
        end
    end
    wx = sum(wxi.*trxi,2);
    wy = sum(wyi.*tryi,2);
    wt = sum(wti.*trti,2);
    [wX,wY,wT] = meshgrid(wx,wy,wt);
    W = wX.*wY.*wT;
else
    ifsint = mod(o(2),2);
    o2 = s*floor((o(2)+1)/2);
    for i=0:k(2)
        if k(2)==0
            coeff = 1;
        else
            coeff = nchoosek(k(2),i);
        end
        wti(:,i+1) = coeff*weight_poly(var.t,F(2),k(2)-i);
        sign = (-1)^(floor((i+1-ifsint)/2));
        if mod(i+ifsint,2)==0
           trti(:,i+1)=sign*o2^i*cos(o2*var.t);
        else
           trti(:,i+1)=sign*o2^i*sin(o2*var.t); 
        end
    end
    wx = sum(wxi.*trxi,2);
    wt = sum(wti.*trti,2);
    [wX,wT] = meshgrid(wx,wt);
    W = wX'.*wT';
end

end
