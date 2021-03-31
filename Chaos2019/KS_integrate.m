% kursiv.m - solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
%   u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0,32*pi]
%   computation is based on v = fft(u), so linear term is diagonal
%   compare p27.m in Trefethen, "Spectral Methods in MATLAB", SIAM 2000
%   AK Kassam and LN Trefethen, July 2002

function [uu,tt] = KS_integrate(Nx,Nt,par,if_plot,if_save)
% Spatial grid and initial condition:
Lx = 32*pi;                         % physical length of domain
x = Lx*(1:Nx)'/Nx;
dx = x(2) - x(1);
u = cos(x/16).*(1+sin(x/16));       % initial condition
v = fft(u);
% Precompute various ETDRK4 scalar quantities:
h = 5/100;                          % time-step
k = [0:Nx/2-1 0 -Nx/2+1:-1]'/16;    % wave numbers
L = par(1)*k.^2 - par(2)*k.^4;      % Fourier multipliers
E = exp(h*L); E2 = exp(h*L/2);
M = 16;                             % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);        % roots of unity
LR = h*L(:,ones(M,1)) + r(ones(Nx,1),:);
Q  = h*real(mean(           (exp(LR/2)-1)./LR              ,2));
f1 = h*real(mean(  (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3  ,2));
f2 = h*real(mean(    (2+LR+exp(LR).*(-2+LR))./LR.^3    ,2));
f3 = h*real(mean(  (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3  ,2));

% Main time-stepping loop:
uu = u; tt = 0;
tskip = 50;
tmax = 550;                         % time-length
dt = 0.25;
nmax = round(tmax/h); nplt = floor((tmax/Nt)/h*dt);
g = -0.5i*k;
for n = 1:nmax
    t = n*h;
    Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v + Q.*Nv;
    Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v + Q.*Na;
    Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = g.*fft(real(ifft(c)).^2);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    if mod(n,nplt)==0
        u = real(ifft(v));
        uu = [uu,u]; tt = [tt,t];
    end
end
%dt = tt(2)-tt(1);
uu = uu(:,floor(tskip/dt)+1:end-1);
tt = tt(:,floor(tskip/dt)+1:end-1);
uu = uu(:,1:2000);
tt = tt(:,1:2000);
% Plot results:
if(if_plot)
    surf(tt,x,uu), shading interp, lighting phong, axis tight
    view([-90 90]), colormap(autumn); set(gca,'zlim',[-5 50])
    light('color',[1 1 0],'position',[-1,2,2])
    material([0.30 0.60 0.60 40.00 1.00]);
end

% Save results
if(if_save)
    
    [Nx,Nt] = size(uu);
    save_name = ['KS_',num2str(par(1)),'_',num2str(par(2)),'.mat'];
    save(save_name,'tt','uu','par','Nx','Nt','dx','dt','-v7.3')
    
end

end
