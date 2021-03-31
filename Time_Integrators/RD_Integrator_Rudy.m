%%%%%% Rudy's R-D Integrator %%%%%%

%{
Acquired from 

https://github.com/snagcliffs/PDE-FIND/tree/master/Datasets

through the paper 

Samuel H. Rudy, Steven L. Brunton, Joshua L. Proctor, J. Nathan Kutz
Science Advances26 Apr 2017 : e1602614 

%}

% lambda-omega reaction-diffusion system
%  u_t = lam(A) u - ome(A) v + d1*(u_xx + u_yy) = 0
%  v_t = ome(A) u + lam(A) v + d2*(v_xx + v_yy) = 0
%
%  A^2 = u^2 + v^2 and
%  lam(A) = 1 - A^2
%  ome(A) = -beta*A^2


t=0:0.05:2;
d1=0.1; d2=0.1; beta=1.0;
L=20; n=512; N=n*n;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n); y=x;
kx=(2*pi/L)*[0:(n/2-1) -n/2:-1]; ky=kx;

% INITIAL CONDITIONS

[X,Y]=meshgrid(x,y);
[KX,KY]=meshgrid(kx,ky);
K2=KX.^2+KY.^2; K22=reshape(K2,N,1);

m=2; % number of spirals

U_t = zeros(length(x),length(y),length(t));
V_t = zeros(length(x),length(y),length(t));

U_t(:,:,1)=tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
V_t(:,:,1)=tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));

% REACTION-DIFFUSION
uvt=[reshape(fft2(U_t(:,:,1)),1,N) reshape(fft2(V_t(:,:,1)),1,N)].';
[t,uvsol]=ode45(@(t,uvt) reaction_diffusion_rhs(t,uvt,[],K22,d1,d2,beta,n,N),t,uvt);


for j=1:length(t)-1
ut=reshape((uvsol(j,1:N).'),n,n);
vt=reshape((uvsol(j,(N+1):(2*N)).'),n,n);
U_t(:,:,j+1)=real(ifft2(ut));
V_t(:,:,j+1)=real(ifft2(vt));

figure(1)
pcolor(x,y,V_t(:,:,j+1)); shading interp; colormap(hot); colorbar; drawnow; 
end

dx = x(2) - x(1);
dt = t(2) - t(1);
save('reaction_diffusion_big.mat','t','x','y','U_t','V_t','dx','dt','-v7.3')

%%
load reaction_diffusion_big
pcolor(x,y,U_t(:,:,end)); shading interp; colormap(hot)


function rhs = reaction_diffusion_rhs(t,uvt,dummy,K22,d1,d2,beta,n,N)

% Calculate u and v terms
ut=reshape((uvt(1:N)),n,n);
vt=reshape((uvt((N+1):(2*N))),n,n);
u=real(ifft2(ut)); v=real(ifft2(vt));

% Reaction Terms
u3=u.^3; v3=v.^3; u2v=(u.^2).*v; uv2=u.*(v.^2);
utrhs=reshape((fft2(u-u3-uv2+beta*u2v+beta*v3)),N,1);
vtrhs=reshape((fft2(v-u2v-v3-beta*u3-beta*uv2)),N,1);

rhs=[-d1*K22.*uvt(1:N)+utrhs
-d2*K22.*uvt(N+1:end)+vtrhs];

end