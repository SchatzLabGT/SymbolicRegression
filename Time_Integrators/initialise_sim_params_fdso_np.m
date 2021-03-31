function initialise_sim_params_fdso_np(I_mks,dt_mks,n,monopole_sense)
%% The code was originally written by Balachandra Suri
%{

INPUT:
------
I_mks           - Driving current in Amps
dt_mks          - Integration time-step 
n               - grid points per unit length (grid spacing: dx = 1/n)
monopole_sense  - forcing profile orientation

OUTPUT:
-------
Directory containing necessary files for integrator to function properly

To recreate PNAS data set:
I_mks           = 19.5/1000
dt_mks          = 1/32
n               = 40 (and then sub-sampled to 10 post creation)
monopole_sense  = 'ccw'

NOTE: Using n = 20 subsampled down to n = 10 should yield comparable
results as the n = 40 sub-sampled used in PNAS paper. Use n = 20 for faster
results.

%}
%%
disp('Initialising Simulation Parameters')

% Load forcing profile
if n == 20
    load('forcing_profile_20.mat','F0x')
elseif n == 40
    load('forcing_profile_40.mat','F0x')
else
    error(['No forcing profile for the grid spacing n = ',num2str(n)])
end

% Create directory for desired forcing current
dir_name = sprintf('%05.2fmA',I_mks*1000);
if (~exist(dir_name,'dir'))
mkdir(dir_name);
end
cd(dir_name);


%% Experimental Parameters Specification in MKS
% Constants
inch = 2.54/100;                                %since we use inch a lot in our length scales

% Fluid Properties
rho_mks = 959.4;                               % Depth Averaged Density
mu_mks = 0.0031;                               % effective dynamic viscosity
nu_mks = mu_mks/rho_mks;                        % kinematic Viscosity
rayfric_mks = 0.0645;                            % Rayleigh friction
advpf = 0.826;                                  % advection prefactor

% here we specitfy the number of magnet pairs
Nmp = 7;

% Basic length scales of the domain
% width of the aluminium box;
Lx_mks = 6*inch + 1*inch;
% height of aluminium box
Ly_mks = Nmp*inch + 2*inch;
% length scale of forcing
L0y_mks = 1/2*inch;

%% here we set the sense of monopole for global circulation
if(strcmp(monopole_sense,'ccw'))
    sense = +1;
elseif(strcmp(monopole_sense,'cw'))
    sense = -1;
end
%% define length, velocity and time scales
% F_mks = get_forcing_info_np(I_mks,sense,'amp');
% A = F_mks/rho_mks;
A = compute_A(I_mks,rho_mks);
Ls_mks = L0y_mks;
Us_mks = sqrt(A*Ls_mks);
Ts_mks = Ls_mks/Us_mks;

%% Non Dimensionalization
% Length scales

% length of box
Lx = Lx_mks/Ls_mks;
% height of box
Ly = Ly_mks/Ls_mks;

% Non dimensional simulation parameters
% viscous Reynolds number
Re_nu = Us_mks*Ls_mks/nu_mks;
% friction Reynolds Number
Re_rf = Us_mks/Ls_mks/rayfric_mks;
% non-dimensional time
dt = dt_mks/Ts_mks;

%% Defining a Grid

% Define the x and y axes.
% the number of intervals we define a unit length into
n0 = n;
% bins along x direction
nx = round(Lx*n0);
% bins along y direction
ny = round(Ly*n0);
% x axis coordinates
x = linspace(0,Lx,nx+1);
% y axis coordinates
y = linspace(0,Ly,ny+1);
dx = Lx/nx;
dy = Ly/ny;

[X,Y] = meshgrid(x,y);

% Define the coordinates of grid , i.e xu, xv etc.

xu = linspace(0,Lx,nx+1);
yu = linspace(-dy/2,Ly+dy/2,ny+2);

xv = linspace(-dx/2,Lx+dx/2,nx+2);
yv = linspace(0,Ly,ny+1);

xp = linspace(-dx/2,Lx+dx/2,nx+2);
yp = linspace(-dy/2,Ly+dy/2,ny+2);

% Define U, V and P matrices: Using interior points only

U = zeros(ny,nx-1);
V = zeros(ny-1,nx);
P = zeros(ny,nx);

% Here we save the size of each matrix
sizeU = size(U);
sizeV = size(V);
sizeP = size(P);

% Define the extended matrices: Including boundary points

Ue = zeros(ny+2,nx+1);
Ve = zeros(ny+1,nx+2);

%% Laplacian Definitions
LapU = compute_Laplacian(ny,nx-1,-2,-3,dx,dy);
LapV = compute_Laplacian(ny-1,nx,-3,-2,dx,dy);
LapP = compute_Laplacian(ny,nx,-1,-1,dx,dy);

%% Crank Nicholson LHS matrices

LHSU = (1+1/2*dt/Re_rf)*speye(size(LapU)) - 1/2*dt/Re_nu*LapU;
LHSV = (1+1/2*dt/Re_rf)*speye(size(LapV)) - 1/2*dt/Re_nu*LapV;
LHSP = -LapP; % the negative sign is necessary for cholsky decomposition
LHSP(1,1) = 3/2*LHSP(1,1);


%% Rearrange the rows and columns to maximize sparse zeros

permU = symamd(LHSU);
permV = symamd(LHSV);
permP = symamd(LHSP);

%% We perform the Cholskey decomposirion here

UTcholU = chol(LHSU(permU,permU));
LTcholU = UTcholU';
UTcholV = chol(LHSV(permV,permV));
LTcholV = UTcholV';
UTcholP = chol(LHSP(permP,permP));
LTcholP = UTcholP';

save('sim_params.mat','Nmp','I_mks','dt_mks','advpf','rayfric_mks','nu_mks','dt',...
    'A','Ls_mks','Us_mks','Ts_mks','Re_nu','Re_rf','monopole_sense');
save('grid_params.mat','Lx','Ly','Ls_mks','n0','nx','ny','dx','dy','dt','sizeU','sizeV','sizeP',...
    'x','y','X','Y','xu','xv','yu','yv','xp','yp');
save('differential_ops.mat','U','V','P','Ue','Ve','permU','permV','permP','LapU','LapV','LapP',...
    'UTcholU','UTcholV','UTcholP','LTcholU','LTcholV','LTcholP');

save('differential_ops.mat','F0x','-append');
fprintf(' Simulation Parameters\n nu = %05.02f cSt. \n rayfric = %05.3f 1/s \n advpf = %05.3f',nu_mks*10E5,rayfric_mks,advpf);
pause(3);
clear
cd ..
end
%% In this function we compute the laplacian matrices
function [Laplacian] = compute_Laplacian(nrows,ncols,xbcval,ybcval,dx,dy)
% here we compute the x derivate part of lapacian
% we create a tri diagonal sparse matrix with coefficients for one column
d2Qdx2 = spdiags(ones(ncols,1)*[1, -2, 1],[-1,0,1],ncols,ncols);
% we use implicit boundary conditions here, and hence this requires setting
% the end of the diagonal elements to match boundary conditions
d2Qdx2(1,1) = xbcval; d2Qdx2(ncols,ncols) = xbcval;
% here we create the full matrix for as many rows and columns
d2Qdx2 = (1/dx^2)*kron(d2Qdx2,speye(nrows));

% here we compute the y derivate part of lapacian
% we create a tri diagonal sparse matrix with coefficients for one row
d2Qdy2 = spdiags(ones(nrows,1)*[1, -2, 1],[-1,0,1],nrows,nrows);
d2Qdy2(1,1) = ybcval; d2Qdy2(nrows,nrows) = ybcval;
d2Qdy2 = (1/dy^2)*kron(speye(ncols),d2Qdy2);

Laplacian = d2Qdx2 + d2Qdy2;
end
%% RAVI's computation of A
function A = compute_A(I_mks,rho_mks)
inch = 2.54/100;
hc = 0.003;
hd = 0.003;
B0 = .2515-25.71*0.9e-3;
B1 = -25.71;
Le = (6.0*inch+1.0*inch);
B_mean = (B0 + B1*hd + B1*hc/2)*hc/(hc+hd);
J = I_mks/Le/hc;
A = B_mean*J*(1.0/rho_mks);
end
