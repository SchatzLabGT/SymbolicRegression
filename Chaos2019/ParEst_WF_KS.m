%%%% Weak Formulation of Kuramoto-Sivashinsky Parameter Estimation %%%%

%{

INPUTS 
------
filename : string for path to -v7.3 .mat file containing data
N_d : number of integration domains
F : weight function powers [alpha, beta]
wts : frequencies of weight functions used in x and t
      (0 : frequency 0, 1-2 : sin/cos of frequency pi, etc.)
D = [Dx, Dt] : size of integration domain
if_track : enable output during running
sig : size of artificial noise
seed : random seed to fix realization of noise and integration domains

OUTPUTS
-------
ksi : estimated parameters in physical units
res : residual of estimation
Q : library matrix

%}
%%
function [ksi,res,Q] = ParEst_WF_KS(filename,N_d,F,wts,D,if_track,sig,seed)

%% INITIALIZE
% Define matfile
traj = matfile(filename);

% Get size of velocity fields
[Lx,Lt] = size(traj,'uu');

% Grid densities
dt = traj.dt;
dx = traj.dx; 

% % Size of local domain
clearvars -global
global var
%var.Dx = D(1);
%var.Dt = D(2);
var.Dx = D(1)-1;
var.Dt = D(2)-1;

var.harX = length(wts{1});
var.harT = length(wts{2});

% Create subdomain
var.x = linspace(-1,1,var.Dx+1); 
var.t = linspace(-1,1,var.Dt+1);

% Define variable conversions
S_x = 2/(dx*var.Dx);
S_t = 2/(dt*var.Dt);

% Set random number generator state
if nargin>7
    rng(seed);
end

% Time sampling scheme
P = zeros(2,N_d);
P(1,:) = randi([1,Lx-var.Dx],N_d,1);
P(2,:) = randi([1,Lt-var.Dt],N_d,1);

% Initialize Target and Library
Q = zeros(N_d*var.harX*var.harT,10); %4

%% FILL TERMS

n_lib = 0;
n_track = 10;

U_full = traj.uu;

% Add noise to u
if sig > 0
    U_full = U_full + sig*std(U_full(:))*randn(size(U_full)); % Gaussian
    %U_full = U_full + sig*std(U_full(:))*sqrt(3)*(-1+2.*rand(size(U_full))); % uniform
end

for i = 1:length(wts{1})
for j = 1:length(wts{2})

end
end

% Compute each library row
for i = 1:length(wts{1})
for j = 1:length(wts{2})
    % Pre-compute weight functions and their derivatives
    ox = wts{1}(i);
    ot = wts{2}(j);
    dA00{i}{j} = weight_full([0,0],F,[ox,ot]);
    dA01{i}{j} = weight_full([0,1],F,[ox,ot]);
    dA10{i}{j} = weight_full([1,0],F,[ox,ot]);
    dA20{i}{j} = weight_full([2,0],F,[ox,ot]);
    dA40{i}{j} = weight_full([4,0],F,[ox,ot]);
    dA30{i}{j} = weight_full([3,0],F,[ox,ot]);

for np = 1:N_d  

n_lib = n_lib + 1;

% Output
if if_track && n_lib == n_track
    if n_lib < 100
        disp(['Library Row # : ',num2str(n_lib)])
        n_track = n_track + 10;
    elseif n_lib < 1000
        disp(['Library Row # : ',num2str(n_lib)])
        n_track = n_track + 100;
    else
        disp(['Library Row # : ',num2str(n_lib)])
        n_track = n_track + 1000;
    end
end

% Indices for integration domain
rx = P(1,np):(P(1,np)+var.Dx);
rt = P(2,np):(P(2,np)+var.Dt);

% u field on integration domain
U = U_full(rx,rt);

% Time Derivative Term
B = -U.*dA01{i}{j}*S_t;

Q(n_lib,1) = trapz(var.x,trapz(var.t,B,2),1);

% Advection Term 
th1 = -(1/2)*U.^2.*dA10{i}{j}*S_x; 
Q(n_lib,2) = trapz(var.x,trapz(var.t,th1,2),1);

% Laplacian Term
th2 = U.*dA20{i}{j}*S_x^2;
Q(n_lib,3) = trapz(var.x,trapz(var.t,th2,2),1);

% Biharmonic Term
th3 = U.*dA40{i}{j}*S_x^4; 
Q(n_lib,4) = trapz(var.x,trapz(var.t,th3,2),1);

% Linear Term
th4 = U.*dA00{i}{j};
Q(n_lib,5) = trapz(var.x,trapz(var.t,th4,2),1);

% First Order Derivative
th5 = -U.*dA10{i}{j}*S_x;
Q(n_lib,6) = trapz(var.x,trapz(var.t,th5,2),1);

% Third Order Derivative
th6 = -U.*dA30{i}{j}*S_x^3;
Q(n_lib,7) = trapz(var.x,trapz(var.t,th6,2),1);

% Quadratic Term
th7 = U.^2.*dA00{i}{j};
Q(n_lib,8) = trapz(var.x,trapz(var.t,th7,2),1);

% Cubic Term
th8 = U.^3.*dA00{i}{j};
Q(n_lib,9) = trapz(var.x,trapz(var.t,th8,2),1);

% Constant Term
th9 = dA00{i}{j};
Q(n_lib,10) = trapz(var.x,trapz(var.t,th9,2),1);

end
            
end
end

% Sparsification parameter
gamma = 1.4;

% Model parameters
ksi = find_coeffs(Q(:,:),gamma);

% Normalize with respect to du/dt term
if ksi(1)~=0
    ksi = ksi/ksi(1);
end

% Compute residual (relative)
res = norm(Q(:,:)*ksi)/norm(Q(:,:))/norm(ksi);

end