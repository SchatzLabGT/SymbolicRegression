%%%% Integral Method Kuramoto-Sivashinsky Parameter Estimation %%%%

%% PSEUDO-CODE
%{

INPUTS 
------
filename : string for path to -v7.3 .mat file containing data
N_d : number of integration domains
D = [fr, Dt] : size of integration domain

OUTPUTS
-------
par : estimated parameters in physical units
ksi : estimated parameters in non-dimensinoal units
res : residual of estimation (utility unclear)

Initialize
    -size of local domain
    -any tracking variables
    -library and target

Fill terms
    -advection term
    -laplacian term
    -rayleigh term
    -time-derivative term

Regression
    -invert library onto target
    -calculate norm of b-\Theta\xi

Reminders
    -disp to-do list

weight_full(k)
    -inputs: k = [kx,ky,kt], order of derivative(s)
    -output: 3D array corresponding to local sub-domain

weight_poly(x,m,k)
    -inputs: absiscae, polynomial order, derivative order
    -output: additive windowing polynomial

%}
%%
function [ksi,res,Q,q0] = ParEst_wf_KS(filename,N_d,D,if_track,if_symreg)

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
Dx = D(1);
Dt = D(2);

% Create subdomain
x = linspace(-1,1,Dx+1); 
t = linspace(-1,1,Dt+1);

% Define variable conversions
S_x = 2/(dx*Dx);
S_t = 2/(dt*Dt);

% Time sampling scheme
P = zeros(2,N_d);
P(1,:) = randi([1,Lx-Dx],N_d,1);
P(2,:) = randi([1,Lt-Dt],N_d,1);

% Initialize Target and Library
q0 = zeros(N_d,1);
Q = zeros(length(q0),8);

%% FILL LIBRARY AND TARGET
n_lib = 0;
n_track = 10;
            
% Pre-make derivatives of windowing functions
p = [4,3];                          % weight function exponents
dA00 = weight_full([0,0],p,x,t);
dA01 = weight_full([0,1],p,x,t);
dA10 = weight_full([1,0],p,x,t);
dA20 = weight_full([2,0],p,x,t);
dA40 = weight_full([4,0],p,x,t);
dA30 = weight_full([3,0],p,x,t);

for np = 1:N_d  

    n_lib = n_lib + 1;

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
    rx = P(1,np):(P(1,np)+Dx);
    rt = P(2,np):(P(2,np)+Dt);

    % Velocity fields on integration domain
    U = traj.uu(rx,rt); 

    % Target
    B = U.*dA01*S_t; 
    q0(n_lib,1) = trapz(x,trapz(t,B,2),1);

    % Advection Term 
    th1 = -(1/2)*U.^2.*dA10*S_x; 
    Q(n_lib,1) = trapz(x,trapz(t,th1,2),1);

    % Laplacian Term
    th2 = U.*dA20*S_x^2;
    Q(n_lib,2) = trapz(x,trapz(t,th2,2),1);

    % Biharmonic Term
    th3 = U.*dA40*S_x^4; 
    Q(n_lib,3) = trapz(x,trapz(t,th3,2),1);

    % Linear Term
    th4 = U.*dA00;
    Q(n_lib,4) = trapz(x,trapz(t,th4,2),1);

    % First Order Derivative
    th5 = U.*dA10*S_x;
    Q(n_lib,5) = trapz(x,trapz(t,th5,2),1);

    % Third Order Derivative
    th6 = U.*dA30*S_x^3;
    Q(n_lib,6) = trapz(x,trapz(t,th6,2),1);

    % Quadratic Term
    th7 = U.^2.*dA00;
    Q(n_lib,7) = trapz(x,trapz(t,th7,2),1);

    % Cubic Term
    th8 = U.^3.*dA00;
    Q(n_lib,8) = trapz(x,trapz(t,th8,2),1);

end
            
%% REGRESSION
% Parameters
if if_symreg
    ksi = SINDy(Q, q0);     % sparsify library
    res = norm(q0 - Q*ksi);
else
    ksi = Q(:,1:3) \ q0;
    res = norm(q0 - Q(1:3)*ksi);
end
    
end
