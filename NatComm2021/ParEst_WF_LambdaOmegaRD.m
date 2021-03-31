%%%% Lambda-Omega RD Weak Formulation System Identification %%%%%
function [ksi,res,Q,q0] = ParEst_WF_LambdaOmegaRD(filename,N_d,D,if_track)
%{
INPUTS:
-------
filename                -- string containing path to matfile w/ data
N_d                     -- number of integration domains to saple
D = [Dx,Dt]             -- size of integration domain in space/time
if_track                -- bool for showing progress in command window

OUTPUTS:
--------
ksi     -- non dimensional estimated parameters
res     -- residual = mean(abs(Q*ksi-q0))
Q       -- Assembled library
q0      -- time derivative term
%}

%% INITIALIZE
% Define matfile
traj = matfile(filename);

% Get size of velocity fields
[~,Lx,Lt] = size(traj,'U_t');

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

% Initialize integration domain location
P = zeros(3,1);

% Initialize Target and Library
q0 = zeros(N_d,2);
Q = zeros(length(q0),10,2);

%% FILL TERMS
n_lib = 0;
n_track = 10;
            
% Pre-make derivatives of weight functions
p = [2,2,1]; % weight function exponents
dA000 = weight_full([0,0,0],p,x,t);
dA001 = weight_full([0,0,1],p,x,t);
dA200 = weight_full([2,0,0],p,x,t);
dA020 = weight_full([0,2,0],p,x,t);
        
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

    % Choose random point
    P(1) = randi([1,Lx-Dx]);
    P(2) = randi([1,Lx-Dx]);
    P(3) = randi([1,Lt-Dt]);

    % Indices for integration domain
    rx = P(1):(P(1) + Dx);
    ry = P(2):(P(2) + Dx);
    rt = P(3):(P(3) + Dt);

    % U in integration domain
    U = traj.U_t(ry,rx,rt); 
    V = traj.V_t(ry,rx,rt);

    %%%%% U %%%%%
    % du/dt
    B = -U.*dA001*S_t; 
    q0(n_lib,1) = trapz(x,trapz(x,trapz(t,B,3),2),1);

    % \nabla^2 u
    th1 = U.*(dA020 + dA200)*S_x^2;
    Q(n_lib,1,1) = trapz(x,trapz(x,trapz(t,th1,3),2),1);

    % u
    th2 = U.*dA000;
    Q(n_lib,2,1) = trapz(x,trapz(x,trapz(t,th2,3),2),1);

    % u^2
    th3 = U.^2.*dA000;
    Q(n_lib,3,1) = trapz(x,trapz(x,trapz(t,th3,3),2),1);

    % u^3
    th4 = U.^3.*dA000;
    Q(n_lib,4,1) = trapz(x,trapz(x,trapz(t,th4,3),2),1);

    % v
    th5 = V.*dA000;
    Q(n_lib,5,1) = trapz(x,trapz(x,trapz(t,th5,3),2),1);

    % v^2
    th6 = V.^2.*dA000;
    Q(n_lib,6,1) = trapz(x,trapz(x,trapz(t,th6,3),2),1);

    % v^3
    th7 = V.^3.*dA000;
    Q(n_lib,7,1) = trapz(x,trapz(x,trapz(t,th7,3),2),1);

    % uv
    th8 = U.*V.*dA000;
    Q(n_lib,8,1) = trapz(x,trapz(x,trapz(t,th8,3),2),1);

    % u^2v
    th9 = U.^2.*V.*dA000;
    Q(n_lib,9,1) = trapz(x,trapz(x,trapz(t,th9,3),2),1);

    % uv^2
    th10 = U.*V.^2.*dA000;
    Q(n_lib,10,1) = trapz(x,trapz(x,trapz(t,th10,3),2),1);

    %%%%% V %%%%%
    % dv/dt
    B = -V.*dA001*S_t;
    q0(n_lib,2) = trapz(x,trapz(x,trapz(t,B,3),2),1);

    % \nabla^2 v
    th1 = V.*(dA020 + dA200)*S_x^2;
    Q(n_lib,1,2) = trapz(x,trapz(x,trapz(t,th1,3),2),1);

    % u
    th2 = U.*dA000;
    Q(n_lib,2,2) = trapz(x,trapz(x,trapz(t,th2,3),2),1);

    % u^2
    th3 = U.^2.*dA000;
    Q(n_lib,3,2) = trapz(x,trapz(x,trapz(t,th3,3),2),1);

    % u^3
    th4 = U.^3.*dA000;
    Q(n_lib,4,2) = trapz(x,trapz(x,trapz(t,th4,3),2),1);

    % v
    th5 = V.*dA000;
    Q(n_lib,5,2) = trapz(x,trapz(x,trapz(t,th5,3),2),1);

    % v^2
    th6 = V.^2.*dA000;
    Q(n_lib,6,2) = trapz(x,trapz(x,trapz(t,th6,3),2),1);

    % v^3
    th7 = V.^3.*dA000;
    Q(n_lib,7,2) = trapz(x,trapz(x,trapz(t,th7,3),2),1);

    % uv
    th8 = U.*V.*dA000;
    Q(n_lib,8,2) = trapz(x,trapz(x,trapz(t,th8,3),2),1);

    % u^2v
    th9 = U.^2.*V.*dA000;
    Q(n_lib,9,2) = trapz(x,trapz(x,trapz(t,th9,3),2),1);

    % uv^2
    th10 = U.*V.^2.*dA000;
    Q(n_lib,10,2) = trapz(x,trapz(x,trapz(t,th10,3),2),1);
        
end % np
            
%% REGRESSION
% Parameters
ksi_u = SINDy(Q(:,:,1), q0(:,1));
ksi_v = SINDy(Q(:,:,2), q0(:,2));

fields = ["Laplacian";"u";"u2";"u3";"v";"v2";"v3";"uv";"u2v";"uv2"];

ksi.u = cell2struct(num2cell(ksi_u),fields);
ksi.v = cell2struct(num2cell(ksi_v),fields);

% How good was the estimation
res_u = mean(abs(q0(:,1) - Q(:,:,1)*ksi_u));
res_v = mean(abs(q0(:,2) - Q(:,:,2)*ksi_v));

res = [res_u, res_v];

end

