%%%% Integral Method Q2D-NS Parameter Estimation %%%%
%% ----------------------------------------------------------------------
%                           PREAMBLE
% -----------------------------------------------------------------------
%{
INPUTS 
------
filename        : string for path to -v7.3 .mat file containing data
    - U_t,V_t   : non-dim velocity fields
    - dx        : non-dim spatial grid spacing
    - dt_mks    : dimensional temporal grid spacing (seconds)
    - Ts_mks    : time-scale in seconds
N_d             : number of integration domains
N_h = [Q R S] 
    - S         : max harmonic in time (>=1)
    - Q,R       : max harmonic in x,y resp. (>=0)
D = [fr, Dt] 
    - fr        : spatial-area of int domain as fraction of full area
    - Dt        : time-length of integration domain in grid points
track           : boolean to track library filling in cmd window

OUTPUTS
-------
ksi : estimated parameters in non-dimensinoal units
res : residual of estimation (utility unclear)
Q   : matrix containing integral evaluations of gov eq terms
q0  : vector containing integral evaluation of time-derivative term
P   : 3xN_d matrix containing integration domain locations 

SECTIONS
--------
Initialize
    -size of local domain
    -how many harmonics to sample
    -any tracking variables
    -library and target

Construct Library
    -advection term
    -laplacian term
    -rayleigh term
    -time-derivative term

Regression
    -invert library onto target or use SINDy
    -calculate norm of b-\Theta\xi

Reminders
    -disp to-do list

FUNCTIONS
---------
Construct 3D Weight Function
    - use General Liebniz Rule to find the k^th oder product rule of base
        weight and harmonics
    - use meshgrid to combine 3 1D composite fields into final 3D weight
        function

Construct 1D Base Weight Function
    - find base poly coefficients based on exponent m
    - convert poly coefficients based on derivative order k
    - compute full polynomial

Construct 1D Weight harmonics
    - Calculate derivatives of sin(pi*q*x) in 1D

%}
%%
function [ksi,res,Q,q0,P] = ParEst_WF_Kolmo(filename,N_d,N_h,D,track)
%% ----------------------------------------------------------------------
%                           INITIALIZE
% -----------------------------------------------------------------------
% Define matfile
traj = matfile(filename);

% Get size of velocity fields
[Ly,Lx,Lt] = size(traj,'U_t');

% Grid densities
dt = traj.dt;
dx = traj.dx; 

% Size of local domain
fr = D(1);
Dx = round(fr*Lx); % sizes of the integration domain
Dy = round(fr*Lx);
Dt = D(2);

% Sampling scheme
P = zeros(3,N_d); % Starting corner of integration domain
P(1,:) = randi([1,Lx-Dx],N_d,1);
P(2,:) = randi([1,Ly-Dy],N_d,1);
P(3,:) = randi([1,Lt-Dt],N_d,1);

% Define derivative conversions
S_x = 2/(dx*Dx);
S_y = 2/(dx*Dy);
S_t = 2/(dt*Dt);

% Create subdomain
x = linspace(-1,1,Dx); 
y = linspace(-1,1,Dy);
t = linspace(-1,1,Dt);

% Initialize Target and Library
q0 = zeros(N_d*(N_h(1)+1)*(N_h(2)+1)*(N_h(3)+1),1);
q1 = q0; q2 = q0; q3 = q0;

%% ----------------------------------------------------------------------
%                           CONSTRUCT LIBRARY
% -----------------------------------------------------------------------
n_lib = 0;
n_track = 10;
p = [3,3,0]; % weight function exponents

for np = 1:N_d  

    % Indices for integration domain
    rx = P(1,np):(P(1,np)+Dx-1);
    ry = P(2,np):(P(2,np)+Dy-1);
    rt = P(3,np):(P(3,np)+Dt-1);

    % Velocity fields on integration domain
    U = traj.U_t(ry,rx,rt);
    V = traj.V_t(ry,rx,rt);

    for q = 0:N_h(1)
    for r = 0:N_h(2)
    for s = 1:N_h(3)
                
            n_lib = n_lib + 1;

            if track && n_lib == n_track
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
            
            % Make wave numbers global for use in weight_full()
            q_vec = [q,r,s];

            % Make derivatives of windowing functions
            dA100 = weight_full_harm([1,0,0],p,q_vec,x,y,t); 
            dA010 = weight_full_harm([0,1,0],p,q_vec,x,y,t);           
            dA101 = weight_full_harm([1,0,1],p,q_vec,x,y,t); 
            dA011 = weight_full_harm([0,1,1],p,q_vec,x,y,t);            
            dA110 = weight_full_harm([1,1,0],p,q_vec,x,y,t);
            dA020 = weight_full_harm([0,2,0],p,q_vec,x,y,t); 
            dA200 = weight_full_harm([2,0,0],p,q_vec,x,y,t);                         
            dA210 = weight_full_harm([2,1,0],p,q_vec,x,y,t); 
            dA120 = weight_full_harm([1,2,0],p,q_vec,x,y,t);              
            dA030 = weight_full_harm([0,3,0],p,q_vec,x,y,t); 
            dA300 = weight_full_harm([3,0,0],p,q_vec,x,y,t);                                   
            % --> the 3 digits correspond to derivatives, NOT to harmonics

            % Time-derivative
            b = (V.*dA101*S_x - U.*dA011*S_y)*S_t;
            q0(n_lib,1) = trapz(x,trapz(y,trapz(t,b,3)));

            % Advection Term (incompressible)
            th1 = U.*V.*(dA020*S_y^2 - dA200*S_x^2) + ...
                  (U.^2 - V.^2).*dA110*S_x*S_y;
            q1(n_lib,1) = trapz(x,trapz(y,trapz(t,th1,3)));   

            % Laplacian Term
            th2 = U.*(dA210*S_x^2*S_y + dA030*S_y^3) - ...
                V.*(dA300*S_x^3 + dA120*S_x*S_y^2);
            q2(n_lib,1) = trapz(x,trapz(y,trapz(t,th2,3)));

            % Rayleigh Term
            th3 = V.*dA100*S_x - U.*dA010*S_y;
            q3(n_lib,1) = trapz(x,trapz(y,trapz(t,th3,3)));           

    end % s
    end % r
    end % q
end % np

% Assemble Library
Q = [q1 q2 q3];

%% ----------------------------------------------------------------------
%                               REGRESSION
% -----------------------------------------------------------------------
% Compute parameters
ksi = Q \ q0;

% Compute residual
res = mean(abs(q0 - Q*ksi));

end
