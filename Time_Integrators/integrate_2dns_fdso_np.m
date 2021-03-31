function integrate_2dns_fdso_np(I_mks,ti_mks,tf_mks,if_continue,save_for_rec,save_rate,if_plot)
%{
The code was written by Balachandra Suri
Crank Nicholson for the Linear term
2nd order Adam's Bashforth for the NonLinear term
%}

%% Input Arguments
% ti_mks is the initial time
% tf_mks is how long we want to simulate for, usually 1000 seconds
% if_continue 0: not continuing 1:continuing from previous data
% save_for_rec: 0: don't save along the way 1: save 
% save_rate: times per second to save
% if_plot: 0:don't plot 1: plot results

% SAMPLE: integrate_2dns_fdso_np(2000,2100,1,1,4,0)
% -> start from t=2000s from previous run, and continue until 2100s, saving
% every 1/4 s.

%% Some matlab commands to clear display etc...
clc
close all
pause(2);

%% Check directory for desired current

dir_name = sprintf('%05.2fmA',I_mks*1000);
if (~exist(dir_name,'dir'))
    error(['Directory, ' dir_name ', does not exist. Initialize with' ...
           ' initialise_sim_params_fdso_np.m'])
end
cd(dir_name)

%% Here we initialise some variables
U = []; V = []; P = []; F0x = []; X = []; Y = []; omega_t = []; dt = [];
U_t = []; V_t = []; P_t = [];

%% here we load sim_params.mat
if(if_continue == 1)
    file_name = input('enter file name: ','s');
    load(file_name,'U','V','P');
    % we set the arbitraty constant on pressure to zero
    P = P-mean(P(:));
elseif(if_continue == 0)
    display('starting a simulation from random initial conditions');
    load('differential_ops.mat', 'U','V','P');
    P = P-mean(P(:));
end
load('sim_params.mat','Ts_mks','dt_mks','advpf','Ls_mks','Us_mks','Re_nu','Re_rf','I_mks');
dt = dt_mks/Ts_mks;
load('grid_params.mat','dx','dy','X','Y','Ly','Lx');
x = 0:dx:Lx;
y = 0:dy:Ly;
load('differential_ops.mat','Ue','Ve','permU','permV','permP',...
    'UTcholU','UTcholV','UTcholP','LTcholU','LTcholV','LTcholP','F0x');

%% Here we set the time variable for integration
nt = round((tf_mks-ti_mks)/dt_mks);
tf = (ti_mks+nt*dt_mks)/Ts_mks;

%% Here we initialise some counter variables
gamma_cnt = 0;
% rate at which we want toe vorticity field to be saved
omega_save_fps = save_rate;
% total number of omega that are to be saved
omega_save_n = round(tf_mks-ti_mks)*omega_save_fps;
% one in how many omega do we intend to save. in terms of numbers
omega_save_dn = round(1/dt_mks/omega_save_fps);
% lets us know how many omegas have been saved
omega_save_cnt = 0;
% marks when is the next instant to save omega
omega_save_next = omega_save_dn;

%% Allocate Additional memory for temporary variables
u = zeros(size(U(:)));
v = zeros(size(V(:)));
phi = zeros(size(P(:)));
Phie = zeros(size(P)+[2,2]);

[lU,wU] = size(U);
[lV,wV] = size(V);
[lP,wP] = size(P);

if(save_for_rec)
    omega_t = zeros(round(Ly/dy)+1,round(Lx/dx)+1,omega_save_n);
    U_t = zeros(lU,wU,omega_save_n);
    V_t = zeros(lV,wV,omega_save_n);
    P_t = zeros(lP,wP,omega_save_n);
end
rms_u = zeros(nt,1);
rms_v = zeros(nt,1);
err_u = zeros(nt,1);
err_v = zeros(nt,1);
gamma_arr = zeros(nt,1);

%% This is the loop where we step forward in time

for k = 1:nt
    
    %% Extend the Matrices to include boundary points too
    % we read the interior values into the extended matrices
    Ue(2:end-1,2:end-1) = U;
    Ve(2:end-1,2:end-1) = V;
    
    % we now impose the boundary conditions on Ue using boundary points
    Ue(:,1) = 0;
    Ue(:,end) = 0;
    Ue(1,:) = -Ue(2,:);
    Ue(end,:) = -Ue(end-1,:);
    
    % we now impose the boundary conditions on Ve using boundary points
    Ve(:,1) = -Ve(:,2);
    Ve(:,end) = -Ve(:,end-1);
    Ve(1,:) = 0;
    Ve(end,:) = 0;
    
    %% here we choose the transaition parameter for hybrid difference
    % set gamma = 0 for central difference
    % gamma = min(1.2*dt*max([max(abs(U(:)))/dx,max(abs(V(:)))/dy]),1);
    gamma = 0;
    if(gamma==0 && gamma_cnt == 0)
        disp('Using a central difference scheme for NLT');
        gamma_cnt = 1;
    elseif(gamma_cnt == 0)
        disp('Using a hybrid difference scheme for NLT');
        gamma_cnt = 1;
    end
    gamma_arr(k) = gamma;
    
    %% Here we compute the Nonlinear terms for  U,V
    Uax = (Ue(:,2:end)+Ue(:,1:end-1))/2;
    Uay = (Ue(2:end,:)+Ue(1:end-1,:))/2;
    Vax = (Ve(:,2:end)+Ve(:,1:end-1))/2;
    Vay = (Ve(2:end,:)+Ve(1:end-1,:))/2;
    
    % here we compute the x and y gradients
    dUx = (Ue(:,2:end)-Ue(:,1:end-1))/2;
    dUy = (Ue(2:end,:)-Ue(1:end-1,:))/2;
    dVx = (Ve(:,2:end)-Ve(:,1:end-1))/2;
    dVy = (Ve(2:end,:)-Ve(1:end-1,:))/2;
    
    % here we compute the nonlinear terms for U equation
    %U advecting U along x
    UadvU = Uax.^2 - gamma*abs(Uax).*dUx;
    dU2dx = (UadvU(:,2:end) - UadvU(:,1:end-1))/dx;
    %V advecting U along y
    VadvU = Vax.*Uay - gamma*abs(Vax).*dUy;
    dVUdy = (VadvU(2:end,:)-VadvU(1:end-1,:))/dy;
    
    % add both parts of nltx using interior points only
    nltU = advpf*(dU2dx(2:end-1,:) +  dVUdy(:,2:end-1));
    
    % here we compute the nonlinear terms for V equation
    % V advecting V along y
    VadvV = Vay.^2 - gamma*abs(Vay).*dVy;
    dV2dy = (VadvV(2:end,:) - VadvV(1:end-1,:))/dy;
    % U advectng V along x
    UadvV = Uay.*Vax - gamma*abs(Uay).*dVx;
    dUVdx = (UadvV(:,2:end)-UadvV(:,1:end-1))/dx;
    
    % add both parts of nltV using only interior points
    nltV = advpf*(dV2dy(:,2:end-1) + dUVdx(2:end-1,:));
    
    % This would run only for the first Adam Bashforth iteration
    % since the estimate for the previous step is not available
    if(~exist('nltU_prev','var'))
        nltU_prev = nltU;
    end
    if(~exist('nltV_prev','var'))
        nltV_prev = nltV;
    end
    
    %% Here we compute the Linear term, Pressure term and force of RHS
    % This is the linear term for the U, V components
    
    % The following lines compute the terms of laplacian
    d2Udx2 = (Ue(2:end-1,3:end)+Ue(2:end-1,1:end-2)-2*Ue(2:end-1,2:end-1))/dx^2;
    d2Udy2 = (Ue(3:end,2:end-1)+Ue(1:end-2,2:end-1)-2*Ue(2:end-1,2:end-1))/dy^2;
    d2Vdx2 = (Ve(2:end-1,3:end)+Ve(2:end-1,1:end-2)-2*Ve(2:end-1,2:end-1))/dx^2;
    d2Vdy2 = (Ve(3:end,2:end-1)+Ve(1:end-2,2:end-1)-2*Ve(2:end-1,2:end-1))/dy^2;
    
    ltU = 1/Re_nu/2*(d2Udx2(:)+d2Udy2(:)) - 1/Re_rf/2*U(:);
    ltV = 1/Re_nu/2*(d2Vdx2(:)+d2Vdy2(:)) - 1/Re_rf/2*V(:);
    
    if(save_for_rec && k == omega_save_next)
        omega_save_cnt = omega_save_cnt+1;
        dUdy = (Ue(2:end,:)-Ue(1:end-1,:))/dy;
        dVdx = (Ve(:,2:end)-Ve(:,1:end-1))/dx;
%         omega_t(:,:,omega_save_cnt) = (dVdx-dUdy)/Ts_mks;
         U_t(:,:,omega_save_cnt) = U(:,:);
         V_t(:,:,omega_save_cnt) = V(:,:);
         P_t(:,:,omega_save_cnt) = P(:,:);
        omega_save_next = omega_save_next + omega_save_dn;
    end
    
    % here we compute the Gradient of the pressure
    dPdx = (P(:,2:end) - P(:,1:end-1))/dx;
    dPdy = (P(2:end,:) - P(1:end-1,:))/dy;
    
    % here we compute the rhs of both the U,V equations
    rhsu = U(:) + dt*(ltU(:) - 3/2*nltU(:) + 1/2*nltU_prev(:) - dPdx(:) + F0x(:));
    rhsv = V(:) + dt*(ltV(:) - 3/2*nltV(:) + 1/2*nltV_prev(:) - dPdy(:));
    
    %here we solve the linear equations using matrix inversion
    u(permU) = UTcholU\(LTcholU\rhsu(permU));
    v(permV) = UTcholV\(LTcholV\rhsv(permV));
    
    % here we reshape the vectors back into matrix form
    U = reshape(u,size(U));
    V = reshape(v,size(V));
    
    if(isnan(rms(U(:))) || isnan(rms(V(:))))
        display('Simulation diverged to NAN, try reducing dt');
        break
    end
    %% Here we re-form the extended matrices and impose boundary conditions
    % we read the interior values into the extended matrices
    Ue(2:end-1,2:end-1) = U;
    Ve(2:end-1,2:end-1) = V;
    
    % we now impose the boundary conditions on Ue using boundary points
    Ue(:,1) = 0;
    Ue(:,end) = 0;
    Ue(1,:) = -Ue(2,:);
    Ue(end,:) = -Ue(end-1,:);
    
    % we now impose the boundary conditions on Ve using boundary points
    Ve(:,1) = -Ve(:,2);
    Ve(:,end) = -Ve(:,end-1);
    Ve(1,:) = 0;
    Ve(end,:) = 0;
    
    %% here we compute the pressure correction, the inverse laplacian of div.Vel
    dUdx = (Ue(2:end-1,2:end)-Ue(2:end-1,1:end-1))/dx;
    dVdy = (Ve(2:end,2:end-1)-Ve(1:end-1,2:end-1))/dy;
    DivVel = dUdx+dVdy;
    rhsp = -DivVel(:)/dt;
    phi(permP) = UTcholP\(LTcholP\rhsp(permP));
    % This sets the arbitrary constant in pressure computation to zero
    Phi = reshape(phi,size(P))-mean(phi(:));
    
    % we compute the gradient of the pressure correction
    dPhidx = (Phi(:,2:end) - Phi(:,1:end-1))/dx;
    dPhidy = (Phi(2:end,:) - Phi(1:end-1,:))/dy;
    
    % here we compute the pressure correction
    Phie(2:end-1,2:end-1) = Phi;
    Phie(1,:) = Phie(2,:);
    Phie(end,:) = Phie(end-1,:);
    Phie(:,1) = Phie(:,2);
    Phie(:,end) = Phie(:,end-1);
    
    % The following lines compute the terms of laplacian of Phi
    d2Phidx2 = (Phie(2:end-1,3:end)+Phie(2:end-1,1:end-2)-2*Phie(2:end-1,2:end-1))/dx^2;
    d2Phidy2 = (Phie(3:end,2:end-1)+Phie(1:end-2,2:end-1)-2*Phie(2:end-1,2:end-1))/dy^2;
    
    %% Here we update the final correct velocity and the pressure
    U = U - dt*dPhidx;
    V = V - dt*dPhidy;
    P = P + Phi - dt/2/Re_nu*(d2Phidx2+d2Phidy2) + dt/2/Re_rf*Phi;
    
    %% here we save the non-linear term computed using old velocity
    nltU_prev = nltU;
    nltV_prev = nltV;
    % we just store the value of the norm to be sure of time dependence
    rms_v(k) = rms(V(:));
    rms_u(k) = rms(U(:));
end

% we read the interior values into the extrended matrices
Ue(2:end-1,2:end-1) = U;
Ve(2:end-1,2:end-1) = V;

% we now impose the boundary conditions on Ue
Ue(:,1) = 0;
Ue(:,end) = 0;
Ue(1,:) = -Ue(2,:);
Ue(end,:) = -Ue(end-1,:);

% we now impose the boundary conditions on Ve
Ve(:,1) = -Ve(:,2);
Ve(:,end) = -Ve(:,end-1);
Ve(1,:) = 0;
Ve(end,:) = 0;

dUdy = (Ue(2:end,:)-Ue(1:end-1,:))/dy;
dVdx = (Ve(:,2:end)-Ve(:,1:end-1))/dx;
omega = (dVdx-dUdy)/Ts_mks;

if(save_for_rec && k == omega_save_next)
    omega_save_cnt = omega_save_cnt+1;
%     omega_t(:,:,omega_save_cnt) = (dVdx-dUdy)/Ts_mks;
     U_t(:,:,omega_save_cnt) = Ue;
     V_t(:,:,omega_save_cnt) = Ve;
     P_t(:,:,omega_save_cnt) = P;
end

U0 = (Ue(2:end,:)+Ue(1:end-1,:))/2*Us_mks;
V0 = (Ve(:,2:end)+Ve(:,1:end-1))/2*Us_mks;

dt = save_rate/Ts_mks;      % time-step between saved snapshot (non-dim)
dt_int_mks = dt_mks;        % integration time-step (seconds)
file_name = sprintf('state_%05.2fmA_%.4d.mat',I_mks*1000,tf_mks);
save(file_name,'-v7.3','omega','U','V','P','U0','V0','Us_mks','Ts_mks','Ls_mks','dx','dt','dt_int_mks','tf_mks','ti_mks',...
    'I_mks','rms_u','rms_v','gamma_arr','x','y','nltU_prev','nltV_prev','omega','err_u','err_v','F0x');
if(save_for_rec)
    save(file_name,'U_t','V_t','P_t','omega_save_dn','advpf','Re_nu','Re_rf','-append');
end

if(if_plot)
    max_omega = max(omega(:));
    min_omega = min(omega(:));
    hfig = figure(1);
    
    contourf(X,Y,omega,20,'LineStyle','None');
    daspect([1 1 1]);
    caxis([min_omega,max_omega]);
    colorbar
    hold on
    quiver(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),U0(1:4:end,1:4:end),V0(1:4:end,1:4:end))
    xlabel({'X'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    ylabel({'Y'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    title({sprintf('Simulation  at %05.2fmA',I_mks*1000)},'Fontsize',12,'FontWeight','bold','FontName','Times');
    set(hfig, 'Position', [50 50 750 fix(Ly/Lx*750)]);
    export_fig(sprintf('%05.2fmA_flow.pdf',I_mks*1000),'-transparent');
    close(hfig)
    
    hfig = figure(2);
    plot((1:nt)*dt_mks,rms_u*Us_mks,'b*-');
    hold on
    plot((1:nt)*dt_mks,rms_v*Us_mks,'r*-');
    title({sprintf('Root Mean Square of U (blue) and V (red), I = %05.2f mA',I_mks*1000)},'Fontsize',12,'FontWeight','bold','FontName','Times');
    xlabel({'t (seconds)'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    ylabel({'rms velocity (m/s)'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    set(hfig, 'Position', [50 50 750 fix(Ly/Lx*750)]);
    export_fig(sprintf('%05.2fmA_norm.pdf',I_mks*1000),'-transparent');
    close(hfig)
    
     hfig = figure(3);
    plot((1:nt)*dt_mks,gamma_arr,'k*-');
    title({sprintf('CFL condition for I = %05.2f mA and 1/dt = %03d',I_mks*1000,1/dt_mks)},'Fontsize',12,'FontWeight','bold','FontName','Times');
    xlabel({'t (seconds)'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    ylabel({'(1/dt)(dx)U_{max}'},'Fontsize',12,'FontWeight','bold','FontName','Times');
    set(hfig, 'Position', [50 50 800 300]);
    export_fig(sprintf('%05.2fmA_cfl.pdf',I_mks*1000),'-transparent');
    close(hfig)
    
end
cd ..
clear all
end
