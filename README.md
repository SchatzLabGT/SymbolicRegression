# PDE_Discovery_Weak_Formulation
MATLAB codes for papers: 
## Robust and Optimal Sparse Regression for Nonlinear PDE Models
References and authorship are attributed both in the MATLAB codes (provided in the Chaos2019 directory), as well in the paper for which this directory was created. Summary of included files below:

KS_Integrate.m: Integrator for Kuramoto-Sivashinsky equation by Trefethen.

ParEst_WF_KS.m: Implementation of the weak formulation method for the Kuramoto-Sivashinsky equation.

find_coeffs.m: A general solver for sparse regression from the linear system ΘX = 0.

weight_full.m: Helper method for generating weight functions and their derivatives.

weight_poly.m: Helper method for generating weight functions' components in 1D. 

## Using Noisy or Incomplete Data to Discover Models of Spatiotemporal Dynamics

In Time_Integrators all MATLAB codes used to generate data sets can be found, and the data set for the Kuramoto-Sivashinsky equation is given directly (in addition to its integrator). References and authorship are attributed both in the MATLAB codes, as well in the paper for which this directory was created. Brief instructions for use are given below:

KS_Integrator_Trefethen.m : Running this directly should generate the data set in the same directory; some modification to the variable "nplt" on line 31 is necessary for different time resolutions of the stored trajectory.

RD_Integrator_Rudy.m : The length of the trajectory and the desired temporal resolution can be changed in the variable "t" on line 24. Otherwise, running directly should yield the same dataset used in the paper.

initialise_sim_params_fdso_np.m : Run this first to create the necessary files required to run the main integrator. There are two grid resolutions available, 20 grid points per unit length and 40 grid points per unit length, and the forcing profiles for these associated grid resolutions are stored in the files "forcing_profile_20.mat" and "forcing_profile_40.mat", respectively. This code is where spatial grid resolution, integration time-step, as well as physical parameters are determined.

integrate_2dns_fdso_np.m : Main integrator for the 2D Navier-Stokes equation with Kolmogorov forcing. Running this will create a trajectory in the directory that the initialization code creates. This code is where the length of the trajectory and the temporal sampling rate are determined.

In the other folder, Weak_Formulation_Codes, are the functions used to determine the governing equations for the systems detailed in the paper. The titles indicate for which system the function is intended for. Note that the principle input for these MATLAB codes is a string containing the filename for the .mat file that contains the trajectory and relevant parameters (usually only dx & dt, the spatial spacing and time between snapshots, respectively). 

**Important Note 1**
The trajectories must be saved with the switch '-v7.3', which enables them to be loaded in pieces with the matfile formalism in MATLAB. This should be handled automatically in the integrators.

ParEst_wf_KS.m : System identifier for the Kuramoto-Sivashinsky equation. Comes with symbolic regression functionality, though it is not used in the paper.

ParEst_WF_Kolmo.m : System identifier for the 2D Navier-Stokes data.

ParEst_WF_RD.m : System identifier for the lambda-omega reaction diffusion system. Comes with symbolic regression functionality.
