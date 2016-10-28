%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Author: Mattia Serra                  %
                    % Email:  serram@ethz.ch                %
                    % Date:   26/10/2016                    %
                    % Web:  http://www.zfm.ethz.ch/~serra/  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The present code computes closed null-geodesics in two-dimensianal
% dynamical systems using the theoretical results derived in [1].
% Specifically (cf. section 4.3 of [1]), this code computes elliptic LCSs 
% or coherent Lagrangian  vortices [2].
% With minor modifications, it can be used to compute any null-geodesics
% (see. [1]).

% Reference:
%[1] Mattia Serra and George Haller, "Efficient Computation of Null-Geodesic with 
%    Applications to Coherent Vortex Detection",  sumbitted, (2016).
%[2] G. Haller, FJ. Beron-Vera, Coherent Lagrangian vortices: the black
%     holes of turbulence.  J. Fluid Mech. 731 (2013) R4
%% Step0: Load the Tensor field (Right Cauchy Green) entries and their spatial derivatives 
clear all; clc; close all; tic
load('TensorFieldData.mat','x1_g','x2_g','lam2','C11','C12','C22','C11x1','C11x2','C12x1','C12x2','C22x1','C22x2')

% x1_g : x1 component of the spatial grid 
% x2_g : x2 component of the spatial grid 
% lam2 : highest Cauchy-Green (C) eigenvalue (Only needed for the final FTLE plot)
% Cijx1: x1-derivatice of the Cij entry 
% Cijx2: x2-derivatice of the Cij entry 
% Cij  : ij entries of the C strain tensor 

% Note: The velocity field used in this work are derived using the altimetry produced by SSALTO/DUACS and
% distributed by AVISO, with support from CNES (http://www.aviso.oceanobs.com/duacs). 
%% Step1: Compute r_\lambda(0)=(x_0,\phi0) (cf. eq. (39) of [1])   
    %Define the desired set of \lambda values 
    lamV = linspace(0.9,1.1,7);
                                                       
    %Compute r_\lambda(0)=(x_0,\phi0): the \lambda-dependent initial conditions for the reduced 3D null-geodesic flow (eq. (39) of [1]) 
    [x0lam,y0lam,phi0lam]=r0_lam(lamV,C11,C12,C22,x1_g,x2_g);
%% Step2: Compute \phi'
    % Returns a function handle for \phi'(x,y,\phi) (cf. eq.(38) of [1])
    % and the (x)-dependent functions needed to evaluate its domain of existence V. (cf. eq (37) of [1])
    [phiPrGr,C22mC11Gr,C12Gr]=Phi_prime(C11,C11x1,C11x2,C12,C12x1,C12x2,C22,C22x1,C22x2,x1_g,x2_g);    
   
        % Clear unnecessary variables 
        clear C11 C11x1 C11x2 C12 C12x1 C12x2 C22 C22x1 C22x2
%% Step3: Solve the Initial value null-geodesic problem (cf. eqs. (38-39) of [1])
    % Parameterization of r(s)
    sVec = [0:0.001:12];
    % ODE solver options
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);
    % Select the number of cores for Parallel computing 
    NCores=25;                   
    % Compute closed Null-Geodesics  
    [x1Psol,x2Psol] = FindClosedNullGeod(C22mC11Gr,C12Gr,phiPrGr,x0lam,y0lam,phi0lam,x1_g,x2_g,lamV,sVec,options,NCores);   
    
    % Plot all the closed null-geodesics on FTLE
    PlotAllClosedNullGeodesics(x1Psol,x2Psol,x1_g,x2_g,lamV,lam2)
 %% Step4: Find outermost closed null-geodesics 

    % Find Outermost closed null-geodesics 
    [x1LcOutM,x2LcOutM,LamLcOutM] = FindOutermost(x1Psol,x2Psol,lamV,sVec);
   
    % Plot Outermost Closed null-geodesics on FTLE
    PlotOutmost(x1LcOutM,x2LcOutM,LamLcOutM,lamV,x1_g,x2_g,lam2)    