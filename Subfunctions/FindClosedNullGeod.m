%% References:
%[1] Mattia Serra and George Haller, "Efficient Computation of Null-Geodesic with 
%    Applications to Coherent Vortex Detection",  sumbitted, (2016).
%%
% function [x1Psol,x2Psol]=FindClosedNullGeod(C22mC11Gr,C12Gr,phiPrGr,x1_g,x2_g,lamV,sVec,options);   

% Input arguments:
    % C22mC11Gr, C12Gr, phiPrGr : see step 2
    % x1_g, x2_g, lamV          : see steps 0-1
    % sVec, options             : see step 3
    % NCores                    : Number of cores for Parallel computing 

% Output arguments:
    %   x1Psol    : x1-component of closed null-geodesics  
    %   x2Psol    : x2-component of closed null-geodesics  

%--------------------------------------------------------------------------
% Author: Mattia Serra  serram@ethz.ch
% http://www.zfm.ethz.ch/~serra/
%--------------------------------------------------------------------------

function [x1Psol,x2Psol]=FindClosedNullGeod(C22mC11Gr,C12Gr,phiPrGr,x0lam,y0lam,phi0lam,x1_g,x2_g,lamV,sVec,options,NCores);   

    % Initialize the variables containing the periodic solutions of the
    % initial value problem 
    x1Psol = cell(1,length(lamV));
    x2Psol = x1Psol;
    phiPsol = x1Psol;

    % Define the limits of the (x-y) domain to stop particles at the
    % boundary 
    x1_glim = [min(x1_g),max(x1_g)];
    x2_glim = [min(x2_g),max(x2_g)];

tic
    % Compute closed orbits of the Initial Value Problem (cf. eqs. (38-39) of [1])
    for kklam = 1:1:length(lamV)
        kklam

        % Extract the r0_lam for the current value of \lambda
        x0 = x0lam{kklam};
        y0 = y0lam{kklam};
        phi0 = phi0lam{kklam};
        lam=lamV(kklam);

        %% Opening MATLAB Pool %%
        Np=size(x0,1);
        cpu_num = min(NCores,Np);
        id = ceil( linspace(0,Np,cpu_num+1) );

        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)                                          % if parpool is not open
            parpool('local',cpu_num)
        elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
            delete(gcp)
            parpool('local',cpu_num)
        end
        %% Integrate the ODE (38) in [1]
        tic
        spmd
            Range = id(labindex)+1:id(labindex+1);
             Z0loc = [x0(Range);y0(Range);phi0(Range)];  
             [~,xxfTot,yyfTot,zzfTot] = Advect_r(phiPrGr,C22mC11Gr,C12Gr,x1_glim,x2_glim,sVec,Z0loc,options);     
        end
        toc
        
        % Put the trajectories of  ODE (38) in matrix form 
        X_Vf = cat(2,xxfTot{:});
        Y_Vf = cat(2,yyfTot{:});
        Z_Vf = cat(2,zzfTot{:});

        %Warning in case of dimensionality mismatch
        if size(X_Vf,2)~=length(x0)
            disp('smpd dim. mismatch')
        end
        
        % Find periodic solutions 
        [X1lco,X2lco,philco] = PeriodicSolutions(X_Vf,Y_Vf,Z_Vf);

        % Final curves in a cell 
        x1Psol{kklam} = X1lco;
        x2Psol{kklam} = X2lco;
        phiPsol{kklam} = philco;

    end
end

