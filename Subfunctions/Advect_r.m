%% References:
%[1] Mattia Serra and George Haller, "Efficient Computation of Null-Geodesic with 
%    Applications to Coherent Vortex Detection",  sumbitted, (2016).
%%
% function [~,xxx,yyy,zzz]=Advect_r(phiPrGr,C22mC11Gr,C12Gr,x_glim,y_glim,sVec,Z0,options)  

% Input arguments:
    % C22mC11Gr, C12Gr, phiPrGr : see step 2
    % x1_glim, x2_glim,         : domain limits to stop particles if they
    %                             reach the boundaries 
    % sVec, options             : see step 3
    % Z0                        : Initial conditions 

    % Output arguments:
    %   xxx    : x1-component of trajectories of the ODE (38) in [1]
    %   yyy    : x2-component of trajectories of the ODE (38) in [1]
    %   zzz    : \phi-component of trajectories of the ODE (38) in [1]

%--------------------------------------------------------------------------
% Author: Mattia Serra  serram@ethz.ch
% http://www.zfm.ethz.ch/~serra/
%--------------------------------------------------------------------------
function     [Time,xxx,yyy,zzz] = Advect_r(phiPrGr,C22mC11Gr,C12Gr,x1_glim,x2_glim,sVec,Z0,options) 
    % Shared variables 
    x1m = x1_glim(1);
    x1M = x1_glim(2);
    x2m = x2_glim(1);
    x2M = x2_glim(2);
     Np = numel(Z0)/3;

    % Ode solver
    [Time,XXXpartcl] = ode45(@(t,zVec)r_prime(zVec),sVec,Z0,options); 
        %% Function to evaluate r'(s) (cf. ODE (38) in [1])
        
        function V_intrp = r_prime(zVecs)
        XXxx=zVecs(1:Np);
        YYyy=zVecs(1+Np:2*Np);
        ZZzz=zVecs(1+2*Np:3*Np);
        
        % Freeze particles at the boundaries of the doamin or when the ODE
        % (38) is not defined 
        Bll=1+0*XXxx;
        Bll(XXxx>x1M)=0;Bll(XXxx<x1m)=0;
        Bll(YYyy>x2M)=0;Bll(YYyy<x2m)=0;
        
        % Freeze particles at the boundaries of the domain of definition V (cf eq. (37) of [1]) of
        % the ODE (38) of [1].
            % Domain of existence 
            DoE=2*C12Gr(XXxx,YYyy).*cos(2*ZZzz)+sin(2*ZZzz).*C22mC11Gr(XXxx,YYyy);
            Bll(abs(DoE)<1e-2)=0;
            
        % Evaluate r'(s)
        utemp_scp=cos(ZZzz);
        vtemp_scp=sin(ZZzz);
        wtemp_scp=phiPrGr(XXxx,YYyy,ZZzz);
        Norma=sqrt(utemp_scp.^2+vtemp_scp.^2+wtemp_scp.^2);
        V_intrp=[utemp_scp./Norma.*Bll;vtemp_scp./Norma.*Bll;wtemp_scp./Norma.*Bll];
        end
    xxx=XXXpartcl(:,1:Np);
    yyy=XXXpartcl(:,1+Np:2*Np);
    zzz=XXXpartcl(:,1+2*Np:3*Np);

end









