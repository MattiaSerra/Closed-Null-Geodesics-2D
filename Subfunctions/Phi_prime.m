%% References:
%[1] Mattia Serra and George Haller, "Efficient Computation of Null-Geodesic with 
%    Applications to Coherent Vortex Detection",  sumbitted, (2016).
%%
% [phiPrGr,C22mC11Gr,C12Gr]=Phi_prime(C11,C11x1,C11x2,C12,C12x1,C12x2,C22,C22x1,C22x2,x1_g,x2_g)

% Input arguments:
    % Cij   : ij entries of the C strain tensor 
    % Cijx1 : x1-derivatice of the Cij entry 
    % Cijx2 : x2-derivatice of the Cij entry 
    % x1_g  : x1 component of the spatial grid 
    % x2_g  : x2 component of the spatial grid 

% Output argument:
%   C22mC11Gr: gridded interpolant object for C22(x,y)-C11(x,y)
%   C12Gr    : gridded interpolant object for C12(x,y)
%   phiPrGr  : gridded interpolant object for \phi'(x,\phi)


%--------------------------------------------------------------------------
% Author: Mattia Serra  serram@ethz.ch
% http://www.zfm.ethz.ch/~serra/
%--------------------------------------------------------------------------

function [phiPrGr,C22mC11Gr,C12Gr]=Phi_prime(C11,C11x1,C11x2,C12,C12x1,C12x2,C22,C22x1,C22x2,x1_g,x2_g)

%Define the \phi component of the spatial grid 
phi_g = linspace(0,2*pi,180);

% Constract 3D arrays to build the gridded interpolant for \phi'
C113df = repmat(C11,1,1,numel(phi_g));
C123df = repmat(C12,1,1,numel(phi_g));
C223df = repmat(C22,1,1,numel(phi_g));
C11x3d = repmat(C11x1,1,1,numel(phi_g));
C11y3d = repmat(C11x2,1,1,numel(phi_g));
C12x3d = repmat(C12x1,1,1,numel(phi_g));
C12y3d = repmat(C12x2,1,1,numel(phi_g));
C22x3d = repmat(C22x1,1,1,numel(phi_g));
C22y3d = repmat(C22x2,1,1,numel(phi_g));
[~,~,Z3d]=meshgrid(x1_g,x2_g,phi_g);

% phi' (cf. eq. (38) of [1])
phiPr = -((C11x3d.*cos(Z3d)+C11y3d.*sin(Z3d)).*cos(Z3d).^2+(C12x3d.*cos(Z3d)+C12y3d.*sin(Z3d)).*sin(2*Z3d)+(C22x3d.*cos(Z3d)+C22y3d.*sin(Z3d)).*sin(Z3d).^2)./(sin(2*Z3d).*(C223df-C113df)+2*cos(2*Z3d).*C123df);
phiPrGr = griddedInterpolant ({x1_g,x2_g,phi_g},permute(phiPr,[2 1 3]),'linear');


% Compute the (x)-dependent functions needed to define the domain of
% existence V (cf. eq. (37) of [1]) of the reduced 3D null-geodesic flow (cf. eq. (38) of [1]).
C22mC11Gr = griddedInterpolant ({x1_g,x2_g},permute(C22-C11,[2 1]),'linear');
C12Gr = griddedInterpolant ({x1_g,x2_g},permute(C12,[2 1]),'linear');  
end

