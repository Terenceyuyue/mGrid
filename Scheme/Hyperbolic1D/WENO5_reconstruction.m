function q = WENO5_reconstruction(U,ind)
% 5th order WENO reconstruction for 1-D scalar or vectorial problems with
% periodic condition.
% 
% U: nv-by-Nx, where nv is the length of the variables of vectorial
%    problems and Nx is the number of cells. That is, each row of U
%    corresponds to a variable.
%
% positive flux at x_{i+1/2}: ind = [0 +1]
% negative flux at x_{i-1/2}: ind = [0, -1]
%
% Copyright (C)  Terence Yu
%
% Reference: https://github.com/wme7/WENO


% periodic permutation
vmm = circshift(U,ind*2);
vm  = circshift(U,ind);
v = U;
vp  = circshift(U,-ind);
vpp = circshift(U,-ind*2);

% reconstruction polynomials
p0 = (2*vmm-7*vm+11*v)/6;
p1 = (-vm+5*v+2*vp)/6;
p2 = (2*v+5*vp-vpp)/6;

% Smooth parameters
IS0 = 13/12*(vmm-2*vm+v).^2 + 1/4*(vmm-4*vm+3*v).^2; 
IS1 = 13/12*(vm-2*v+vp).^2 + 1/4*(vm-vp).^2;
IS2 = 13/12*(v -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% constants
c0 = 1/10; c1 = 3/5; c2 = 3/10; epsilon = 1e-6;

% alpha 
alpha0 = c0./(epsilon+IS0).^2;
alpha1 = c1./(epsilon+IS1).^2;
alpha2 = c2./(epsilon+IS2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% Non-linear weights
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% WENO5 polynomial
q = w0.*p0 + w1.*p1 + w2.*p2;