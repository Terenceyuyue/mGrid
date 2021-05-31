function [q2,q1] = WENO5_reconstruction(FU2,FU1)
% 5th order WENO reconstruction for 1-D scalar or vectorial problems with
% periodic condition.
%
% FU2, q2: positive flux at x_{j+1/2}
% FU1, q1: negative flux at x_{j+1/2}
% 
% FU: nv-by-Nx, where nv is the length of the variables of vectorial
%    problems and Nx is the number of cells. That is, each row of FU
%    corresponds to a variable.
%
% Copyright (C)  Terence Yu
%
% Reference: https://github.com/wme7/WENO

%% Positive flux at x_{j+1/2}
% periodic permutation
vmm = circshift(FU2,[0,2]);
vm  = circshift(FU2,[0,1]);
v = FU2;
vp  = circshift(FU2,[0,-1]);
vpp = circshift(FU2,[0,-2]);

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
q2 = w0.*p0 + w1.*p1 + w2.*p2;

%% Negative flux at x_{j+1/2}
% periodic permutation
um = circshift(FU1,[0,1]);
u  = FU1;
up = circshift(FU1,[0,-1]);
upp = circshift(FU1,[0,-2]);
uppp = circshift(FU1,[0,-3]);

% reconstruction polynomials
p0 = (-um+5*u+2*up)/6;
p1 = (2*u+5*up-upp)/6;
p2 = (11*up-7*upp+2*uppp)/6;

% Smooth parameters
IS0 = 13/12*(up-2*u+um).^2 + 1/4*(3*up-4*u+um).^2; 
IS1 = 13/12*(upp-2*up+u).^2 + 1/4*(upp-u).^2;
IS2 = 13/12*(uppp-2*upp+up).^2 + 1/4*(uppp-4*upp+3*up).^2;

% constants
c0 = 3/10; c1 = 3/5; c2 = 1/10; epsilon = 1e-6;

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
q1 = w0.*p0 + w1.*p1 + w2.*p2;