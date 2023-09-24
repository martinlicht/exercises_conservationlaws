% Solution12
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The one-dimensional linear advection equation is
% solved using a WENO reconstruction of the 
% solution at the cell interfaces (using cell-averages), 
% followed by an SSP-RK3 integration in time

clc
clear all
close all

% Initial data set
data = 2;
switch data
    case 1
        u0 =@(x) sin(pi*x);
        bc  = 'Periodic';
        Tfinal = 1;
    case 2
        u0 =@(x) 1*(abs(x)<0.5)+ -1*(abs(x)>=0.5);
        bc  = 'Periodic';
        Tfinal = 0.5;
end

% Set time parameters
CFL    = 0.5;
a      = 1;

% Set number of cells 
N = 150:50:250;

% Set order of accuracy 2k-1
k = 3;

% Solve
SOLN = solver(u0,a,bc,N,Tfinal,CFL,k);

% Plot error
if(length(N) > 1)
    find_err(u0,a,N,SOLN,Tfinal,k);
end





