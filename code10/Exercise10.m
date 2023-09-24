% Solution 10
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The one-dimensional linear advection equation is
% solved using a 3rd order reconstruction of the 
% solution at the cell interfaces (using cell-averages), 
% followed by an SSP-RK3 integration in time

clc
clear all
close all

% Initial data set
data = 1;
switch data
    case 1
        u0 =@(x) sin(pi*x);
        bc  = 'Periodic';
        Tfinal = 0.5;
    case 2
        u0 =@(x) 1*(x<0)+ 0*(x>=0);
        bc  = 'Open';
        Tfinal = 0.5;
end

% Choose stencil to reconstruct the values in cell i
%     DW : uses cells i, i+1, i+2
%     CEN: uses cells i-1, i, i+1
%     UW : uses cells i-2, i-1, i
stencil = 'UW';
switch stencil
    case 'DW'
        r = 0;
    case 'CEN'
        r = 1;
    case 'UW'
        r = 2;
end

% Set time parameters
CFL    = 0.2;
a      = 1;

% Set number of cells 
N = [100,200,400,800];

% Solve
SOLN = solver(u0,a,bc,N,Tfinal,CFL,r);
find_err(u0,a,N,SOLN,Tfinal);




