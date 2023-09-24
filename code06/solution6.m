% Solution06: Problem 3
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The invisvid Burgers equation is solved using the
% Roe flux

clc
clear all
close all

% Initial left and right states for the Riemann Problem
U0 = @(x) -1*(x<1) + 1*(x>1);

FinalTime = 0.25;

% Roe flux
RoeFlux =@(UL,UR) 0.5*UL.^2.*( (UL+UR) >=0) + 0.5*UR.^2.*( (UL+UR) <0);


% Discretization
l = 50:10:100;
%l = 51:10:101;
%l = 50:5:100;

M = length(l);

for i=1:M
    h = 2/l(i);
    k = 0.5*h;
    xf = 0:h:2;
    xc = 0.5*h:h:2-0.5*h;
    N  = length(xc);
    
    % Finding cell-average values of initial condition
    % Since we are looking at first-order schemes, it is enough to
    % approximate the integrals using the mid-point quadrature rule. This
    % essentially means we take the cell-center values as the cell-averages
    U = U0(xc);
    
    time = 0;
    while(time<FinalTime)
        
        if(time+k>FinalTime)
            k = FinalTime-time;
        end
        
        U_ext = [U(1),U,U(end)];
        
        Flux  = RoeFlux(U_ext(1:end-1),U_ext(2:end));
        
        U = U - (k/h)*(Flux(2:end) - Flux(1:end-1));
        
        time = time + k;
    end
    
    figure(1)
    plot(xc,U,'LineWidth',2)
    title(['Solution on a mesh with h = 2/',num2str(l(i))])
    hold all
    pause(0.1)
    
end


