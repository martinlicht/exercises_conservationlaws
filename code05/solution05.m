% Solution05: Problem 3, 4 and 5
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The invisvid burgers equation is solved using two 
% different sets of initial condtions using an
% upwind inspired scheme, EO flux and the generalized Lax-Friedrichs scheme. 
% The first method is non-conservative, while the next two are conservative. 

clc
clear all
close all

% Initial left and right states for the Riemann Problem
% u0(x) = ul if x<0.0, ur if x>0.0
ul = 1;  ur = 0; 


% Discretization	
h  = 0.01;
k  = 0.5*h;
x  = -0.5:h:1;
nX = numel(x);
t  = 0:k:0.5;
nT = numel(t);

% Initial condition
U0 = ul*(x<0) + ur*(x>=0);

% U1 --> solution of Upwind inspired scheme
% U2 --> solution with EO flux
% U3 --> solution of Generalized Lax-Friedrichs scheme
U1 = U0; 
U2 = U0;
U3 = U0;

% EO flux (minima of Burgers flux at u=0)
EOflux =@(u,v) max(u,0).^2/2 + min(v,0).^2/2 - 0; 


% Shock speed
s  = (ul+ur)/2; 

for i = 1:nT
    
    % Create extended solution arrays
    U1_ext = [U1(1),U1,U1(end)];
    U2_ext = [U2(1),U2,U2(end)];
    U3_ext = [U3(1),U3,U3(end)];
    
	% Update solutions
	U1 = U1 - (k/h) * U1_ext(2:end-1).* (U1_ext(2:end-1) - U1_ext(1:end-2));
    
    EOfluxeval = EOflux(U2_ext(1:end-1),U2_ext(2:end));
    U2 = U2 - (k/h) * (EOfluxeval(2:end) - EOfluxeval(1:end-1)); 
    
	U3 = ( U3_ext(3:end) + U3_ext(1:end-2))/2 - ...
         (0.5*k/h) * ( U3_ext(3:end).^2 - U3_ext(1:end-2).^2 )/2;
	
    % The true solution, shock speed calculated using Huginiot jump condition
	U_exact = ul*(x<s*t(i)) + ur*(x>=s*t(i));
    
	% Visualize the numerical approximation and the correct solution
	plot(x,U1,'-r','LineWidth',2)
    hold all
    plot(x,U2,'-b','LineWidth',2)
    plot(x,U3,'-m','LineWidth',2)
    plot(x,U_exact,'-k','LineWidth',2)
	ylim([-0.25 1.75]);xlim([-0.5 1]);
	set(gca,'XTick',-0.5:0.5:1);set(gca,'YTick',0:0.5:1.5,'FontSize',15)
	grid on;
	legend('Upwind inspired','EO','Generalized Lax-Friedrichs','True solution','Location','best')
    hold off
	drawnow	
end


