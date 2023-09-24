clc
clear all
close all

% Discretization	
h  = 0.01;
k  = 0.5*h;
x  = -1:h:1;
nX = numel(x);
t  = 0:k:0.5;
nT = numel(t);

% Initial condition
U0 = 1 + 0.2*sin(2*pi*x) ;

% U1 --> Non-Conservative Scheme
% U2 --> Conservative Scheme
U1 = U0; 
U2 = U0;

% Right-hand side
F = @(x,t) -0.2*2*pi*cos(2*pi*(x-t)) + (1+0.2*sin(2*pi*(x-t))).*(0.2*2*pi*cos(2*pi*(x-t))) ;

for i = 1:nT
    % Create extended solution arrays
    U1_ext = [U1(nT-1),U1,U1(2)] ;
    U2_ext = [U2(nT-1),U2,U2(2)] ;
    
    % Right-hand side
    F_kt = k* F(x,t(i));
    
	% Update solutions
	U1 = U1 - (k/h)*U1_ext(2:end-1).*(U1_ext(2:end-1)-U1_ext(1:end-2)) + F_kt;
    U2 = U2 - 0.5*(k/h)*((U2_ext(2:end-1)).^2-(U2_ext(1:end-2)).^2) + F_kt;
	
    % Exact solution
	U_exact = 1 + 0.2*sin(2*pi*(x-t(i)));
    
	% Visualize the numerical approximation and the correct solution
	plot(x,U1,'-r','LineWidth',2)
    hold all
    plot(x,U2,'-b','LineWidth',2)
    plot(x,U_exact,'-k','LineWidth',2)
	ylim([-0.25 1.75]);xlim([-1 1]);
	set(gca,'XTick',-1:0.5:1);set(gca,'YTick',0:0.5:1.5,'FontSize',15)
	grid on;
	legend('Non-Conservative Form','Conservative Form','True solution','Location','best')
    hold off
	drawnow	
end


