% Solution09
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The one dimensional linearized acoustic 
% equations are solved using the MUSLC approach, with periodic/open
% boundary conditions, and initial data as given in the exercise.

clc
clear all
close all

% Initial data set
    data = 3;
switch data
    case 1
        pIC =@(x) sin(6*pi*x);
        vIC =@(x) 0*x;
        bc  = 'Periodic';
    case 2
        pIC =@(x) -1.5*(x<0)+ 1.5*(x>=0);
        vIC =@(x) 0*x;
        bc  = 'Open';
    case 3
        pIC =@(x) 1*(x<-0.5) + sin(2*pi*x).*(x>=-0.5).*(x<0.1) + 0*(x>=0.1);
        vIC =@(x) 0*x;
        bc  = 'Open';
end


% Physical constants [u0,k0,p0]
u0 = 1/2;
K0 = 1;
p0 = 1;
c0 = sqrt(K0/p0);   

% Find various matrices
A      = [u0,K0;1/p0,u0];
S      = [-p0*c0, p0*c0; 1,1];
Sinv   = [-1, p0*c0; 1, p0*c0]/(2*p0*c0);
Lambda = [u0-c0,0;0,u0+c0];
absA   = S*abs(Lambda)*Sinv;

% Define Discretization and time parameters
h  = 0.02;
xf = -1:h:1;
xc = (-1+0.5*h):h:(1-0.5*h);
N  = length(xc);
Tfinal = 0.4;
CFL    = 0.5;


% Averaging initial conditions
% Cell-center values sufficient for second-order schemes (can you see why?)
Uavg = [pIC(xc);vIC(xc)];

% U1 --> solution with no slope limiting
% U2 --> solution with minmod limiting
% U3 --> solution with MUSCL limiting
U1 = Uavg;
U2 = Uavg;
U3 = Uavg;

time = 0; iter = 0;
plot_every = 10;

% Solve 
while time < Tfinal
    
    k = CFL*h/(abs(u0) + c0);
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    
    % Obtain RHS for time update
    RHS1 = evalRHS(U1,A,absA,k,h,N,bc,'NONE');
    RHS2 = evalRHS(U2,A,absA,k,h,N,bc,'MINMOD');
    RHS3 = evalRHS(U3,A,absA,k,h,N,bc,'MUSCL');
	
    % Update solution to k
    U1 = U1 + RHS1;
    U2 = U2 + RHS2;
    U3 = U3 + RHS3;
    
    time = time + k;
    iter = iter + 1;
    
    if(mod(iter,plot_every) == 0 || time == Tfinal)
        % Finding exact solution at the current time
        Uexact = find_exact(pIC,vIC,S,Sinv,Lambda,xc,time);
        
        % Visualize the solution
        figure(1)
        subplot(2,1,1)
        plot(xc,Uexact(1,:),'-k','LineWidth',2);
        hold all
        plot(xc,U1(1,:),'.r','LineWidth',2);
        plot(xc,U2(1,:),'-.b','LineWidth',2);
        plot(xc,U3(1,:),'--m','LineWidth',2);
        xlim([-1 1]);
        legend('Exact','NONE','MINMOD','MUSCL','Location','Best')
        grid on;
        title(['Pressure at time = ',num2str(time)])
        hold off
        
        subplot(2,1,2)
        plot(xc,Uexact(2,:),'-k','LineWidth',2);
        hold all
        plot(xc,U1(2,:),'.r','LineWidth',2);
        plot(xc,U2(2,:),'-.b','LineWidth',2);
        plot(xc,U3(2,:),'--m','LineWidth',2);
        xlim([-1 1]);
        legend('Exact','NONE','MINMOD','MUSCL','Location','Best')
        grid on;
        title(['Velocity at time = ',num2str(time)])
        hold off
        
    end
 
end


