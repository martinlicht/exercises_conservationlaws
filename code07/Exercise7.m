% Solution07
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws.
% The one dimensional linearized acoustic equations are solved with
% periodic/open boundary conditions, and initial data
% as given in the exercise.

clc
clear all
close all

% Initial data set
data = 3;
switch data
    case 1
        pIC =@(x) sin(2*pi*x);
        vIC =@(x) 0*x;
        bc  = 'Periodic';
    case 2
        pIC =@(x) 1*(x>0);
        vIC =@(x) 0*x;
        bc  = 'Open';
    case 3
        pIC =@(x) 1*(x<0) + sin(2*pi*x).*(x>=0);
        vIC =@(x) 0*x;
        bc  = 'Open';
end

% Define Discretization and time parameters
h  = 0.01;
xf = -1:h:1;
xc = (-1+0.5*h):h:(1-0.5*h);
N  = length(xc);
Tfinal = 0.4;
CFL    = 0.5;

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


% Averaging initial conditions
% Cell-center values sufficient for first-order schemes
U = [pIC(xc);vIC(xc)];

time = 0; iter = 0;
plot_every = 10;

% Solve
while time < Tfinal
    
    k = CFL*h/(abs(u0) + c0);
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    
    % Applying boundary conditions to obtain extended vector
    U_ext = apply_bc(U,bc);
    
    Flux = GodunovFlux(A,absA,U_ext(:,1:end-1),U_ext(:,2:end));
    U = U - k/h*(Flux(:,2:end) - Flux(:,1:end-1));
    
    time = time + k;
    iter = iter + 1;
    
    if(mod(iter,plot_every)==0 || time ==Tfinal)
        % Finding exact solution at the current time
        Uexact = find_exact(pIC,vIC,S,Sinv,Lambda,xc,time);
        
        % Visualize the solution
        figure(1)
        subplot(2,1,1)
        plot(xc,U(1,:),'-r','LineWidth',2);
        hold all
        plot(xc,Uexact(1,:),'--k','LineWidth',2);
        ylim([-2 2]);xlim([-1 1]);
        legend('Numerical Pressure','Exact Pressure','Location','Best')
        grid on;
        title(['Time = ',num2str(time)])
        hold off
        
        subplot(2,1,2)
        plot(xc,U(2,:),'-r','LineWidth',2);
        hold all
        plot(xc,Uexact(2,:),'--k','LineWidth',2);
        ylim([-2 2]);xlim([-1 1]);
        legend('Numerical Velocity','Exact Velocity','Location','Best')
        grid on;
        hold off
        
    end
    
end


