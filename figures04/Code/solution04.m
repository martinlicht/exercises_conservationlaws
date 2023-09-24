% Solution03: Problem 5
% This script was written for EPFL MATH459, Numerical Methods for
% Conservation Laws. The scalar advection equation du/dx + du/dt = 0
% with riemann initial data: u = 1 if x<0,  0 if x>0.
% The problem is solved using the following schemes:
%       1. Upwind
%       2. Lax-Friedrichs
%       3. Lax-Wendroff
%       4. Beam-Warming
% The solution is visualized and the accuracy of the scheme tested.

clc
clear all
close all

% Scheme Options: UW --> Upwind
%                 LF --> Lax-Friedrichs
%                 LW --> Lax-Wendroff
%                 BW --> Beam-Warming
% Task Options  : Solve --> To simply run the scheme with h=0.0025 and k/h = 0.5
%                 Acc   --> Test accuracy of the scheme by generating Log-Log
%                           plots
Scheme 	= 'BW';
Task 	= 'Acc';

% Advection speed and final time
a  = 1;
Tf = 0.5;

% Initial condition
U0 =@(x,t) 1*(x-a*t < 0);

% Set the function H for the scheme
%   u^{n+1} = H(u^n)
% as well as the number of ghost cells Ng needed at each end of the mesh (for
% boundary conditions) and the scheme name for saving solutions.
% Here lam = a*k/h and U is the extended solution vector
if(strcmp(Scheme,'UW'))
    Ng    = 1;
    Hfunc =@(U,lam) (1 - lam)*U(2:end-1) + lam*U(1:end-2);
    name  = 'Upwind';
elseif(strcmp(Scheme,'LF'))
    Ng    = 1;
    Hfunc =@(U,lam) (1 - lam)/2*U(3:end) + (1 + lam)/2*U(1:end-2);
    name  = 'Lax-Friedrichs';
elseif(strcmp(Scheme,'LW'))
    Ng    = 1;
    Hfunc =@(U,lam) (1 - lam^2)*U(2:end-1) + (lam^2-lam)/2*U(3:end)...
                     + (lam^2+lam)/2*U(1:end-2);
    name  = 'Lax-Wendroff';
elseif(strcmp(Scheme,'BW'))
    Ng    = 2;
    Hfunc =@(U,lam) (1 - 3*lam/2 + lam^2/2)*U(3:end-2) ...
                    + (2*lam - lam^2)*U(2:end-3)...
                    + (-lam + lam^2)/2*U(1:end-4);
    name  = 'Beam-Warming';
else
    error('Unknown Scheme selected!!');
end


% Performing Tasks
% Run scheme
if strcmp(Task,'Solve')
    fprintf('Solving problem with %s scheme\n',name)
    
    % Discretization
    h  = 0.0025;
    k  = 0.5*h;
    x  = -1:h:1;
    t  = 0:k:Tf;
    N  = numel(x);
    Nt = numel(t);
    
    % Set initial condition
    U = U0(x,0);
    
    % Solve using the numerical scheme. We use open boundary conditions
    lam = a*k/h;
    plot_every = 10; % Plot after every 10 time instances
    
    for i = 2:Nt
        Uext   = [ones(1,Ng)*U(1),U,ones(1,Ng)*U(N)];
        U      = Hfunc(Uext,lam);
        
        if(mod(i,plot_every)==0 || i==Nt)
            Uexact = U0(x,t(i));
            figure(1)
            plot(x,U,'-r','LineWidth',2)
            hold all
            plot(x,Uexact,'-k','LineWidth',2)
            ylim([-0.5,1.5])
            grid on;
            legend('Numerical','Exact','Location','best')
            title([name,', time = ',num2str(t(i))]);
            hold off
            set(gca,'XTick',-0.5:0.5:1,'FontSize',20)
        end
        
    end
    
    % Save plot
    FIG = figure(1);
    fname = sprintf('%s_Profile.pdf',name);
    set(FIG,'Units','Inches');
    pos = get(FIG,'Position');
    set(FIG,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(FIG,fname,'-dpdf','-r0')
    
% Accuracy test
elseif strcmp(Task,'Acc')
    fprintf('Finding accuracy of %s scheme\n',name)
    H = 0.05*2.^(0:-1:-8);
    E = zeros(numel(H),1);
    n = zeros(numel(H),1);
    % Find Error
    for i = 1:numel(H)
        
        % Discretization
        h  = H(i);
        fprintf('... Solving for h = %f\n',h)
        k  = 0.5*h;
        x  = -1:h:1;
        t  = 0:k:Tf;
        N  = numel(x);
        Nt = numel(t);
        
        % Set initial condition
        U = U0(x,0);
        
        % Exact solution at final time
        Uexact = U0(x,Tf);
        
        % Solve using the numerical scheme. We use open boundary conditions
        lam = a*k/h;
        
        for j = 2:Nt
            Uext   = [ones(1,Ng)*U(1),U,ones(1,Ng)*U(N)];
            U      = Hfunc(Uext,lam);
        end
        
        % Measure error in 1 norm
        E(i) = sum(abs(U-Uexact))*h;
        n(i) = N;
    end
    % Make loglog plot
    p = polyfit(log(n),log(E),1);
    FIG = figure(1);
    loglog(n,E,'-ok');grid on;
    xlabel('Resolution','FontSize',20);
    ylabel('Error','FontSize',20);
    title([name,'. Slope = ',num2str(p(1))],'FontSize',20);
    set(gca,'FontSize',20)
    set(FIG,'Units','Inches');
    pos = get(FIG,'Position');
    set(FIG,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fname = sprintf('%s_Accuracy.pdf',name);
    print(FIG,fname,'-dpdf','-r0')
end




