% This is the main solver that approximates the solution

function SOLN = solver(u0,a,bc,N,Tfinal,CFL,k)

% Initialize polynomial coefficients (of degree k-1)
Crec = zeros(k+1,k);
for r=-1:k-1;
    Crec(r+2,:) = ReconstructWeights(k,r);
end;

% Initialize linear
dw = LinearWeights(k,0);

% Compute smoothness indicator matrix
beta = zeros(k,k,k);
for r=0:k-1
    xld = -1/2 + [-r:1:k-r];
    beta(:,:,r+1) = betarcalc(xld,k);
end

print_iter = 1000;

SOLN = {length(N),2};

for i = 1:length(N)
    
    fprintf('... solving on mesh with %d cells \n',N(i));
    
    % DOMAIN AND TIME PARAMETERS:
    h  = 2/N(i);
    xf = -1:h:1;
    xc = -1+0.5*h:h:1-0.5*h;
    
    U = zeros(1,N(i));
    for j = 1:N(i)
        U(j) = integral(u0,xf(j),xf(j+1),'AbsTol',1e-14)/h;
    end
    
    time = 0;
    iter = 0;
    
    while(time<Tfinal)
        dt = CFL*h/abs(a);
        if(time + dt > Tfinal)
            dt = Tfinal - time;
        end
        
        Uold = U;
        
        % SSP-RK3 stage 1
        RHS = evalRHS(U,a,N(i),bc,k,Crec,dw,beta);
        U   = Uold + (dt/h)*RHS;
        
        % SSP-RK3 stage 2
        RHS = evalRHS(U,a,N(i),bc,k,Crec,dw,beta);
        U   = 3*Uold/4.0 + (U + (dt/h)*RHS)/4.0;
        
         % SSP-RK3 stage 3
        RHS = evalRHS(U,a,N(i),bc,k,Crec,dw,beta);
        U   = Uold/3.0 + 2.0*(U + (dt/h)*RHS)/3.0;
        
        time = time + dt;
        
        if(mod(iter,print_iter)==0 || time == Tfinal)
            
            title_text = sprintf('t = %0.3f, N=%d',time,N(i));
            figure(1)
            plot(xc,U,'-','LineWidth',2)
            xlabel('x','FontSize',20)
            title(title_text,'FontSize',20);
            set(gca,'FontSize',20)
            
            pause(0.1)
        end
        iter = iter + 1;
    end
    
    SOLN{i,1} = xc;
    SOLN{i,2} = U;
    
end

EX = find_exact(u0,a,1000,Tfinal);

for i = 1:length(N)
figure(2)
   plot(SOLN{i,1},SOLN{i,2},'LineWidth',2)
   hold all
end
plot(EX{1},EX{2},'k-','LineWidth',2)

legendcell = cellstr(num2str(N', 'N = %-d'));
strcat(legendcell,'Exact');
legend(legendcell,'Location','best')
xlabel('x')
ylabel('u')


return;
