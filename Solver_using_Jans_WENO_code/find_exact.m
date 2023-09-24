% Function evaluates the (averaged) exact solution at a given time

function EXSOLN = find_exact(u0,a,N,time)

h = 2.0/N;
xf = -1:h:1;
xc = -1+0.5*h:h:1-0.5*h;
U = zeros(1,N);
for j = 1:N
    U(j) = integral(u0,xf(j)-a*time,xf(j+1)-a*time,'AbsTol',1e-14)/h;
end
        
EXSOLN{1} = xc;
EXSOLN{2} = U;