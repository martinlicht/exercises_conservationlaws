% Function evaluates the RHS using a 3-pt stencil to obtain cell-interface
% values

function RHS = evalRHS(U,a,bc,r)

U = apply_bc(U,bc,3);

% Obtain reconstructed states
% UL and UR are left and right reconstructed values for the cells
% i = 0,...,N+1
% Do not confuse them with the left and right values at interfaces
UR = rec3(U,r,r);
UL = rec3(U,r-1,r);

% Evaluate flux
flux = 0.5*a*(UL(2:end) + UR(1:end-1)) - 0.5*abs(a)*(UL(2:end)-UR(1:end-1));

% Find flux difference
RHS =  - (flux(2:end) - flux(1:end-1));


return
