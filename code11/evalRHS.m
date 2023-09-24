% Function evaluates the RHS using a 3-pt stencil to obtain cell-interface
% values

function RHS = evalRHS(U,a,N,bc,crj,k)

% UL and UR are left and right reconstructed values for the cells
% i = 0,...,N+1
% Do not confuse them with the left and right values at interfaces

U = apply_bc(U,bc,k);

% Obtain reconstructed states
UL = zeros(1,N+2);
UR = zeros(1,N+2);

for i = 1:N+2
    [UL(i),UR(i)] = eno_recon(U(i:i+2*k-2),crj,k);
end

% Evaluate flux
flux = 0.5*a*(UL(2:end) + UR(1:end-1)) - 0.5*abs(a)*(UL(2:end)-UR(1:end-1));
% Update solution to k
RHS =  - (flux(2:end) - flux(1:end-1));


return
