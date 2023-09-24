% Function evaluates the RHS corresponding to a second-order MUSCL scheme
% with a mid-point rule for time-integration

function RHS = evalRHS(U,A,absA,k,h,N,bc,limiter)

% Need to extend by 2 ghost cells on either side
U_ext = apply_bc(U,bc,2);

% Obtain limited slope for N+2 cells
dU      = zeros(2,N+2);
dU(1,:) = SlopeLimiter(U_ext(1,:),limiter);
dU(2,:) = SlopeLimiter(U_ext(2,:),limiter);

% Obtain cell solution at k/2 with f'(u) = A
Unph = U_ext(:,2:end-1) - 0.5*k/h*A*dU;

% Obtain interface values at k/2
UL = Unph(:,1:end-1) + 0.5*dU(:,1:end-1);
UR = Unph(:,2:end)   - 0.5*dU(:,2:end);

% Evaluate flux
flux = GodunovFlux(A,absA,UL,UR);

% Update solution to k
RHS =  - k/h*(flux(:,2:end) - flux(:,1:end-1));


return
