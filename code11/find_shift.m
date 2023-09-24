% Finds the shift r for choosing the approriate ENO stencil of size k

function r = find_shift(U,k)

eps = 1.0e-12;
N = length(U);
assert(N == 2*k-1)

Dmat = zeros(N,N);
Dmat(1,:) = U;

for i = 2:N
    Dmat(i,1:N-i+1) = Dmat(i-1,2:N-i+2) - Dmat(i-1,1:N-i+1);
end

r = 0;
for i=1:k-1
    if(abs(Dmat(i+1,k-1-r)) < (abs(Dmat(i+1,k-r)) - eps))
        r = r + 1;
    end
end

return