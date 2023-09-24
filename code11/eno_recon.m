function [UL,UR] = eno_recon(U,crj,k)

r  = find_shift(U,k);
UL = sum(crj(r+1,:).*U(k-r:2*k-1-r));
UR = sum(crj(r+2,:).*U(k-r:2*k-1-r));


return