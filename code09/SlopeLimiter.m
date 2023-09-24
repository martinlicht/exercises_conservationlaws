function dU = SlopeLimiter(U,limiter)

dUL = U(2:end-1)-U(1:end-2);
dUR = U(3:end)-U(2:end-1);
switch limiter
    case 'NONE'
        dU = 0*dUL;
    case 'MINMOD'
        dU = minmod([dUL',dUR'])';
    case 'MUSCL'
        dU = minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
    otherwise
        error('Unknown limiter function requested!')
        
end
        

return