function Urec = rec3(U,r,r1)

% Coefficients for interpolation
% row 1 --> r = -1
% row 2 --> r = 0
% row 3 --> r = 1
% row 4 --> r = 2
Crj = [11.0/6.0, -7.0/6.0, 1.0/3.0;  
       1.0/3.0,  5.0/6.0, -1.0/6.0;
       -1.0/6.0, 5.0/6.0, 1.0/3.0;
       1.0/3.0, -7.0/6.0, 11.0/6.0];

Urec = Crj(r+2,1)*U(3-r1:end-2-r1) + Crj(r+2,2)*U(4-r1:end-1-r1) + Crj(r+2,3)*U(5-r1:end-r1);
        
        

return