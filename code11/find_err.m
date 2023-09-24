function find_err(u0,a,N,U,Tfinal,k)

Err = zeros(1,length(N));
for i = 1:length(N)
   EX = find_exact(u0,a,N(i),Tfinal); 
   Err(i) = sum(abs(EX{2} - U{i,2}))*2/N(i);
end

figure(3)
plot(log(2./N),log(Err),'-r','LineWidth',2)
r = polyfit(log(2./N),log(Err),1);
title(sprintf('Slope = %f',r(1)));
xlabel('log(h)')
ylabel('log(Err)')
grid on

