function find_err(u0,a,N,U,Tfinal)

Err = zeros(1,length(N));
for i = 1:length(N)
    EX = find_exact(u0,a,N(i),Tfinal);
    Err(i) = sum(abs(EX{2} - U{i,2}))*2/N(i);
end

p = polyfit(log(N),log(Err),1);
figure(3)
loglog(N,Err,'-ok');grid on;
xlabel('N','FontSize',20);
ylabel('Error','FontSize',20);
title(['Slope = ',num2str(p(1))],'FontSize',20);
set(gca,'FontSize',20)

