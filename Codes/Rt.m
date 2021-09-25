
%lstar = (mu1*-1)+(mu2*1);
lstar =mu'*xstar;
RT(1) = 0;
for i=1:T
RT(i+1) = RT(i) + ([mu1,mu2]*x(:,:,i)- lstar);
end
plot(RT1)
hold on
plot(RT)
title ('\Delta > 0, T=100000')
%title ('\Delta = 0, x^* on the tangent line')
xlabel('t')
ylabel('Cumulative Regret')
legend('without randomness','with randomness')