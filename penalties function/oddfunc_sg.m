function y = oddfunc_sg(x,gamma,lambda)

x = abs(x) ;
% y = lambda*(1 - (-gamma*x)./exp(x)) ;
y = lambda*(10*(gamma*x+1)-gamma)./exp(10*x);
