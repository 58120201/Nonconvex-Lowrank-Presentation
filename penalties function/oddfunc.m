function y = oddfunc(x,gamma,lambda)

x = abs(x) ;
% y = lambda*(1 - (-gamma*x)./exp(x)) ;
y = lambda*(1-((gamma*x+1)./exp(10*x)));
