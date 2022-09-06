function X = prox_nonconvex_plus(Y,lambda,fun,gamma)

% min_X lambda*sum_ g(sigmai(X)) * sigma(X) + 0.5*||X-Y||_F^2
% g(t) = 1 - (?1|x|+1)/e^(?2|x|)(fun)

J = 50;
fun = str2func([fun '_sg']);
% fun = str2func(fun);

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);

[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
Z = zeros(length(S),1);
Z1 = zeros(length(S),1);
for iter = 1:J
    punish = fun(Z,gamma,lambda);
    Z1 =  max(S - punish,0) ;% S(j) - lambda * (1+gamma)*gamma/((abs(S(j))+gamma)^2);
    err = sum((Z1-Z).^2);
    Z = Z1;
    if err < 1e-6
        break
    end
end
X(:,:,1) = U*diag(Z)*V';

for i = 2:round(n3/2)
    [U, S, V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    Z = zeros(length(S),1);
    Z1 = zeros(length(S),1);
    for iter = 1:J
        punish = fun(Z,gamma,lambda);
        Z1 =  max(S - punish,0) ;% S(j) - lambda * (1+gamma)*gamma/((abs(S(j))+gamma)^2);
        err = sum((Z1-Z).^2);
        Z = Z1;
        if err < 1e-6
            break
        end
    end
    X(:,:,i) = U*diag(Z)*V';
    X(:,:,n3+2-i) = conj(X(:,:,i));
end

if mod(n3,2) == 0
    i = round(n3/2)+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    Z = zeros(length(S),1);
    Z1 = zeros(length(S),1);
    for iter = 1:J
        punish = fun(Z,gamma,lambda);
        Z1 =  max(S - punish,0) ;% S(j) - lambda * (1+gamma)*gamma/((abs(S(j))+gamma)^2);
        err = sum((Z1-Z).^2);
        Z = Z1;
        if err < 1e-6
            break
        end
    end
    X(:,:,i) = U*diag(Z)*V';
    
end
X = ifft(X,[],3);





% maxiter = 100;
% mu = 1.1;
% tol = 1e-5;
% hfun_sg = str2func([fun '_sg']) ;
% Y = fft(Y,[],3);
% [n1,n2,n3] = size(Y);
% M = opRestriction(prod([n1,n2]),ones(n1*n2,1));
%
% for i = 1:n3
%     y = M(Y(:,:,i),1);
%     x = zeros(n1*n2,1);
%     lambdaInit = 0.9*max(abs(M(y,2)));
%     lambda = lambdaInit;
%     f_cur = norm(y-M(x,1)) + lambda*norm(x,1);
%     while lambda > lambdaInit*tol%ÍâÑ­»·
%         for j = 1:maxiter
%             f_pre = f_cur;
%             x = x + (1/mu)*M(y-M(x,1),2);
%             [U,S,V] = svd(x,'econ');
%             sigma = diag(S);
%
%             w = hfun_sg(sigma,gamma,lambda);
%             sigma = sigma - w/mu;
%             svp = length(find(sigma>0)) ;
%
%             sigma = sigma(1:svp) ;
%             X = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
%             x = X(:);
%             f_cur = norm(y-M(x,1)) + lambda*norm(x,1);
%
%             if norm(f_cur-f_pre)/norm(f_cur + f_pre)<tol
%                 break;
%             end
%         end
%         if norm(y-M(x,1))<err
%             break;
%         end
%     end
%     lambda = 0.9*lambda;
% end