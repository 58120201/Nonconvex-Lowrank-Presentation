function [L, S] = trpca_explr_plus(X,lambda,opts)
%TRPCA_EXPLR 此处显示有关此函数的摘要
%   此处显示详细说明 IRNN
%   min_{L,S}  sum_i( g(sigmai(L))*sigmai(L) + lambda*||S||_1
%   s.t.         X = L + S

tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

[m,n,p] = size(X);
% Z = zeros(min(m,n),p);

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

dim = size(X);
L = zeros(dim);
S = L;
Y = L;

for iter = 1 : max_iter
    Lk = L;
    Sk = S;     
    % update L
    L = prox_nonconvex_plus(-S+X-Y/mu,1/mu,opts.fun,opts.gamma);
    % update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);
    % S = prox_l21(-L+X+Y/mu,lambda/mu);
    
    dY = L+S-X;
    % dY = X-L-S;
    chgL = max(abs(Lk(:)-L(:))); 
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            err = norm(dY(:));
            disp(['iter=' num2str(iter) ',mu=' num2str(mu) ',err=' num2str(err)]);
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);   
end

