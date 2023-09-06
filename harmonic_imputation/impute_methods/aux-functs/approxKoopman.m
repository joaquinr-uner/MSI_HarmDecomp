function [Xi,mu,phi_end] = approxKoopman(X,Y,sigma2)
% APPROXKOOPMAN Evaluate Koopman modes, eigen values, and eigen function (from paper of Hua:
% High-dimensional time series prediction using kernel-based Koopman mode regression)

% Usage:	[Xi,mu,phi_end] = approxKoopman(X,Y,sigma2)
%
% Input:
%   X: input dataset
%   Y: output dataset
%   sigma2: shape parameter
%  
% Output:
%   Xi: Koopman modes
%   Mu: Koopman eigen values
%   phi_end: Koopman eigen functions

M = size(X,1) ;

Uxy = [X; Y(end,:)] ; % Eq. (11)

tmp = pdist(Uxy) ;
Uga = exp( -(1/sigma2) * squareform(tmp) ) ; %
Uga = Uga(:,1:M);

Ghat = Uga( 1:M, : ) ;
Ahat = Uga( 2:(M+1), : ) ;

[Q, Sigma2] = eig(Ghat) ;
sigmaPINV = 1./sqrt(diag(Sigma2)) ;

Mtmp = bsxfun(@times, Q, sigmaPINV.') ; % <=> Q*pinv(Sigma)
KoopMat = Mtmp.' * Ahat * Mtmp ; % Eq. (18) 
[Vhat, Mu] = eig(KoopMat) ; % eigen values
mu = diag(Mu);

Phixy = Uga * Mtmp * Vhat ; % Eq. (24)

Xi = pinv(Phixy) * Uxy ; % After Eq. (25)

phi_end = Phixy(end,:);
