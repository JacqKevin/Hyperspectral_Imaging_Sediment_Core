function [S_msc]=msc(S,Ref)
% Pretreated your spectra with MSC (Multiplicative Scatter Correction).

% input:
% x data to pretreat
%
% output:
% x_msc pretreated data

% If the user do not choose a reference vector, that will be the mean of x
if nargin<2
    Ref=mean(S);
end

[m, n]=size(S); % Size of x
mz=[ones(1,n); Ref]';
wmz=mz.*ones(n,2);
wz=S.*ones(m,n);
z=wmz'*wmz;
[u,s,v]=svd(z); % Singular value decomposition
sd=diag(s)';
cn=10^12;
ms=sd(1)/sqrt(cn);
cs=max(sd,ms);
zi=u*(diag(cs))*v';
B=(zi\(wmz'*wz'))';
S_msc=(S-(B(:,1)*ones(1,n)))./(B(:,2)*ones(n,1)'); % MSC correction
end