function [S_snv] = snv(S)
% Pretreated your spectra with SNV (Standard Normal Variate).

% input:
% x data to pretreat
%
% output:
% x_snv pretreated data

[~,n]=size(S); % Size of x
xm=mean(S,2); % Mean of x
xc=S-repmat(xm,1,n); % x center on their mean
S_snv=xc./repmat(std(xc,[],2),1,n); % xc standardised by the standard deviation
end