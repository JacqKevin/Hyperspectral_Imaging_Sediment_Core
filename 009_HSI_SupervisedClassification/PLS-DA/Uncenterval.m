function [ Sc ] = Uncenterval( Sval, Scal )
% Uncentered the set with the mean reference.

% Size of the set to uncentered
[nval, ~]=size(Sval);

% Uncentered
Sc=Sval+repmat(mean(Scal),nval,1);

end

