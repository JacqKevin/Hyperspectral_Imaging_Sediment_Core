function [ Sc ] = Uncenterval( Sval, Scal )
% Uncentered the prediction with the calibration set.

% Size of the set to uncentered
[nval, mval]=size(Sval);
[~, mcal]=size(Scal);

% Uncentered
if mval==mcal
    Sc=Sval+repmat(mean(Scal),nval,1);
else
    Sc=Sval+repmat(mean(Scal),nval,mval);
end

end