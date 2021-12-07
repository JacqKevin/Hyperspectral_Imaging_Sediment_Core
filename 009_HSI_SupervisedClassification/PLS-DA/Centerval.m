function [ Sc ] = Centerval( Sval, Scal )
% Function that centered the validation set on the calibration set.

% Size of the validation set
[nval, ~]=size(Sval);

% Centering
if size(Scal,1)>1
    Sc=Sval-repmat(mean(Scal),nval,1);
else
    Sc=Sval-repmat(Scal,nval,1);
end

end