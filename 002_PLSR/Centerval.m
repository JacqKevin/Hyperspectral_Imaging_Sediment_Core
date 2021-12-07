function [ Sc ] = Centerval( Sval, Scal )
% Function that centered the validation set with the calibration set.

% Size of the validation set
[nval, mval]=size(Sval);
[~, mcal]=size(Sval);

% Centering
if mval==mcal
    Sc=Sval-repmat(mean(Scal,'omitnan'),nval,1);
else
    Sc=Sval-repmat(mean(Scal,'omitnan'),nval,mval);
end

end