function [ Sc ] = Center( S )
% Function that centered the data.

[n, ~]=size(S);
Sc=S-repmat(mean(S,'omitnan'),n,1);

end

