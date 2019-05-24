function b = squeeze2(a)
%SQUEEZE Remove ALL singleton dimensions.
%   B = SQUEEZE(A) returns an array B with the same elements as
%   A but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(A,dim)==1.
% A row vector becomes a column vector
%
%   For example,
%       squeeze(rand(2,1,3))
%   is 2-by-3.
%
%   Variante par Alice
%   See also SHIFTDIM.

%   Copyright 1984-2010 The MathWorks, Inc.

if nargin==0
    error(message('MATLAB:squeeze:NotEnoughInputs'));
end
siz = size(a);
siz(siz==1) = []; % Remove singleton dimensions.
siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
b = reshape(a,siz);

end
