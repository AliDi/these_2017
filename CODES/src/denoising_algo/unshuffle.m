function Gunshuffle = unshuffle(G,varargin)
    % Function called in CanonicalCoherences
    % 
    % Grearrange = UNSHUFFLE(G,ind_shuffle1,ind_shuffle2,...)
    % G results from the shuffle following the ind_shuffle index order : G = A(ind_shuffle1,ind_shuffle2,...)
    % Function UNSHUFFLE restore the orginal order such as Gunshuffle = A.
    %
    %The number of vertor ind_shufflei must be equal to the dimension of G
    %
%
for i=1:(nargin-1)
      [~, ind_unshuffling{i}] = sort(varargin{i});
end
Gunshuffle = G(ind_unshuffling{:});
    

end