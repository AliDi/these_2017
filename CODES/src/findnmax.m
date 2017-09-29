function [val subscripts]=findnmax(A, nb_max)
%[subscripts val]=findnmax(A, nb_max)
%
%subscripts : (nb_max x 2)
%returns nb_max maximum values of A in val


A_vect = reshape(A,1,numel(A));

for n=1:nb_max	
	[val(n) ind]= max(A_vect);
	[subscripts(n,1) subscripts(n,2)] = ind2sub(size(A), ind);
	A_vect(ind) = A_vect(ind) - val(n);
end
