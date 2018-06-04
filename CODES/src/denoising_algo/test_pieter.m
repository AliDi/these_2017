function [Spp_out , d1, L ] = test_pieter( Spp)
%calculates diagonal elements d1 to reduce noise on the diagonal of Spp, kipping Spp hermitian nonnegative semi-definite.
%M is the number of microphones (Spp : (MxM))
%d1 : (M x 1)
%from Hald 2016

M=size(Spp,1);
CSM=double(Spp);
M=double(M);
cvx_precision('low')
%cvx_solver sedumi
cvx_begin
	%Pieter
	%variable d1(M,M) nonnegative
	%minimize( norm_nuc(CSM-d1)+0.1*norm(d1(:),1) )	
	variable d1(1,M);
	variable L(M,M) complex;
	%L == hermitian_semidefinite(M);
	%d1 == hermitian_semidefinite(M);
	minimize( norm_nuc(L)+0.5*norm(d1(:),1) )
	subject to
		norm(L + diag(d1) - CSM,2) <= 0.0001;
		%L+d1==CSM;
cvx_end


Spp_out = CSM-(d1);

end

%cvx_begin
%	variable d1(M);
%	CSM + diag(d1,0) == hermitian_semidefinite(M);
%	minimize(sum(d1));
%	%minimize(norm(diag(CSM + diag(d1,0)),2));
%cvx_end
