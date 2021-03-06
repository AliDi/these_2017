function [Spp_out , d1, cvx_slvitr ] = CSMReconstruction( Spp)
%calculates diagonal elements d1 to reduce noise on the diagonal of Spp, kipping Spp hermitian nonnegative semi-definite.
%M is the number of microphones (Spp : (MxM))
%d1 : (M x 1)
%from Hald 2016

M=size(Spp,1);
CSM=double(Spp);
M=double(M);
cvx_precision('high')
cvx_begin
	variable d1(M);
	CSM + diag(d1,0) == hermitian_semidefinite(M);
	minimize(sum(d1));
	%minimize(norm(diag(CSM + diag(d1,0)),2));
cvx_end

Spp_out = CSM+diag(d1);

cvx_status;
cvx_slvtol;
cvx_optval;
cvx_slvitr;

%equivalent to :
%cvx_begin
%	variable d2(M)
%	CSM - diag(d2) == hermitian_semidefinite(M);
%	maximize( sum(d2) )
%cvx_end
%(CSM-diag(d2))==(CSM+diag(d1))

end
