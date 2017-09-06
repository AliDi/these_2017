function S=InterspectraMatrix(MSignals)
%S=CSM(MSignals) : calculates the matrix of auto/inter spectra of Msignals
%
%MSignals : M x 1 x Nsamples matrix
%M : number of signals (ie number of mics or sources)
%Nsamples : numbers of samples
%
%S : M x M x Nsamples matrix
%

    [M , ~, Nsamples] = size(MSignals);
    for f=1:Nsamples
        S(:,:,f) = MSignals(:,1,f) * MSignals(:,1,f)';
    end
    S=S/Nsamples;
end