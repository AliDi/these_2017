function [CSM freqs]=CSMtemp(data_temp , Nfft , Fs)
%Estimates CSM from temporal signals
%
%data_temp (Nmic x Nsamples) : temporal signals
%Nfft : number of points for each segmentation for the FFT (50% overlap)
	
	[Nmic] = size(data_temp,1);
	
	for i=1:Nmic
		[data_spec(i,:,:) freqs]=spectrogram(data_temp(i,:), hamming(Nfft),floor(Nfft*0.66),Nfft,Fs); %Nfft : number of points per segmentation. + default values : 50% overlap + Hamming window
	end
	
	[Nmic F Mw] = size(data_spec);
	for f=1:F
		CSM(:,:,f) = squeeze(data_spec(:,f,:)) * squeeze(data_spec(:,f,:))' /Mw; %mic signals
	end
	Mw
end
