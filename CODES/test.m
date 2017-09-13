for f=1:nb_F	
	%%% Direct problem : generate mic pressures
	Spp_retro(:,:,f) = G_src_mic(:,:,f)*Sqq_est(:,:,f)*G_src_mic(:,:,f)';
end