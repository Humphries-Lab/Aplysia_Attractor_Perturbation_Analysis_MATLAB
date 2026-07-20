function spike_conv = fn_gaussConvSpikes(peaks, sigma_s, fs)
if sigma_s<(1/fs), sigma_s=1/fs; end
sigma_f=sigma_s*fs; half_width=ceil(3*sigma_f);
k_frames=-half_width:half_width;
gauss_kernel=single(exp(-0.5*(k_frames/sigma_f).^2));
gauss_kernel=gauss_kernel/sum(gauss_kernel);
[nN,nF]=size(peaks); spike_conv=zeros(nN,nF,'single');
CHUNK=20;
for c0=1:CHUNK:nN
    c1=min(c0+CHUNK-1,nN);
    for n=c0:c1
        spike_conv(n,:)=conv(single(peaks(n,:)),gauss_kernel,'same');
    end
end
spike_conv(spike_conv<0)=0;
end
