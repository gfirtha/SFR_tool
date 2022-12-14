function [hrtf_extrap] = extrapolate_hrtf(HRTF_meas,fs,x_extrap,x_measurement)

r0 = mean(sqrt(sum(x_measurement.^2,2)));
Nt = size(HRTF_meas,3);

f = (0:Nt-1)/Nt*fs;

HRTF_meas = HRTF_meas(:,:,1:end/2+1);
f = f(1:end/2+1);

[ amp, delay, focused, AA] = get_wfs_driving_function( x_extrap, x_measurement, -x_measurement/r0, fs, 340, 'on');
[~,ix] = max(abs(amp));
%delay = delay - delay(ix)+Nt/5/fs;

AA = fft(AA,Nt,1).';
AA = AA(:,1 : end/2+1);


w = 2*pi*f;
H_prefilter = (-focused*1i*w/340).^0.5;

HRTF_prefiltered = bsxfun(@times, HRTF_meas, reshape(H_prefilter,[1,1,length(H_prefilter)]));
AA = reshape(AA,[size(amp,1),1,size(w,2)]);

H_wfs =  reshape(bsxfun(@times, exp(-1i*bsxfun(@times,w,delay)), amp),[size(amp,1),1,size(w,2)]);
H_wfs = bsxfun( @times, AA, H_wfs );

HRTF_extrap = squeeze(sum(bsxfun(@times, H_wfs, HRTF_prefiltered),1 ));
hrtf_extrap = [HRTF_extrap conj(fliplr(HRTF_extrap(:,2:end-1)))];

end

