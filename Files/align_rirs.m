function [rir0_full] = align_rirs(theta0,r0,rir0,rmax,c,fs)
[~,ix0] = min(abs(theta0));
[~,ix_peak] = max(sum(squeeze(rir0(ix0,:,:)),1));
d_peak = ix_peak/fs*c;
n_corr_l = round((r0-d_peak)/c*fs);
n_max = round(rmax/c*fs);
n_corr_u  = n_max - size(rir0,3) - n_corr_l;
rir0_full = zeros( size(rir0,1),size(rir0,2),n_max );
for n = 1 : size(rir0,1)
    rir0_full(n,:,:) = [zeros(size(rir0,2),n_corr_l),squeeze(rir0(n,:,:)),zeros(size(rir0,2),n_corr_u)];
end

end

