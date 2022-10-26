function [theta_out, hrtf_out] = append_hrtf(theta_orig, hrtf_orig, hrtf_interp, theta_interp )

[~,ix] = min( abs( theta_orig - theta_interp ) );
ix = ix + (sign(theta_interp-theta_orig(ix))>0);

hrtf_out = zeros(size(hrtf_orig,1)+1,size(hrtf_orig,2),size(hrtf_orig,3));

hrtf_out(1:ix-1,:,:) = hrtf_orig(1:ix-1,:,:);
hrtf_out(ix,:,:) = hrtf_interp;
hrtf_out(ix+1:end,:,:) = hrtf_orig(ix:end,:,:);
theta_out = [theta_orig(1:ix-1);theta_interp;theta_orig(ix:end)];

end

