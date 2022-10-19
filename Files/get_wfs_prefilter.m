function h = get_wfs_prefilter(length,fs)
            N = size(obj.virtual_source.source_signal.time_series,1);
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
            H = sqrt(1i*w/obj.c);
            H(end/2+1) = real(H(end/2+1));
            h = fftshift(ifft(H));
            h = h(round(end/2)-50+1:round(end/2)+50).*hann(100);
    %        h = 1;

end

