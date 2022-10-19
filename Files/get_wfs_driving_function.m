function [ amp, delay, focused, AAfilt ] = get_wfs_driving_function( xs, x0, n0, fs, c, antialiasing )

k = bsxfun( @minus, x0, xs );
rho_P = sqrt(sum(k.^2,2));
rho_G = sqrt(sum(x0.^2,2));
R0 = mean(rho_G);

k_P = bsxfun(@times, k,1./rho_P);
k_n = sum(k_P.*n0,2);
k_t = sqrt(1-k_n.^2);
d0 = R0*(k_n + 1i*k_t);
dl = sqrt(sum((x0-circshift(x0,1)).^2,2));

if ~all(k_n<0)

    focused = -1;
    % Non-focused case
    win = double(k_n>=0);
    dref2 = d0.*rho_P./(rho_P + d0);
    % dref = rho_G.*rho_P./(rho_P + rho_G);
    amp = -focused*win.*k_n/sqrt(2*pi).*sqrt(dref2)./rho_P.*dl;
    delay = rho_P / c ;
else
    focused = 1;
    center = [0, 0];
    k_s = (xs-center)/norm(xs-center);

    win1 = k_P*k_s'.*((k_P*k_s')>0);
    win2 = (k_P*k_s'+1)/2;

    dref1 = rho_P.*rho_G./(rho_G - rho_P);
    dref2 = rho_P.*d0./(d0 + rho_P);

    amp = -win1.*k_n/sqrt(2*pi).*sqrt(dref1)./rho_P.*dl;
    amp = -win2.*k_n/sqrt(2*pi).*sqrt(dref2)./rho_P.*dl;

    delay = - rho_P / c ;
end

switch antialiasing
    case 'on'
        if focused == -1
            wc = pi./dl.*340./abs(sqrt(1-k_n.^2));
            N = 960;
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
            Nbut = 8;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = ifft(transfer,[],1);
        else
            wc = pi./dl.*340./abs(sqrt(1-k_n.^2));
            N = 960;
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
            Nbut = 8;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = ifft(transfer,[],1);
        end
    case 'off'
        AAfilt = [];
end

