function  [d,t,mx,my] = linear_events_3d_irregular_LRZ(mx,my,nt,dt,px,py,t0,A,f0,snr,L,seed,I);
% Computes 3D data cube composed of planes


nfft = (2^nextpow2(nt));

% mx = [0:1:nmx-1]*dmx;
% my = [0:1:nmy-1]*dmy;

% nmx = length(mx);
% nmy = length(my);

[nmx,nmy] = size(mx);

t = [0:1:nt-1]*dt;


if I==2
    mx = (mx/max(abs(mx(:)))).^I;  
    my = (my/max(abs(my(:)))).^I;
end


% Make steering vectors
D = zeros(nfft,nmx,nmy);
nevents = length(A);
wavelet = ricker(f0,dt);
nwavelet = length(wavelet);
tdelay = nwavelet*dt/2;
W = fft(wavelet,nfft);
for iw=2:nfft/2
    w = 2*pi*(iw-1)/(nfft*dt);
    T= zeros(nmx,nmy);
    for ie=1:nevents
%         smx = exp(-1i*w*mx*px(ie));
%         smy = exp(-1i*w*my*py(ie));
%         [MX,MY] = ndgrid(smx,smy);

        MX = exp(-1i*w*mx*px(ie));
        MY = exp(-1i*w*my*py(ie));

        T = T + A(ie)*exp(-1i*w*(t0(ie)-tdelay))*W(iw)*MX.*MY;
    end
    D(iw,:,:) =  T;
    D(nfft-iw+2,:,:) = conj(T);
end
d = ifft(D,[],1);
d = real(d(1:nt,:,:));

% Add noise

op = hamming(L);
randn('state',seed);
Noise = convn(randn(size(d)),op,'same');
sd = std(d(:));
sn = std(Noise(:));
a = (sd/sn)/sqrt(snr);
n  = a*Noise; d0 = d;
d = d + n;
snr = var(d0(:))/var(n(:));
return;
