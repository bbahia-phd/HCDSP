function [dout] = epocs_lrz(d,dtrue,dt,f_low,f_high, option, perc_i, perc_f, x_irreg,x_reg,N);
%POCS: 3D (x,y,t) seismic data regularization using POCS
%
% [dout,converg] = pocs(d,dtrue,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);
%
%  IN   d:      data  a t,x,y cube)
%       dtrue:  true data prior to decimation (insert d is you do not have dtrue)
%       T:      sampling operator
%       f_low:  min freq in Hz in the data
%       f_high: max freq in Hz inthe data
%       option:    = 1 linear threshold schedulde
%                    2 exponential threshold schedulde
%                    3 data adaptive threshold schedule
%       perc_i: percentage of max amplitude for initial threshold
%       perc_f: percentage of max amplitude for final threshold (perc_f<<perc_i)

%  OUT  dout:  reconstructed cube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I try to modify the POCS code for reconstruction in time, in order to avoid  %%%
%% the interpolation operator calculation for each frequency.                   %%%
%% However it doesn't work very well under the POCS form                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 0.01;

[nt,nx1,nx2] = size(dtrue);
[nt,nk] = size(d);

nf  = 2^nextpow2(nt);
nk1 = 2^nextpow2(nx1);
nk2 = 2^nextpow2(nx2);


% we don't need the sampling operator, we use the interpolation operator
% for reconstruction

% S = 1.-T;

% Min and max freq indeces for the reconstruction

k_low = max(floor(f_low*dt*nf) + 1,2);
k_high = min(floor(f_high*dt*nf) + 1,nf/2);

Dout = zeros(nf,nx1,nx2);

% from irregular observation to regular 
d_reg = interpolator_Kaiser_sinc_3D(d,x_irreg,x_reg,3,10,1);
dold  = d_reg;

Dtrue = fft(dtrue,nf,1);


for iter = 1:N
    
    D = fft(dold,nf,1);

    for k = k_low:k_high

        xadj = squeeze(D(k,:));

        th = th_schedule(option,xadj,perc_i,perc_f,N);
        y = xadj;
     
        Y = fft2(y,nk1,nk2);
        A = abs(Y);
        Angle = angle(Y);
        I = find(A<th(iter));
        A(I) = 0;
        Y = A.*exp(1i*Angle);

        y = ifft2(Y);

        y = y(1:nx1,1:nx2);

  
        Dout(k,:,:) = y;
        Dout(nf-k+2,:,:)  = conj(y);
    end
    
    yout = ifft(Dout,[],1);
    yout = yout(1:nt,:,:);
    
    yy1 = interpolator_Kaiser_sinc_3D(yout,x_irreg,x_reg,3,10,0);
    yy2 = interpolator_Kaiser_sinc_3D(yy1,x_irreg,x_reg,3,10,1);
    
    
    dnew = d_reg + yout -yy2;
    dold = dnew;
    dif1 = dnew - yout; 
    c1 = sum( (abs(dif1(:))).^2);
    c  = sum( (abs(yout(:))).^2);
    
    fprintf('Iteration number = %6.0f, error = %10.4f\n',iter, c1/c);
    if c1/c < tol
        break
    end
    
end
       

dout = yout;


  
return

function th = th_schedule(option,x,perc_i,perc_f,N)
% Function used by pocs to define the threshold schedule based
% on parameter option
%
%  option == 1 --> linear
%  option == 2 --> exponential
%  option == 3 --> use real amplitude to define curve
%

[nx1,nx2] = size(x);
nk1 = 2^nextpow2(nx1);
nk2 = 2^nextpow2(nx2);
X = fft2(x,nk1,nk2);
A = abs(X);
Amax = max(A(:));
th_i = perc_i*Amax/100;
th_f = perc_f*Amax/100;

% th_i = perc_i/100;
% th_f = perc_f/100;

k = [1:1:N];

if option==1
    th = th_i + (th_f-th_i)*(k-1)/(N-1);
%     amp = sort(A(:));
%     th  = amp(round(floor(th*nk1*nk2)+1));
end

if option==2
    b  = -log(th_f/th_i);
    th = th_i*exp(-b*(k-1)/(N-1));
end

if option==3
    Amax = max(A(:));
    avec = sort(reshape(A,nk1*nk2,1),'descend');
    I = find(avec>perc_i*Amax/100);
    avec(I)=[];
    I = find(avec<perc_f*Amax/100);
    avec(I)=[];
    th = zeros(N,1);
    th(1) = avec(1);
    th(N) = avec(end);
    for j=2:N-1
        th(j)=avec(ceil((j-1)*length(avec)/(N-1)));
    end
end

return;
