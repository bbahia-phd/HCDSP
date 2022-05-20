% Compare the irregular reconstruction with POCS and extended POCS(EPCS)
% Linear Events --------------------------------------------------
% ----------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
%% without random noise 

close all; clear all; clc;

nx = 30;                    % number of traces in in-line direction
ny = 30;                    % number of trace in x-line direction
nt = 200;                   % number of time samples

dx = 20.;                   % intervals in space and time
dy = 20.;
dt = 4/1000;

%% add random jitter for X and Y direction
mx = [0:1:nx-1]*dx;
my = [0:1:ny-1]*dy;

[MX,MY] = ndgrid(mx,my);

rand('seed',6)

jx = (rand(size(MX))*2-1) * 1 *dx;
MX2 = MX + jx;
MX2(1,:) = abs(MX2(1,:)); % avoid the first coordinate is a negative value

jy = (rand(size(MY))*2-1) * 1 *dy;
MY2 = MY + jy;
MY2(:,1) = abs(MY2(:,1)); % avoid the first coordinate is a negative value

%% Linear events
px = [1/4000, 1/9000, -1/4000];          % ray parameters of 3D plane waves
py = [1/4000,-1/4000, -1/6000];

t0 = [0.1, 0.3, 0.6];                    % intercepts
A = [1,-1,1];                            % amplitudes

f0   =  20;                              % central frequency of source                          % S/N of additive noise
L    =  3;                               % Amount of smoothing of noise in time
seed =  21;
I    =  1.;                              % I:1 for linear events

%% Clean regular data----------------------
[d0,t,mx,my]    = linear_events_3d_LRZ(mx,my,nt,dt,px,py,t0,A,f0,99999999,L,seed,I);

%% -----------------Clean irregular data----------------------
[d1,t,mxx,myy] = linear_events_3d_irregular_LRZ(MX2,MY2,nt,dt,px,py,t0,A,f0,99999999,L,seed,I);
% --------------------------------------------------------

dobs0 = d0;
dobs = d1;     % dobs is 3D
               % d1 is 3D
               % d0 is 3D

%% Decimation 
perc1  = 0.5;
ndec1 = floor(perc1*nx*ny);
rand('seed',6)
indx1 = sort(randsample(nx*ny,ndec1));

dobs1 = reshape(dobs,nt,nx*ny,1);        % dobs1 is 2D
dobs0 = reshape(dobs0,nt,nx*ny,1);       % dobs1 is 2D

dobs1(:,indx1) = zeros;                  % dobs1 is 2D
dobs0(:,indx1) = zeros;                  % dobs1 is 2D

[T]  = GetSamplingOp(dobs1);

indxobs = setdiff([1:1:nx*ny]',indx1);   % define the index for binning

%% Erratic noise (if needed)
% percentage of erratic noise
perc2 = 0.3;

indxobs_er = randsample(length(indxobs),floor(length(indxobs)*perc2)); % define the index for adding erratic noise
 
% add erratic noise to onserved data
% dobs1 = SeisErraticNoise(dobs1,dt,floor(length(indxobs)*perc2));

dobs2   = dobs1(:,indxobs);       % extract the real irregular without binning

% ----------------------------------------------------------------------------------
% % ---------------------------------------------------------------------------------
% dobs2 series: real observed irregular data (dobs2,dobs22, etc) (dobs2 is 2D, dobs22 is 3D)
% dobs1 series: observed irregular data after binning (dobs1,dobs11, etc) (dobs1 is 2D, dobs11 is 3D)
% % ---------------------------------------------------------------------------------
% % ---------------------------------------------------------------------------------


%%
dobs1(:,indxobs) = dobs2; % binning dobs2

% reshape to a 2D matrix.
MX2obs = reshape(MX2,nx*ny,1);
MY2obs = reshape(MY2,nx*ny,1);

% find the observed irregular coordinate
% get rid of the traces with zero coordinate
mxobs = MX2obs(indxobs);
myobs = MY2obs(indxobs);

h = zeros(nx,ny,2);
h(:,:,1) = MX;
h(:,:,2) = MY;

htt = zeros(length(mxobs),2);
htt(:,1) = mxobs;
htt(:,2) = myobs;


% reshape the bined irregular data into a 3D cube for plotting

dobs11 = reshape(dobs1,nt,nx,ny);   % dobs11 is 3D
dobs00 = reshape(dobs0,nt,nx,ny);   % dobs11 is 3D



%% Spectrum Plot

d11slice = squeeze(dobs11(:,20,:));
[S,k,f] = fk_spectra(d11slice,dt,dx,6); 
figure;imagesc(k,f,S);
title('observed')
xlabel('Wavenumber');
ylabel('Frequency');

d00slice = squeeze(dobs00(:,20,:));
[S,k,f] = fk_spectra(d00slice,dt,dx,6); 
figure;imagesc(k,f,S);
title('observed')
xlabel('Wavenumber');
ylabel('Frequency');

%%
T11 = reshape(T,nt,nx,ny);
TT = squeeze(T11(1,:,:));
flow  =  1; 
fhigh = 80;
option = 1;
perc_i = 99;
perc_f = 1;
N = 100;
a = 1;
tol = 0.001;


dtrue = zeros(size(dobs00)); 

% interpolator_Kaiser_sinc_3D(U,x_irreg,x_reg,N,a,type)
% x0 = randn(size(dobs));
% value = power_method(x0,@interpolator_Kaiser_sinc_3D,htt,h,N,a);
% value = power_method(x0,Hop,PARAM);

%% For Regular POCS interpolator, the estimation is 3D


%% 1- FOR POCS
% [dout1,e1,e2,freq] = pocs_lrz(dobs11,dtrue,TT,dt,flow,fhigh, option, perc_i, perc_f, N, a, tol);
% 
% QS11 = quality(dout1,d0)
% 
% % 1- POCS error
% edout1 = dout1-d0;
% 
% pp =0.5;
% hFig = figure;
% set(hFig, 'Position', [100 100 1200 600]);
% ntr = 32;
% subplot(141);imagesc(squeeze(d0(:,ntr,:)),[-pp,pp]);colormap(gray)
% subplot(142);imagesc(squeeze(dobs11(:,ntr,:)),[-pp,pp]);colormap(gray)
% subplot(143);imagesc(squeeze(dout1(:,ntr,:)),[-pp,pp]);colormap(gray)
% subplot(144);imagesc(squeeze(d0(:,ntr,:)-dout1(:,ntr,:)),[-pp,pp]);colormap(gray)
% title('Regular POCS reconstruction (X direction)')


%% 2- For EPOCS-freq
tic

[dout2,ee1,ee2,freq2] = epocs_freq_lrz(dobs2,dtrue,dt,flow,fhigh, option, perc_i, perc_f, htt, h,N);

time_epocs = toc

QS12 = quality(dout2,d0)

% 2 - EPOCS error
edout2 = dout2-d0;

pp = 0.5;
ta = (0:1:nt-1)*dt;
ya = [1:1:ny]*dy;

hFig = figure;
set(hFig, 'Position', [100 100 1200 600]);
ntr = 12;
subplot(141);imagesc(squeeze(d0(:,ntr,:)),[-pp,pp]);colormap(gray)
subplot(142);imagesc(squeeze(dobs11(:,ntr,:)),[-pp,pp]);colormap(gray)
subplot(143);imagesc(squeeze(dout2(:,ntr,:)),[-pp,pp]);colormap(gray)
subplot(144);imagesc(squeeze(d0(:,ntr,:)-dout2(:,ntr,:)),[-pp,pp]);colormap(gray)
title('Irregular EPOCS reconstruction (X direction)')


hFig = figure;
set(hFig, 'Position', [100 100 1200 600]);
subplot(141);wigb(squeeze(d0(:,ntr,:)),1,ya,ta,pp);colormap(gray)
subplot(142);wigb(squeeze(dobs11(:,ntr,:)),1,ya,ta,pp);colormap(gray)
subplot(143);wigb(squeeze(dout2(:,ntr,:)),1,ya,ta,pp);colormap(gray)
subplot(144);wigb(squeeze(d0(:,ntr,:)-dout2(:,ntr,:)),1,ya,ta,pp);colormap(gray)
title('Irregular EPOCS reconstruction (X direction)')




