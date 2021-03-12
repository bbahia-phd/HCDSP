using Revise
using HCDSP
using Images, TestImages, Test

# Get color image
img = testimage("mandrill");

# split in channels
ch = float(channelview(img));

# create a quaternion
Q = quaternion(ch[1,:,:,],ch[2,:,:],ch[3,:,:,]);

# Define eigenaxis for transformation
ax = normalize(quaternion(1.f0,1.f0,1.f0))

# Define a side
side = "left";

# qft
Qf = qfft(Q,ax,side);

# iqft
Qi = iqfft(Qf,ax,side);

@show norm(Qi-Q)^2
@test isapprox(Qi,Q)

# qft
Q1 = qfft(Q,ax,side,1:2);
Q2 = qfft(Q,ax,side)

@show norm(Q1-Q2)^2
@test isapprox(Q1,Q2)

# iqft
Q1i = iqfft(Q2,ax,side,1:2);
Q2i = iqfft(Q1,ax,side);

@show norm(Q1i-Q2i)^2
@test isapprox(Q1i,Q2i)

@show norm(Q1i-Q)^2
@test isapprox(Q1i,Q)

@show norm(Q2i-Q)^2
@test isapprox(Q2i,Q)
