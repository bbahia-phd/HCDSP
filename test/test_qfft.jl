using Revise
using HCDSP
using Images, TestImages

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

@test isapprox(Qi,Q)
