# HyperComplex Digital Signal Processing (HCDSP)
I am basing my tests on 

                  https://github.com/JuliaGeometry/Quaternions.jl

An Official Quaternion Julia module that is part of JuliaGeometry to handle quaternions and etc. Used mainly for rotations, the way such a package is developed does not appeal to the way I am intending to use Julia. This other (unofficial) module

                  https://github.com/peakbook/Quaternions.jl

has some closer ideas, too. I will try to build upon them and the Quaternion Toolbox for MATLAB (Sangwine and Le Bihan) and develop a Julia-based Quaternion Tool Box that might be integrated with SeismicJulia.

I will do so by means of a Quaternion module called HCDSP. In HCDSP, we find the definition of a quaternion and some overwrites for Julia built-in functions to work for Quaternion typed variables. I will traceback the quaternion functions I need to translate, and work with them at first. For instance, taking the Fourier transform of a quaternion vector.
