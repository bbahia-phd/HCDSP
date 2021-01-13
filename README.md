Multicomponent Seismic Data Processing using Quaternion Algebra

IT CAN BE HELPFUL TO ENSURE YOUR CODE STILL WORKS AFTER YOU MAKE CHANGES,
AND CAN BE USED WHEN DEVELOPING A WAY OF SPECIFYING THE BEHAVIORS YOUR CODE
SHOULD HAVE WHEN COMPLETE.

I am basing my tests on 

                  https://github.com/JuliaGeometry/Quaternions.jl

An Official Quaternion Julia module that is part of JuliaGeometry to handle quaternions and etc. 

Used mainly for rotations, the way it was coded does not appeal to the way I am intending to use Julia. This other (unofficial) module

                  https://github.com/peakbook/Quaternions.jl

has some closer ideas too. I will try to build upon them and the Quaternion Toolbox for MATLAB (Sangwine and Le Bihan) and develop a Julia-based Quaternion Tool Box that can be integrated with SeismicJulia.

I will do so by means of a Quaternion module called QSEIS. In QSEIS, we find the definition of a quaternion and some overwrites for Julia built-in functions to work for Quaternion typed variables. I will traceback the quaternion functions I need to translate, and work with them at first. For instance, taking the Fourier transform of a quaternion vector.


# I am learning this now.
# See https://docs.julialang.org/en/v1/stdlib/Test/
#
# This is a testset: A simple unit testing functionality, a way to see
#if your is correct by checking if the results are what you expect.

#IT CAN BE HELPFUL TO ENSURE YOUR CODE STILL WORKS AFTER YOU MAKE CHANGES,
#AND CAN BE USED WHEN DEVELOPING A WAY OF SPECIFYING THE BEHAVIORS YOUR CODE
#SHOULD HAVE WHEN COMPLETE.

#I am also basing my tests 

#https://github.com/JuliaGeometry/Quaternions.jl

#An Official Quaternion Julia module that is part of JuliaGeometry to handle quaternions and etc. Used mainly for rotations, the way it was coded does not appeal to the way I am intending to use Julia. This other (unofficial) module

#https://github.com/peakbook/Quaternions.jl

#has some closer ideas too. 

# I will try to build upon them and the Quaternion Toolbox for MATLAB (Sangwine and Le Bihan) and develop a Julia-based Quaternion Tool Box that can be integrated with SeismicJulia. I will do so by means of a Quaternion module called QSEIS (the one called Quaternions already exists, and I intended to call it QSEIS anyways HAHA). In QSEIS, we find the definition of a quaternion and some overwrites for Julia built-in functions to work for Quaternion typed variables. 

# I will traceback the quaternion functions I need to translate, and work with them at first. For instance, taking the Fourier transform of a quaternion vector.

# I need to test my module. This is where I do it. Prolixo.
