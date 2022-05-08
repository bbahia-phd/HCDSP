cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate("./dev/hcdsp/.")
Pkg.status()

using Revise

using PyPlot
using SeisMain, SeisPlot

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";  
# data dir home
data_path = "/dev/Breno_GOM/iso_vsp01/sgy/";  

# Convert to Seis files
SegyToSeis( joinpath(data_path,"sgy/iso_vsp01_zx.sgy"),
            joinpath(data_path,"iso_vsp01_zx.seis"));

SeisWindow( "./../iso_vsp01/sgy/iso_vsp01_zx.seis",
            "./../iso_vsp01/sgy/iso_vsp01_zx_crg1350.seis",
            key=["gelev"],
            minval=[-1350.0],
            maxval=[-1350.0]);            

# SeisWindow( joinpath(data_path,"iso_vsp01_zx.seis"),
#             joinpath(data_path,"iso_vsp01_zx_crg1350.seis"),
#             key=["gelev"],
#             minval=[-1350.0],
#             maxval=[-1350.0]);

SegyToSeis( joinpath(data_path,"sgy/iso_vsp01_zy.sgy"),
            joinpath(data_path,"iso_vsp01_zy.seis"));

SeisWindow( "./../iso_vsp01/sgy/iso_vsp01_zy.seis",
            "./../iso_vsp01/sgy/iso_vsp01_zy_crg1350.seis",
            key=["gelev"],
            minval=[-1350.0],
            maxval=[-1350.0]);  

# SeisWindow(joinpath(data_path,"iso_vsp01_zy.seis"),
#             joinpath(data_path,"iso_vsp01_zy_crg1350.seis"),
#             key=["gelev"],
#             minval=[-1350.0],
#             maxval=[-1350.0]);
 
SegyToSeis( joinpath(data_path,"sgy/iso_vsp01_zz.sgy"),
            joinpath(data_path,"iso_vsp01_zz.seis"));

SeisWindow( "./../iso_vsp01/sgy/iso_vsp01_zz.seis",
            "./../iso_vsp01/sgy/iso_vsp01_zz_crg1350.seis",
            key=["gelev"],
            minval=[-1350.0],
            maxval=[-1350.0]);  
            
# SeisWindow(joinpath(data_path,"iso_vsp01_zz.seis"),
#            joinpath(data_path,"iso_vsp01_zz_crg1350.seis"),
#            key=["gelev"],
#            minval=[-1350.0],
#            maxval=[-1350.0]);