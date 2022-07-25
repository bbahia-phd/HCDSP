pwd()
cd(joinpath(homedir(),"projects/HCDSP/test/hyper_events"))

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using SeisMain, SeisPlot, SeisProcessing

