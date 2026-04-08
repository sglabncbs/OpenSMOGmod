from OpenSMOGmod import SBM
import sys

print ("This example code uses an addition input_addon.xml file for OpenSMOGmod forces")
#setting MD params
simul_prefix = "Output"
dt = 0.0010             #ps
collision_rate = 1.0    #ps^-1
r_cutoff = 3.0          #nm 
T = 0.5                 #reduced units 
sbm = SBM(name=simul_prefix, time_step=dt, collision_rate=collision_rate, r_cutoff=r_cutoff, temperature=T,pbc=True)

#platform="cuda"` or ="HIP" or ="opencl" or ="cpu"
sbm.setup_openmm(platform='opencl',GPUindex='default')

sbm.saveFolder("Output")
grofile = "Input/input_struc.gro"    
topfile = "Input/input_topol.top"
xmlfile = "Input/input_topol.xml"
modfile = "Input/input_addon.xml"
sbm.loadSystemFiles(Grofile=grofile, Topfile=topfile, Xmlfile=xmlfile,Modfile=modfile)

sbm.createSimulation()

trjformat = "xtc"
sbm.createReporters(trajectory=True,trajectoryFormat=trjformat, energies=True, energy_components=True, interval=1000, checkpoint=True, checkpointInterval=10000)

#report is for verbose
sbm.run(nsteps=50000000, report=True, interval=1*10**3)

