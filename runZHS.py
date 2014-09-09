#!/usr/bin/env python

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# The code relies on changing an input parameter in   #
# a template file. The problem is I wanted to use sed #
# in order to modify the template file, but there is  #
# a random number involved. I don't know the validity #
# or how much 'randomness' properties exist from bash #
# $Random, so I will use TRandom3 object from ROOT in #
# order to submit jobs.                               #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

from ROOT import TRandom3
from subprocess import Popen
import os
import stat

#---------------------------------#
# User specifications here
#---------------------------------#

#NEvents   = 100 #00
NEvents   = 50 #00
NPrimary  = 1 #5000
#NPrimary  = 1
#beamE     = "10e3" #"1e6" # MeV
beamE     = "100e3" # MeV -- 100 GeV
#beamE     = "1e6" # MeV
outDir    = "beam"+beamE+"MeV_"+str(NPrimary)+"Prim_"+str(NEvents)+"NEvt/"
#outDir    = "testing"+beamE+"MeV/"
jobOutput = "output_beam"+beamE+".txt"
radLengths= 0.92 * 40 * 40

print radLengths

#
## Useful method
#
def n():
    return "\n"

#
## Some variable defs
#
random    = TRandom3()
tempName  = "template_"+beamE+".inp"
beamName  = "beam_"+beamE+".inp"

#---------------------------------#
# Method to build beam template.
# Right now just single particle. 
#---------------------------------#
def beamTemplate(energy):
    out = open(beamName,"w")
    #for i in range(0,1000):
    #for i in range(0,5000):
    for i in range(0,NPrimary):
        out.write("1. ")
        out.write("1 ")
        out.write(energy+" ")   # particle energy
        out.write("0. ")              # x-position
        out.write("0. ")              # y-position
        out.write("0. ")              # z-position
        out.write(n())
    out.close()
    return

#---------------------------------#
# Define method to build template
# taking as input whatever needs  
# to be configurable. This is based
# on template_ice.inp
#---------------------------------#
def buildTemplate(randNum, outName):
    out = open(tempName,"w")
    out.write("ice.med"+n())
    out.write(str(randNum)+n())
    out.write("-1"+n())
    out.write("0 0 0"+n())
    out.write("1 0 1"+n())
    out.write("1 1"+n())
    out.write("0.01"+n())
    out.write("100000."+n())
    out.write("1"+n())
    out.write("0. 0."+n())
    out.write("1"+n())
    out.write(str(radLengths) +n()) 
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write(outName+n())
    out.write("grid_frequency.inp"+n())
    #out.write("7 \t0. \t0.1 \t-0.1 \t0.3 \t-0.3 \t5.0 \t-5.0"+n())
    #out.write("11 \t0. \t-55.8 \t-45.8 \t-35.8 \t-25.8 \t-15.8 \t-5.8 \t4.2 \t14.2 \t24.2 \t33"+n())
    out.write("11 \t0. \t0.1 \t0.2 \t0.3 \t0.4 \t0.5 \t0.6 \t0.7 \t0.8 \t0.9 \t1"+n())
    out.write("1 0."+n())
    out.write("1"+n())
    #out.write("0.12  0.12 0.0"+n())
    out.write("554.4  0.0 924.0"+n())
    out.write("0.0  0.0  1.0"+n())
    out.write("-1.0  0.0  0.0"+n())
    #out.write("template_beam.inp"+n())
    out.write(beamName + n())
    out.close()
    #os.chmod(tempName,stat.S_IEXEC)
    print "Built template: ", tempName
    return

#---------------------------------#
# Loop and execute jobs
#---------------------------------#

# Check for job output location
if os.path.exists(jobOutput):
    os.remove(jobOutput)

# Make sure out directory exists
if not os.path.exists(outDir):
    os.mkdir(outDir)

# Generate Beam template
beamTemplate(beamE)

# Loop
for evt in range(0,NEvents):

    if evt % 20 == 0:
        print "Processing event: ", evt

    # Get event name and seed
    eventName = "e"
    if evt < 10:
        eventName += "00"
    elif evt < 100:
        eventName += "0"
    eventName += str(evt)

    seed = int( 1000000 * random.Rndm() )
    
    # Build the template
    buildTemplate(seed, eventName)
    
    # Execute
    command = "./ZHS_time_freq_Fresnel_beam < " + tempName + " >> " + jobOutput
    run     = Popen([command],shell=True)
    run.wait()
    os.remove(tempName)

    # Move output files
    mv = Popen(["mv *.EMP " + outDir + "&& mv *.out " + outDir],shell=True)
    mv.wait()

# Now remove beam template
os.remove(beamName)
