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

NEvents   = 1 #00
beamE     = "3e12" # MeV
outDir    = "beam"+beamE+"MeV/"
jobOutput = "output_beam"+beamE+".txt"
radLengths= 0.92 * 40 * 100

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
    out.write("0.1"+n())
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
    out.write("5 \t0. \t0.1 \t-0.1 \t0.3 \t-0.3"+n())
    out.write("1 0."+n())
    out.write("1"+n())
    #out.write("0.12  0.12 0.0"+n())
    out.write("3.696  3.696 0.0"+n())
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
