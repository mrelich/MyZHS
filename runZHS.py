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

#---------------------------------#
# User specifications here
#---------------------------------#

NEvents = 10
random  = TRandom3()
jobOutput = "output.txt"

#
## Useful method
#
def n():
    return "\n"


#---------------------------------#
# Define method to build template
# taking as input whatever needs  
# to be configurable. This is based
# on template_ice.inp
#---------------------------------#

def buildTemplate(randNum, outName):
    out = open("template.inp","w")
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
    out.write("92.4"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write("1"+n())
    out.write(outName+n())
    out.write("grid_frequency.inp"+n())
    out.write("5  0. 0.1 -0.1 0.2 -0.2"+n())
    out.write("1 0."+n())
    out.write("1"+n())
    out.write("0.12  0.12 0.0"+n())
    out.write("0.0  0.0  1.0"+n())
    out.write("-1.0  0.0  0.0"+n())
    out.write("template_beam.inp"+n())
    out.close()

#---------------------------------#
# Loop and execute jobs
#---------------------------------#

if os.path.exists(jobOutput):
    os.remove(jobOutput)

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
    command = "./ZHS_time_freq_Fresnel_beam < template.inp >> " + jobOutput
    run     = Popen([command],shell=True)
    run.wait()
    os.remove("template.inp")

    # Move output files
    mv = Popen(["mv *.EMP output/ && mv *.out output/"],shell=True)
    mv.wait()
