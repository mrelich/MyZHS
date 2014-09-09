
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# Plot the vector potential for each individual shower #
# as well as save the average over N showers.          #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

import sys
from ROOT import TH1F, TFile, TProfile

#----------------------------------------#
# Input the directory
#----------------------------------------#

args = sys.argv
if len(args) != 3: 
    print "Input the directory to run over"
    print "and the number of events to consider"
    sys.exit()

beamE = args[1].split("/")[0]
nEvts = int(args[2])

#----------------------------------------#
# Specify histograms
#----------------------------------------#
outfile = TFile("rootfiles/"+beamE+"_Angular.root","recreate")
#outfile = TFile("test.root","recreate")

# Binning
nbins = 1000
xmin  = -10
xmax  = 10

# Plot the average VP for various angles.
# Open up one of the files and read in the
# angles that are stored
testFile = open(beamE+"/VPe000_T.EMP","r")
Angles = []
for line in testFile:
    if "Angle" in line:
        Angles.append(float(line.split()[-1]))

V_avgVsAngle = []
for ang in Angles:
    V_avgVsAngle.append( TProfile("VP_avg_"+str(ang),"",nbins,xmin,xmax) )

for prof in V_avgVsAngle:
    prof.Sumw2()
#----------------------------------------#
# Have method to determine which angle
# They aren't exactly 0, 10, 20, etc. but
# rather some delta from the cherenkov 
# angle at 55.8
#----------------------------------------#
def getAngle(line):
    # Get Angle
    angle = float(line.split()[-1])
    
    # Find closest match
    minimum = 1000
    mini = -1
    for i in range(len(Angles)):
        if minimum > abs(Angles[i] - angle):
            minimum = abs(Angles[i] - angle)
            mini    = i
    if mini < 0: 
        print "Error: ", mini
    return mini

#----------------------------------------#
# Now loop over and plot
#----------------------------------------#
skippedFirst = False
for iEvt in range(nEvts):

    # Open input file
    infileName = ""
    if iEvt < 10:
        infileName = "VPe00"+str(iEvt)+"_T.EMP"
    elif iEvt < 100:
        infileName = "VPe0"+str(iEvt)+"_T.EMP"
    else:
        print "Something wrong"
        print iEvt
        continue
    
    infile = open(beamE+"/"+infileName,"r")

    # Loop over the data
    maximum = -999
    eventFound = False
    for line in infile:
        if "Angle [deg]" in line:
            iAng = getAngle(line)
            continue
        if "#" in line: 
            continue

        # Get the data
        pair = line.split()
        if len(pair) == 0: continue
        time = float(pair[0])
        #Potential = abs(float(pair[1]))
        Potential = -1*float(pair[1])

        # Fill Data
        V_avgVsAngle[iAng].Fill(time,Potential)

    # Close input file
    infile.close()
    

# This seems to give very wrong results...
# I need to investigate what is actually going on
# so it turns out that internally the number of entries
# per bin are being stored.  So when I go back and manually
# set the bin entries. if I set it to 1 using setBinEntries(bin,1)
# it works.  This is too complicated.  Just take -1*A above!
#def absProf(prof):
#    nbins = prof.GetNbinsX()
#    for i in range(1,nbins):
#        bc = prof.GetBinContent(i)
#        be = prof.GetBinError(i)
#        prof.SetBinContent(i,abs(bc))
#        prof.SetBinError(i,be)
#    return prof

for prof in V_avgVsAngle:
    thebin = prof.FindBin(0)
    print thebin, prof.GetBinContent(thebin), prof.GetBinEntries(thebin)

# Save output
outfile.Write()
outfile.Close()

