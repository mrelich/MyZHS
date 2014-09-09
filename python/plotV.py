
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
outfile = TFile("rootfiles/"+beamE+"_test.root","recreate")

# Binning
nbins = 1000
xmin  = -3
xmax  = 3

V_perEvent = []
for i in range(nEvts):
    V_perEvent.append( TProfile("VP_evt"+str(i),"",nbins,xmin,xmax) )

V_avg = TProfile("VP_avg","",nbins,xmin,xmax)


#----------------------------------------#
# Now loop over and plot
#----------------------------------------#
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
        if "azimuth" in line: 
            eventFound = True
            continue
        if "#" in line and eventFound: 
            break 
        if "#" in line: 
            continue

        # Get the data
        pair = line.split()
        if len(pair) == 0: continue
        time = float(pair[0])
        Potential = abs(float(pair[1]))

        # Fill Data
        V_perEvent[iEvt].Fill(time,Potential)
        V_avg.Fill(time,Potential)
        if maximum < abs(Potential):
            maximum = abs(Potential)
    # Close input file
    infile.close()
    print "Event: ", iEvt, " has max ", maximum

# Now take the absolute value
#for i in range(1,V_avg.GetNbinsX()):
#    bc = V_avg.GetBinContent(i)
#    V_avg.SetBinContent(i,abs(bc))
#    
#    for hist in V_perEvent:
#        bc = hist.GetBinContent(i)
#        hist.SetBinContent(i,abs(bc))


# Save output
outfile.Write()
outfile.Close()

