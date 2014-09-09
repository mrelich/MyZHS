
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# This script will loop over the runs in the specified output #
# directory and will plot the N(e+p) and N(e-p) distributions #
# as a function of radiation lengths.  This can then later be #
# compared to the distributions obtained from GEANT4          #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

from ROOT import TProfile, TFile, kBlack, kBlue
import sys, os

#---------------------------------#
# Get user inputs which should
# be the directory and number
# of events to consider        
#---------------------------------#

#
## specify specific order
#
argv = sys.argv
if len(argv) != 3:
    print "Not correct number of args"
    print "script.py <directory> <num events>"
    sys.exit()
#
## set directory
#
directory = str(argv[1])

if not os.path.exists(directory):
    print "Directory doesn't exist"
    print "Directory: ", directory
    sys.exit()

#
## Set the number of events
#
nEvent = int(argv[2])


#---------------------------------#
# Open output file to store hists
#---------------------------------#

#outfile = TFile("rootfiles/test.root","RECREATE")
outfile = TFile("rootfiles/"+directory+"_NPart.root","RECREATE")

#---------------------------------#
# Setup profiles here
#---------------------------------#

def makeProfile(name, xtitle, ytitle, nbins, xmin, xmax, color, marker):
    prof = TProfile(name,name,nbins,xmin,xmax)
    prof.GetXaxis().SetTitle(xtitle)
    prof.GetYaxis().SetTitle(ytitle)
    prof.SetLineColor(color)
    prof.SetMarkerColor(color)
    prof.SetMarkerStyle(marker)
    return prof

nbins  = 40 #20
xmin   = 0
xmax   = 40. #20.
xtitle = "Radiation Lengths"
ytitle = "Number of Particles"

p_SUM  = makeProfile("NPartSum",xtitle,ytitle,nbins,xmin,xmax,kBlack,20)
p_DIFF = makeProfile("NPartDiff",xtitle,ytitle,nbins,xmin,xmax,kBlue,20)

#---------------------------------#
# Useful formatting info
#---------------------------------#

def makeFName(base, event):
    if event < 10:
        return base + "00" + str(event) + ".EMP"
    elif event < 100:
        return base + "0" +str(event) + ".EMP"
    
    return base + str(event) + ".EMP"

def formatLine(line):
    l_split = line.split()
    if len(l_split) == 0:
        return [-1,-1,-1]
    return [float(l_split[0]),  # radiation length
            float(l_split[1]),  # N(e+p)
            float(l_split[3])   # N(e-p)
            ]
#---------------------------------#
# Loop over events and fill hists
#---------------------------------#
base = "DTe"
for i in range(0,nEvent):
    
    # Open next file
    newFile = open(makeFName(directory +"/" + base,i),"r")
    
    # Loop over file and save info
    for line in newFile:
        
        # Break out of loop when we get
        # to hybrid information
        if "Hybrid" in line: break

        # Don't include garbage
        if "#" in line: continue

        # Get info
        (rad, n_sum, n_diff) = formatLine(line)
        if rad == -1: continue

        # Now fill histograms
        bw = p_SUM.GetXaxis().GetBinWidth(1)
        p_SUM.Fill(rad-bw, n_sum)
        p_DIFF.Fill(rad-bw, n_diff)


# Save info
outfile.Write()
outfile.Close()
    
