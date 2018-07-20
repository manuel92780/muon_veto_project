#!/usr/bin/env python 
from icecube import icetray, dataclasses, simclasses, dataio, phys_services
from I3Tray import I3Tray
from icecube.icetray import I3Units

import numpy as np
import time
import argparse
parser = argparse.ArgumentParser(description = "filter simulation into manageable data sizes")
parser.add_argument('infiles', nargs='+')
parser.add_argument('-d','--debug', default=False, dest = 'BUG')
parser.add_argument('-s','--simulation', default='corsika', dest = 'SIM')
args = parser.parse_args()

#timing stuff
start_time = time.asctime()
print 'Started:', start_time

#load a long list of files
file_name = args.infiles
infile    = dataio.I3FrameSequence(file_name)
normalization  = len(file_name)
print "loaded your "+args.SIM+" MC in " + str(normalization) +" files!"
if(args.SIM == 'corsika'): normalization = 1


#initialize the arrays
event  = [];
event.append([]); #[0]  weight
event.append([]); #[1]  muon bundle size 

muon_1 = [];
muon_1.append([]);#[0] energy
muon_1.append([]);#[1] zenith
muon_1.append([]);#[2] azimuth
muon_1.append([]);#[3] x pos
muon_1.append([]);#[4] y pos
muon_1.append([]);#[5] z pos
muon_1.append([]);#[6] radial distance

muon_2 = [];
muon_2.append([]);#[0] energy
muon_2.append([]);#[1] zenith
muon_2.append([]);#[2] azimuth
muon_2.append([]);#[3] x pos
muon_2.append([]);#[4] y pos
muon_2.append([]);#[5] z pos
muon_2.append([]);#[6] radial distance

#initialize the loop
event_count = 1
bad_events = 0;
for frame in infile:
    if(event_count%100000 == 0):
        print "Event: "+ str(event_count);
        if(args.BUG): break;
    event_count+=1
    if("EnteringMuon_0" not in frame): continue; #requires at least one muon to enter detector
    weight = frame["Weight"].value/normalization
    number_muons = 0;
    energy_1 = 0; zenith_1 = 0; azimuth_1 = 0; xpos_1 = 0; ypos_1 = 0; zpos_1 = 0; rad_1 = 0;
    energy_2 = 0; zenith_2 = 0; azimuth_2 = 0; xpos_2 = 0; ypos_2 = 0; zpos_2 = 0; rad_2 = 0;
    for num in range(10000): #only count entering muons
        if("EnteringMuon_"+str(num) not in frame): break;
        number_muons +=1;
        #and sort at each step
        particle = frame["EnteringMuon_"+str(num)]; primary = frame["MCPrimary"];
        #particle_track = (frame["MMCTrackList"])[num];
        energy = particle.energy/ I3Units.GeV
        if(energy > energy_1):
            energy_2 = energy_1; zenith_2 = zenith_1; azimuth_2 = azimuth_1;
            xpos_2 = xpos_1; ypos_2 = ypos_1; zpos_2 = zpos_1;
            rad_2 = rad_1; 
            energy_1 = energy; zenith_1 = particle.dir.zenith; azimuth_1 = particle.dir.azimuth;
            xpos_1 = particle.pos.x; ypos_1 = particle.pos.y; zpos_1 = particle.pos.z;
            #rad_1 = phys_services.I3Calculator.closest_approach_distance(particle,primary.pos);
        elif(energy > energy_2):
            energy_2 = energy; zenith_2 = particle.dir.zenith; azimuth_2 = particle.dir.azimuth;
            xpos_2 = particle.pos.x; ypos_2 = particle.pos.y; zpos_2 = particle.pos.z;
            #rad_2 = phys_services.I3Calculator.closest_approach_distance(particle,primary.pos);

    if(weight > 100): 
        bad_events+=1
        continue
    event[0].append(weight); event[1].append(number_muons);
    muon_1[0].append(energy_1); muon_1[1].append(zenith_1); muon_1[2].append(azimuth_1); 
    muon_1[3].append(xpos_1); muon_1[4].append(ypos_1); muon_1[5].append(zpos_1); muon_1[6].append(rad_1);
    muon_2[0].append(energy_2);muon_2[1].append(zenith_2); muon_2[2].append(azimuth_2);
    muon_2[3].append(xpos_2); muon_2[4].append(ypos_2);muon_2[5].append(zpos_2); muon_2[6].append(rad_2);

print "looped over: " + str(event_count) + " events"
print "tossed: " + str(bad_events) + " events due to inf bug"

#dump into binary .npy file
data = [];
data.append(event);
data.append(muon_1);
data.append(muon_2);
np.save(args.SIM+"_array.npy",data)

#timing stuff
end_time = time.asctime()
print 'Ends:', end_time
