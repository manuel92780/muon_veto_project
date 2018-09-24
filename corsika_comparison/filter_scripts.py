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
parser.add_argument('--emin', default=1e3, dest = 'EMIN')
parser.add_argument('--emax', default=1e5, dest = 'EMAX')
args = parser.parse_args()

#timing stuff
start_time = time.asctime()
print 'Started:', start_time

#calculate the radial distance 
def find_distance(x1,y1,z1,x2,y2,z2):
    dist = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return dist

#load a long list of files
file_name = args.infiles
infile    = dataio.I3FrameSequence(file_name)
normalization  = len(file_name)
print "loaded your "+args.SIM+" MC in " + str(normalization) +" files!"
normalization = 1.0


#initialize the arrays
event  = [];
event.append([]); #[0]  weight
event.append([]); #[1]  muon bundle size 
event.append([]); #[2]  sum N muons energy
event.append([]); #[3]  weight^2
event.append([]); #[4]  mistake bin
event.append([]); #[5]  primary type
event.append([]); #[6]  primary energy
event.append([]); #[7]  Radius 12, needs at least 2
event.append([]); #[8]  Radius 13, needs at least 3
event.append([]); #[9]  Radius 23, needs at least 3
event.append([]); #[10]  Radius 14, needs at least 4
event.append([]); #[11]  Radius 24, needs at least 4
event.append([]); #[12]  Radius 34, needs at least 4

muon_1 = [];
muon_1.append([]);#[0] energy
muon_1.append([]);#[1] zenith
muon_1.append([]);#[2] azimuth
muon_1.append([]);#[3] x pos
muon_1.append([]);#[4] y pos
muon_1.append([]);#[5] z pos

muon_2 = [];
muon_2.append([]);#[0] energy
muon_2.append([]);#[1] zenith
muon_2.append([]);#[2] azimuth
muon_2.append([]);#[3] x pos
muon_2.append([]);#[4] y pos
muon_2.append([]);#[5] z pos

muon_3 = [];
muon_3.append([]);#[0] energy
muon_3.append([]);#[1] zenith
muon_3.append([]);#[2] azimuth
muon_3.append([]);#[3] x pos
muon_3.append([]);#[4] y pos
muon_3.append([]);#[5] z pos

muon_4 = [];
muon_4.append([]);#[0] energy
muon_4.append([]);#[1] zenith
muon_4.append([]);#[2] azimuth
muon_4.append([]);#[3] x pos
muon_4.append([]);#[4] y pos
muon_4.append([]);#[5] z pos 

#initialize the loop
event_count = 1
bad_events = 0; no_mu_events = 0;
for frame in infile:
    if(event_count%100000 == 0):
        print "Event: "+ str(event_count);
        if(args.BUG): break;
    event_count+=1
    if("EnteringMuon_0" not in frame): continue; #requires at least one muon to enter detector
    weight = frame["Weight"].value/normalization
    prim_energy = frame["MCPrimary"].energy/ I3Units.GeV;
    prim_type   = 0;#frame["MCPrimary"].type;
    number_muons = 0; sum_energy = 0; 
    rad_12 = 0; rad_13 = 0; rad_14 = 0;
    rad_23 = 0; rad_24 = 0; rad_34 = 0;
    energy_1 = 0; zenith_1 = 0; azimuth_1 = 0; xpos_1 = 0; ypos_1 = 0; zpos_1 = 0;
    energy_2 = 0; zenith_2 = 0; azimuth_2 = 0; xpos_2 = 0; ypos_2 = 0; zpos_2 = 0;
    energy_3 = 0; zenith_3 = 0; azimuth_3 = 0; xpos_3 = 0; ypos_3 = 0; zpos_3 = 0;
    energy_4 = 0; zenith_4 = 0; azimuth_4 = 0; xpos_4 = 0; ypos_4 = 0; zpos_4 = 0;
    for num in range(10000): #only count entering muons
        if("EnteringMuon_"+str(num) not in frame): break;
        number_muons +=1;
        #and sort at each step
        particle = frame["EnteringMuon_"+str(num)];
        energy = particle.energy/ I3Units.GeV
        sum_energy += energy
        if(energy > energy_1):
            energy_4 = energy_1; zenith_4 = zenith_1; azimuth_4 = azimuth_1;
            xpos_4 = xpos_1; ypos_4 = ypos_1; zpos_4 = zpos_1;
            energy_3 = energy_1; zenith_3 = zenith_1; azimuth_3 = azimuth_1;
            xpos_3 = xpos_1; ypos_3 = ypos_1; zpos_3 = zpos_1;
            energy_2 = energy_1; zenith_2 = zenith_1; azimuth_2 = azimuth_1;
            xpos_2 = xpos_1; ypos_2 = ypos_1; zpos_2 = zpos_1;
            energy_1 = energy; zenith_1 = particle.dir.zenith; azimuth_1 = particle.dir.azimuth;
            xpos_1 = particle.pos.x; ypos_1 = particle.pos.y; zpos_1 = particle.pos.z;
        elif(energy > energy_2):
            energy_4 = energy_1; zenith_4 = zenith_1; azimuth_4 = azimuth_1;
            xpos_4 = xpos_1; ypos_4 = ypos_1; zpos_4 = zpos_1;
            energy_3 = energy_1; zenith_3 = zenith_1; azimuth_3 = azimuth_1;
            xpos_3 = xpos_1; ypos_3 = ypos_1; zpos_3 = zpos_1;
            energy_2 = energy; zenith_2 = particle.dir.zenith; azimuth_2 = particle.dir.azimuth;
            xpos_2 = particle.pos.x; ypos_2 = particle.pos.y; zpos_2 = particle.pos.z;
        elif(energy > energy_3):
            energy_4 = energy_1; zenith_4 = zenith_1; azimuth_4 = azimuth_1;
            xpos_4 = xpos_1; ypos_4 = ypos_1; zpos_4 = zpos_1;
            energy_3 = energy; zenith_3 = particle.dir.zenith; azimuth_3 = particle.dir.azimuth;
            xpos_3 = particle.pos.z; ypos_3 = particle.pos.y; zpos_3 = particle.pos.z;
        elif(energy > energy_4):
            energy_4 = energy; zenith_4 = particle.dir.zenith; azimuth_4 = particle.dir.azimuth;
            xpos_4 = particle.pos.z; ypos_4 = particle.pos.y; zpos_4 = particle.pos.z;
    if(sum_energy > args.EMAX or sum_energy < args.EMIN): continue
    if(weight > 100): 
        bad_events+=1
        continue

    if(number_muons > 3):
        rad_12 = find_distance(xpos_1, ypos_1, zpos_1, xpos_2, ypos_2, zpos_2)
        rad_13 = find_distance(xpos_1, ypos_1, zpos_1, xpos_3, ypos_3, zpos_3)
        rad_14 = find_distance(xpos_1, ypos_1, zpos_1, xpos_4, ypos_4, zpos_4)
        rad_23 = find_distance(xpos_2, ypos_2, zpos_2, xpos_3, ypos_3, zpos_3)
        rad_24 = find_distance(xpos_2, ypos_2, zpos_2, xpos_4, ypos_4, zpos_4)
        rad_34 = find_distance(xpos_3, ypos_3, zpos_3, xpos_4, ypos_4, zpos_4)
    elif(number_muons > 2):
        rad_12 = find_distance(xpos_1, ypos_1, zpos_1, xpos_2, ypos_2, zpos_2)
        rad_13 = find_distance(xpos_1, ypos_1, zpos_1, xpos_3, ypos_3, zpos_3)
        rad_23 = find_distance(xpos_2, ypos_2, zpos_2, xpos_3, ypos_3, zpos_3)
    elif(number_muons > 1):
        rad_12 = find_distance(xpos_1, ypos_1, zpos_1, xpos_2, ypos_2, zpos_2)

    event[0].append(weight); event[1].append(number_muons); event[2].append(sum_energy); 
    event[3].append(weight*weight); event[4].append(rad_12);
    event[5].append(prim_type); event[6].append(prim_energy);
    event[7].append(rad_12); event[8].append(rad_13); event[9].append(rad_23);
    event[10].append(rad_14); event[11].append(rad_24); event[12].append(rad_34);
    muon_1[0].append(energy_1); muon_1[1].append(zenith_1); muon_1[2].append(azimuth_1); 
    muon_1[3].append(xpos_1); muon_1[4].append(ypos_1); muon_1[5].append(zpos_1);
    muon_2[0].append(energy_2);muon_2[1].append(zenith_2); muon_2[2].append(azimuth_2);
    muon_2[3].append(xpos_2); muon_2[4].append(ypos_2);muon_2[5].append(zpos_2);

print "looped over: " + str(event_count) + " events"
print "tossed " + str(no_mu_events) + " no muon events"
print "tossed " + str(bad_events) + " events due to inf bug"

#dump into binary .npy file
data = [];
data.append(event);
data.append(muon_1);
data.append(muon_2);
np.save(args.SIM+"_array.npy",data)

#timing stuff
end_time = time.asctime()
print 'Ends:', end_time
