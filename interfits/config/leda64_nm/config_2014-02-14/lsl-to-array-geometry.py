# -*- coding: utf-8 -*-
import pandas as pd
import ujson as json
import numpy as np
SOL = 299792458

a   = pd.read_csv('leda64nm-mapping.csv')
ant = pd.read_csv('leda64nm-antennas.csv')

a_sorted = a.sort('LEDA_input', ascending=1)
ant_sorted = ant.sort('STAND', ascending=1)

stands = a_sorted['stand']
stand_names = ["st_%i"%s for s in stands.values]

stabxyz = []
stabxyz_topo = []
delays = []

for s in stands:
    srow = ant[ant['STAND'] == s].values[0, 1:4]
    x,y,z = np.round(srow,3)
    
    lat2 = 34.070*np.pi/180.0
    trans1 = np.matrix([[0, -np.sin(lat2), np.cos(lat2)],
                        [1,  0,               0],
                        [0,  np.cos(lat2), np.sin(lat2)]])
    xyz = trans1*np.matrix([[x],[y],[z]])
    
    stabxyz.append(list(np.round(np.array(xyz).flatten(), 3)))
    stabxyz_topo.append((x,y,z))
    
    # Adding in delays
    sd = ant[ant['STAND'] == s][' DELAY'].values * SOL
    sd = sd.tolist()
    delays.append(sd)


print "TOPOCENTRIC:"
print "------------"    
for ll in range(len(stabxyz)):
    s = stabxyz_topo[ll]
    print "%s\t%2.3f\t%2.3f\t%2.3f"%(stand_names[ll], s[0], s[1], s[2])

print "\nECEF"
print "----"    
for ll in range(len(stabxyz)):
    s = stabxyz[ll]
    print "%s\t%2.3f\t%2.3f\t%2.3f"%(stand_names[ll], s[0], s[1], s[2])

#print stand_names
#print delays


json.dump(stand_names, open('json/stands-temp.json', 'w'))
json.dump(stabxyz_topo, open('json/stabxyz-topo-temp.json', 'w'))
json.dump(delays, open('json/delays-temp.json', 'w'))