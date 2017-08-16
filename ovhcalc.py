#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ovhcalc.py
"""Calculate overlap volume histogram (OVH) from DICOM RT Structure data."""
# Copyright (c) 2016 gluce

from dicompylercore import dicomparser, dvh, dvhcalc
import numpy as np
import pylab as pl
import math

#import matplotlib.pyplot as plt
def find_nearest_index(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return idx

def point_inside_polygon(x,y,poly):

    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def calculate_ovh(structure, ptv_roi, non_ptv_roi, wplt= False):

    #get the RT structures
    rtss = dicomparser.DicomParser(structure)

    structures = rtss.GetStructures()

    #Calculate ovh for non ptv organs to ptv organs
    #PTV information
    ss = structures[ptv_roi]
    ss['planes'] = rtss.GetStructureCoordinates(ptv_roi)
    keyss = ss['planes'].keys()
    z_ptv_min = min([float(j) for j in keyss])
    z_ptv_max = max([float(j) for j in keyss])
    ptv_points=[]
    for j in keyss:
        ptv_points.extend(ss['planes'][j][0]['data'])
    #Non-PTV Information
    s = structures[non_ptv_roi]
    total_volume = 0
    s['planes'] = rtss.GetStructureCoordinates(non_ptv_roi)
    key = s['planes'].keys()
    ovh = [0]
    min_Dis = 0
    max_Dis = 0
    for j in key:
        z = s['planes'][j][0]['data'][0][2]
	#Find if the point is in the PTV
        if (z >= z_ptv_min) and (z <= z_ptv_max):
            second_check = True
            polyss=[]
            for k in range(0,len(ss['planes'][keyss[find_nearest_index([float(jss) for jss in keyss],z)]][0]['data'])):
                polyss.append([ss['planes'][keyss[find_nearest_index([float(jss) for jss in keyss],z)]][0]['data'][k][0]\
                                  ,ss['planes'][keyss[find_nearest_index([float(jss) for jss in keyss], z)]][0]['data'][k][1]])
        else:
            second_check = False
            polyss = []

        x_min = s['planes'][j][0]['data'][0][0]
        x_max = s['planes'][j][0]['data'][0][0]
        y_min = s['planes'][j][0]['data'][0][1]
        y_max = s['planes'][j][0]['data'][0][1]

        for k in range(0,len(s['planes'][j][0]['data'])):
            if x_min > s['planes'][j][0]['data'][k][0]:
                x_min = s['planes'][j][0]['data'][k][0]
            if x_max < s['planes'][j][0]['data'][k][0]:
                x_max = s['planes'][j][0]['data'][k][0]
            if y_min > s['planes'][j][0]['data'][k][1]:
                y_min = s['planes'][j][0]['data'][k][1]
            if y_max < s['planes'][j][0]['data'][k][1]:
                y_max = s['planes'][j][0]['data'][k][1]
            del s['planes'][j][0]['data'][k][2]
        poly = s['planes'][j][0]['data']

        x_min = int(x_min)
        y_min = int(y_min)
        x_max = int(x_max)
        y_max = int(y_max)

        if (x_min != x_max) and (y_min != y_max):
            for x in range(x_min, x_max):
                for y in range(y_min, y_max):
                    if point_inside_polygon(x,y,poly):
                           #calculate minimum distance Dis to ptv_x
                        Dis = 999999
                        for l in range(0,len(ptv_points)):
                            dist=((x - ptv_points[l][0])**2+(y - ptv_points[l][1])**2+(z - ptv_points[l][2])**2)
                            if (dist < Dis): Dis = dist
                        #deter whether the Dis is minus
			Dis = math.sqrt(Dis)

                        if second_check and point_inside_polygon(x,y,polyss):
                            Dis = - int(Dis)
                        else:
                            Dis = int(Dis)

                        if Dis > max_Dis:
                            ovh.extend(np.zeros((Dis-max_Dis,), dtype=np.int))
                            max_Dis = Dis
                        if Dis < min_Dis:
                            ovh=np.zeros((min_Dis-Dis,), dtype=np.int).tolist() + ovh
                            min_Dis = Dis

                        ovh[Dis - min_Dis] = ovh[Dis - min_Dis] + 1

                        total_volume = total_volume + 1

    if total_volume != 0 :
        for j in range(1, max_Dis-min_Dis+1):
            ovh[j]=ovh[j]+ovh[j-1]

        ovh=[float(j) *100 / total_volume for j in ovh]

        if wplt:
            #print total_volume
            #print ovh
            pl.plot(np.linspace(min_Dis, max_Dis, max_Dis- min_Dis + 1), ovh)
	    pl.title('OVH graph for %s & %s'%(ss['name'],s['name']))
	    pl.xlabel('Distance ()')
	    pl.ylabel('Percentage Volume')
            pl.savefig('OVH_%s_%s.png'%(ss['name'],s['name']))
            print('OVH_%s_%s.png is generated!'%(ss['name'],s['name']))
            
        return {'Graph':ovh, 'max_Dis':max_Dis,'min_Dis':min_Dis}
    else:
	print 'OVH is empty'
        return {'Graph':[], 'max_Dis':0,'min_Dis':0}


# ========================== Test OVH Calculation =========================== #

def main():
    	# Read the example RT structure and RT dose files
        rtssfile = "testdata/rtss.dcm"

	RTss = dicomparser.DicomParser(rtssfile)
	RTstructures = RTss.GetStructures()

	ptv_rois = []
	non_ptv_rois = []


	for i in range(1, len(RTstructures)):
	    if 'ptv' in RTstructures[i]['name'].lower():
		ptv_rois.append(RTstructures[i]['id'])
	    else:
		non_ptv_rois.append(RTstructures[i]['id'])
	print '-------PTV-------'
	for i in range(0, len(ptv_rois)):
	    print ("%d:")%RTstructures[ptv_rois[i]]['id']+ " "+ RTstructures[ptv_rois[i]]['name'] 
	print '-----Non-PTV-----'
	for i in range(0, len(non_ptv_rois)):
	    print ("%d:")%RTstructures[non_ptv_rois[i]]['id'] + " " +RTstructures[non_ptv_rois[i]]['name'] 

	var1 = raw_input("Please enter OVH target id: ")
	var2 = raw_input("Please enter OVH organ id: ")

	print 'select PTV %s and Organs %s'%(RTstructures[int(var1)]['name'], RTstructures[int(var2)]['name'])

	ovhc=calculate_ovh(rtssfile, int(var1), int(var2), True)
if __name__ == "__main__":
    main()
