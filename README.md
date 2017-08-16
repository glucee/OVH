# OVH
Python to generate Overlap Volume Histogram from DICOM file

Filename: ovhcalc.py

Description: the file contains method calculating OVH graph from rtss DICOM file by identifying the ROI number of 2 organs. It also returns a dict with min_Dis and max_Dis for x-axis, and a list for OVH plot.

Main function: 
calculate_ovh(structure, ptv_ori, non_ptv_roi, wplt=False): 
Wplt: True: show OVH graph
                    
The steps in this function are the following:
Distinguish PTV and Organs, because the ovh is the geometric relationship between the OARs and PTV. We store all the PTV in the ptv_points[] and store all the organs in the non_ptv_points[] separately.
Calculate the shortest distance from each point in the Organs to the surface of the target PTV. Please note the Orangs and PTV are all 3-D objects, calculating on only 2-D surface may result in the wrong distance. This script calculate the real distance by considering the 3-D coordinate.
Extract the CDF of the distance array, which is the required OVH plot. The function has 3 iterated loops which could be optimized. 

Main(): 
This function prints all the PTVs and Organs, and outputs OVH graph for selected PTV and Organ.

Usage:

Python ovhcalc.py
Then all the PTVs and Organs will be printed with (ID: Name)

------PTV------
ID: Name

------Non-PTV------
ID: Name

The user could select a PTV ID and a organ ID to calculate, finally a OVH figure will be saved to the same folder.
