import sys
import numpy as np
from molrearr.face2face import *
from molrearr.drone_move import *
from molrearr.quality_check import *
from molrearr.optimization import *
from molrearr.read_xyz2list import *

def merge_xyz(filename, mol_1,mol_2):
    with open(filename,"w") as merge_file:
        merge_file.write(str(int(mol_1[-1])+int(mol_2[-1]))+'\n')
        merge_file.write('\n')
        #print(mol_1[-1])

        for index in range(int(mol_1[-1])):
            merge_file.write("{:<3}{:17.9f}{:17.9f}{:17.9f}\n".format(mol_1[0][index],float(mol_1[1][index][0]),float(mol_1[1][index][1]),float(mol_1[1][index][2])))
        for index in range(int(mol_2[-1])):
            merge_file.write("{:<3}{:17.9f}{:17.9f}{:17.9f}\n".format(mol_2[0][index],float(mol_2[1][index][0]),float(mol_2[1][index][1]),float(mol_2[1][index][2])))

#molecule 1, immobile
xyz_file1=sys.argv[1]
molecule_1=read_xyz2list(xyz_file1)
#molecule 2, dock to molecule 1
xyz_file2=sys.argv[2]
molecule_2=read_xyz2list(xyz_file2)

hotspot_in_file1=sys.argv[3]
hotspot_in_file2=sys.argv[4]

# read setting file
setting_file=open('reactants_docking_setting.txt',"r")
sline=setting_file.readlines()

# hotspots distances setting
restrict_distance=float(sline[16])

# basic parameters of mol1 and mol2 in original position
center1=center_of_geometry(molecule_1)
center2=center_of_geometry(molecule_2)

hot1=hot_coordinate(molecule_1,int(hotspot_in_file1))
hot2=hot_coordinate(molecule_2,int(hotspot_in_file2))

hot_center_vector1=np.subtract(hot1, center1)
hot_center_vector2=np.subtract(hot2,center2)

hot_center_distance1=two_point_distance(hot1, center1)
hot_center_distance2=two_point_distance(hot2, center2)

# move to face to face first

position1=face2faceS(xyz_file1,xyz_file2,int(hotspot_in_file1),int(hotspot_in_file2),restrict_distance)

# basic parameter of mol2 in position1 (face to face)
molecule_2new=read_xyz2list(position1)


# optimization

molecule_2final=optimal_position(xyz_file1,position1,hotspot_in_file1,hotspot_in_file2)
molecule_2finallist=read_xyz2list(molecule_2final)


hot_new=hot_coordinate(molecule_2finallist,int(hotspot_in_file2))
center_new=center_of_geometry(molecule_2finallist)
hot1_to_hot_new=two_point_distance(hot1,hot_new)
center1_to_center_new=two_point_distance(center1,center_new)


#final quality check

#block atoms
strict_filter_level=int(sline[21])
find_block=point_between_hotspots(molecule_1,molecule_2finallist,int(hotspot_in_file1),int(hotspot_in_file2),strict_filter_level)
print('block_number:'+str(find_block[-2]))
print('vec_angle<'+str(strict_filter_level*10)+':'+str(find_block[-1])) 

#other atoms within 1.6 ang
other_bonds=other_near_atoms(molecule_1,molecule_2finallist,int(hotspot_in_file1),int(hotspot_in_file2),sline[23])
print('unused_bonds:'+str(other_bonds))

if hot1_to_hot_new <= center1_to_center_new:
    print('Good start-point!')
else:
    print('Need to check! ')
# merge xyz

merge_process=merge_xyz(str(sys.argv[5]),molecule_1,molecule_2finallist)
