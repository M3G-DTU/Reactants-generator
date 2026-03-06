from molrearr.center_rotate import *
from molrearr.translation_on_sphere import *
from molrearr.two_vectors_angle import *
from molrearr.read_xyz2list import *

# 5anlges for mol2's movement and rotation
def drone_move(mol1,hotspot_in_file1,mol2,hotspot_in_file2,anglelist):
    step_theta, step_phi, step_x, step_y, step_z = [float(x) for x in anglelist]
    molecule_1 = read_xyz2list(mol1)                                                # TODO change this to the list format
    hot1_coord = hot_coordinate(molecule_1,int(hotspot_in_file1))

    molecule_2=read_xyz2list(mol2)                                                  # TODO change this to the list format     
    hot_new=hot_coordinate(molecule_2,int(hotspot_in_file2))
    dum_molecule=translation_on_sphere(mol2,hot_new,hot1_coord,step_theta,step_phi)
    dum_list=read_xyz2list(dum_molecule)
    hot_new2_coord=hot_coordinate(dum_list,int(hotspot_in_file2))
    dum_molecule2=center_rotate(dum_molecule,hot_new2_coord,[float(step_x),float(step_y),float(step_z)])
    
    return dum_molecule2             
