from molrearr.read_xyz2list import *
from molrearr.two_vectors_angle import *

def point_between_hotspots(mol1,mol2,hot_1,hot_2,factor):
    record=[]
    n=0
    m=0
    hot1_coord=hot_coordinate(mol1,hot_1)
    hot2_coord=hot_coordinate(mol2,hot_2)
    dist_hot=two_point_distance(hot1_coord,hot2_coord)
    hot1_hot2_vector=np.subtract(hot2_coord, hot1_coord)
    for index,atom in enumerate(mol1[1]):
        if two_point_distance(atom,hot2_coord)<dist_hot:
            record.append([1,index+1])
            n+=1
            block_hot_vector1=np.subtract(atom, hot1_coord)
            block_vector_angle=angle_between(block_hot_vector1,hot1_hot2_vector)
            if block_vector_angle <10*factor or block_vector_angle>(180-10*factor):
                m+=1

    for index,atom in enumerate(mol2[1]):
        if two_point_distance(atom,hot1_coord)<dist_hot:
            record.append([2,index+1])
            n+=1
            block_hot_vector2=np.subtract(atom, hot2_coord)
            block_vector_angle=angle_between(block_hot_vector2,hot1_hot2_vector)
            if block_vector_angle <10*factor or block_vector_angle>(180-10*factor):
                m+=1

    record.append(n)
    record.append(m)
    return record

def far_atom(mol, hot_coord):
    dist_list=[]
    for atom in mol[1]:
        dist=two_point_distance(atom,hot_coord)
        dist_list.append(float(dist))

    return dist_list.index(max(dist_list))

def other_near_atoms(mol1,mol2,hot_1,hot_2,factor):
    n=0
    for index,atom in enumerate(mol1[1]):
        if index != int(hot_1)-1:
            for index2,atom2 in enumerate(mol2[1]):
                if index2 != int(hot_2)-1:
                    dist=two_point_distance(atom,atom2)
                    if dist <= float(factor):
                        #print(index)
                        #print(index2)
                        n+=1
    
    return n

def other_atom_in_mol2_to_hot1(mol1,mol2,hot_2):
    dist_list=[]
    index2_list=[]
    for atom in mol1[1]:
        for index2,atom2 in enumerate(mol2[1]):
            if index2 != int(hot_2)-1:
                dist=two_point_distance(atom,atom2)
                dist_list.append(dist)
                index2_list.append(index2)
    
    min_dist=min(dist_list)
    #return [min_dist,index2[dist_list.index(min_dist)]]
    return float(min_dist)
    