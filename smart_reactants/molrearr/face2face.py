import numpy as np
import math
from molrearr.fast_call_orient import *
from molrearr.read_xyz2list import *

def rotation_matrix(vec1,vec2):
    a,b=(np.array(vec1)/np.linalg.norm(vec1)).reshape(3), (np.array(vec2)/np.linalg.norm(vec2)).reshape(3)
    v=np.cross(a,b)
    if any(v):  # if not all zeros then
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
 
    else:
        return np.eye(3)  # cross of all zeros only occurs on identical directions

# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6
 
# Calculates rotation matrix to euler angles
# The result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R) :
 
    assert(isRotationMatrix(R))
     
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    #sy = math.sqrt(R[2,1] * R[2,1] + R[2,2] * R[2,2])
     
    singular = sy < 1e-6
 
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
 
    return np.array([math.degrees(x), math.degrees(y), math.degrees(z)])

def rotationMatrixToAxleAngel(R):
    assert(isRotationMatrix(R))

    theta=math.acos((R[0,0]+R[1,1]+R[2,2]-1)*0.5)
    vx=1/(2*math.sin(theta))*(R[2,1]-R[1,2])
    vy=1/(2*math.sin(theta))*(R[0,2]-R[2,0])
    vz=1/(2*math.sin(theta))*(R[1,0]-R[0,1])

    return [math.degrees(theta),vx,vy,vz]

def face2face(mol2,hotspot_in_file2,vector1,vector2,target1):

    rota1=rotation_matrix(vector2,vector1*(-1))
    axleangle1=rotationMatrixToAxleAngel(rota1)
    
    args=' -rv '+str(axleangle1[0])+' '+str(axleangle1[1])+' '+str(axleangle1[2])+' '+str(axleangle1[3])
    
    dum1=orient_file(mol2,args)
    
    
    dum_molecule=read_fakexyz2list(dum1)

    hot2_dum=hot_coordinate(dum_molecule,hotspot_in_file2)
    vectors=np.subtract(target1,hot2_dum)

    args2='-tx '+str(vectors[0])+' -ty '+str(vectors[1])+' -tz '+str(vectors[2])
    final=orient_fakefile(dum1,args2)

    return final

def face2faceS(mol1,mol2,hot1,hot2,restrict):
    mol_1=read_xyz2list(mol1)
    mol_2=read_xyz2list(mol2)
    hotspot1=hot_coordinate(mol_1,hot1)
    hotspot2=hot_coordinate(mol_2,hot2)
    center1=center_of_geometry(mol_1)
    center2=center_of_geometry(mol_2)
    hot_center_vector1=np.subtract(hotspot1, center1)
    hot_center_vector2=np.subtract(hotspot2,center2)

    hot_center_distance1=two_point_distance(hotspot1, center1)
    
    restrict_distance=float(restrict)

    ratio1=restrict_distance/hot_center_distance1
    target1=[hotspot1[0]+hot_center_vector1[0]*ratio1,hotspot1[1]+hot_center_vector1[1]*ratio1,hotspot1[2]+hot_center_vector1[2]*ratio1]

    new=face2face(mol2,hot2,hot_center_vector1,hot_center_vector2,target1)
    return new


if __name__ == "__main__":

    #mol_1=read_xyz2list('benzo.xyz')
    #mol_2=read_xyz2list('benzo.xyz')
    #hotspot1=hot_coordinate(mol_1,2)
    #hotspot2=hot_coordinate(mol_1,4)
    #center1=center_of_geometry(mol_1)
    #center2=center_of_geometry(mol_2)
    #hot_center_vector1=np.subtract(hotspot1, center1)
    #hot_center_vector2=np.subtract(hotspot2,center2)
    #hot_center_distance1=two_point_distance(hotspot1, center1)
    #restrict_distance=2.2
    #ratio1=restrict_distance/hot_center_distance1
    #target1=[hotspot1[0]+hot_center_vector1[0]*ratio1,hotspot1[1]+hot_center_vector1[1]*ratio1,hotspot1[2]+hot_center_vector1[2]*ratio1]
    #new=face2face('benzo.xyz',4,hot_center_vector1,hot_center_vector2,target1)
    #print(mol_1)
    #print(hotspot1)
    new=face2faceS('benzo.xyz','benzo.xyz',2,4,2.2)
    print(new)
    

    with open('new1.xyz',"w") as inter_file:
        inter_file.write(new)