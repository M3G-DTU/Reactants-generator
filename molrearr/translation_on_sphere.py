import numpy as np
from molrearr.find_point_on_sphere import *
from molrearr.fast_call_orient import *
from molrearr.read_xyz2list import *

def translation_on_sphere(mol,start,center,step_theta,step_phi):
    target1=step_on_sphere(start,center,math.pi/180*step_theta,math.pi/180*step_phi)
    vector1=np.subtract(target1,start)
    
    a=float(vector1[0])
    b=float(vector1[1])
    c=float(vector1[2])
    
    args='-tx '+str(a)+' -ty '+str(b)+' -tz '+str(c)

    target2=orient_file(mol,args)

    return target2

if __name__ == "__main__":
    mol_1=read_xyz2list('benzo.xyz')
    hotspot1=hot_coordinate(mol_1,2)
    new=translation_on_sphere('benzo.xyz',hotspot1,[0.5,0.5,0.5],60,60)

    print(mol_1)
    print(hotspot1)
    print(new)
    

    with open('new1.xyz',"w") as inter_file:
        inter_file.write(new)
