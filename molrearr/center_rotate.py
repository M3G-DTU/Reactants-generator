from molrearr.read_xyz2list import *
from molrearr.fast_call_orient import *

def new_cartesian(mol, origin_point):
    x=float(origin_point[0])*(-1)
    y=float(origin_point[1])*(-1)
    z=float(origin_point[2])*(-1)
    args='-tx '+str(x)+' -ty '+str(y)+' -tz '+str(z)
    
    target=orient_file(mol,args)
    
    return target
def center_rotate(mol,origin_point,anglelist):
    origin=new_cartesian(mol, origin_point)
    a=float(anglelist[0])
    b=float(anglelist[1])
    c=float(anglelist[2])

    args2='-rx '+str(a)+' -ry '+str(b)+' -rz '+str(c)
    target=orient_fakefile(origin,args2)

    x=float(origin_point[0])
    y=float(origin_point[1])
    z=float(origin_point[2])

    args3='-tx '+str(x)+' -ty '+str(y)+' -tz '+str(z)
    back_to_mol=orient_fakefile(target,args3)


    return back_to_mol

if __name__ == "__main__":
    mol_1=read_xyz2list('benzo.xyz')
    hotspot1=hot_coordinate(mol_1,2)
    new=new_cartesian('benzo.xyz',hotspot1)
    new2=center_rotate('benzo.xyz',hotspot1,[10,10,10])
    print(mol_1)
    print(hotspot1)
    print(new)
    print(new2)

    with open('new1.xyz',"w") as inter_file:
        inter_file.write(new)
    with open('new2.xyz',"w") as inter_file:
        inter_file.write(new2)