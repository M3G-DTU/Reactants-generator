import io
import numpy as np

def read_fakexyz2list(fakefile):
    out = []
    f = io.StringIO(fakefile)
    
    line = f.readline()
    while line != "":
        natoms = int(line)
        comment = f.readline().rstrip()
        names = []
        coords = []
        #extras = []
        for i in range(natoms):
            line = f.readline()
            data = line.split()
            name, x, y, z = data[0:4]
            #extra = data[4:]
            names.append(name.capitalize())
            coords.append([float(x), float(y), float(z)])
            #if extra:
            #extras.append(extra)
        
        out.append(names)
        out.append(coords)
        out.append(comment)
        out.append(natoms)
        line = f.readline()
            
    return out
def read_xyz2list(filename):
    if len(filename.split()) == 1:
        with open(filename,"r") as file:
            fakefile=file.read()
        return read_fakexyz2list(fakefile)
    else:
        return read_fakexyz2list(filename)


def hot_coordinate(molecule,index):
    return molecule[1][int(index)-1]

def center_of_geometry(molecule):
    n,x,y,z=(0,0,0,0)
    a=molecule
    for index, element in enumerate(a[0]):
        if element != 'W':
            n+=1
            x+=a[1][index][0]
            y+=a[1][index][1]
            z+=a[1][index][2]

        cog=[x/n,y/n,z/n]

    return cog

def two_point_distance(coordinate1, coordinate2):
    return np.sqrt(sum(np.power((np.array(coordinate1)-np.array(coordinate2)),2)))

if __name__ == '__main__':

    f=open('benzo.xyz',"r")
    #print(f)
    read=f.read()
    #print(read)
    f1=read_xyz2list('benzo.xyz')
    f2=read_xyz2list(read)
    print(f1)
    print(f2)

    #orient_script=('python Test/orient.py Test/benzo.xyz -tx 10 -ty 10 -tz 10')
    #xyz_new=os.system(orient_script)
    #print(xyz_new)
    #print(type(xyz_new))