import numpy as np
import math
#def find_point_on_sphere(x,y,r,center):
#    a= float(center[0])
#    b= float(center[1])
#    c= float(center[2])
#    x=float(x)
#    y=float(y)
#    r=float(r)
#    z= (r**2-(x-a)**2-(y-b)**2)**(0.5)+c
#    print(z)
#    return [round(x,5),round(y,5),round(z,5)]
#
#def step_on_sphere(start_point,r,center,step):
#    x=float(start_point[0])-float(step)
#    y=float(start_point[1])+float(step)
#
#    return find_point_on_sphere(x,y,r,center)
def step_on_sphere(start_point,center,step1,step2):
    
    x=float(start_point[0])
    y=float(start_point[1])
    z=float(start_point[2])
    
    a= float(center[0])
    b= float(center[1])
    c= float(center[2])

    r=np.sqrt((x-a)**2+(y-b)**2+(z-c)**2)

    theta=np.arccos((z-c)/r)
    #print(theta)
    theta-=step1
    phi=np.arctan((y-b)/(x-a))
    phi-=step2
    #print(phi)
    #phi+=step
    #print(start_point)
    #print(center)
    #print(theta)
    #print(phi)
    x1=a+r*np.sin(theta)*np.cos(phi)
    y1=b+r*np.sin(theta)*np.sin(phi)
    z1=c+r*np.cos(theta)

    return [x1,y1,z1]



