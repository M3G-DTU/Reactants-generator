from molrearr.orient_mod import *

def orient_fakefile(fakefile,args): #call_orient_from_fakefile
    # if a real file or a fakefile
    arg_for_final='intermediate.xyz '+str(args)
    #args such as: 'intermediate.xyz -rx 90'
    final=call_orient2(fakefile,arg_for_final.split())

    return final

def orient_file(filename,args): #call_orient_from_file
    # if a real file or a fakefile
    if len(filename.split()) == 1:
        with open(filename,"r") as file:
            fakefile=file.read()
        #orient_script1=('python orient_mod.py '+str(filename)+' '+str(args))
        # such as: 'python orient.py benzo.xyz -rz 90'
        intermediate=orient_fakefile(fakefile,args)
    else:
        intermediate=orient_fakefile(filename,args)
    return intermediate



if __name__ == "__main__":
    dum1=orient_file('benzo.xyz','-rz 90')
    print(dum1)
    #with open('new1.xyz',"w") as inter_file:
    #    inter_file.write(dum1)
    dum2=orient_file(dum1,'-rx 90')
    print(dum2)
    #with open('new2.xyz',"w") as final_file:
    #    final_file.write(dum2)
