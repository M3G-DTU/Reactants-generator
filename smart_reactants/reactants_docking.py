import os
import math
import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict

# dependencies can be any iterable with strings, 
# e.g. file line-by-line iterator
dependencies = [
  'numpy',
  'geatpy',
]

# here, if a dependency is not met, a DistributionNotFound or VersionConflict
# exception is thrown. 
pkg_resources.require(dependencies)

#functions
def molecule_info(out_file):
    openold = open(out_file,"r")
    rline = openold.readlines()
    for i in range (len(rline)):
        if "CONDENSED LOCAL ELECTROPHILICITY AND NUCLEOPHILICITY" in rline[i]:
            start = i+5
    for m in range (start, len(rline)):
        if "Total :" in rline[m]:
            end = m-1
            break
    atom_properties=[]
    for line in rline[start:end]:
        ap=line.split()
        ap.remove(':')
        atom_properties.append(ap)

    return atom_properties

def fukui_info(out_file):
    openold = open(out_file,"r")
    rline = openold.readlines()
    for i in range (len(rline)):
        if "ATOMIC DESCRIPTORS: CANONICAL ENSEMBLE" in rline[i]:
            start = i+5
    for m in range (start, len(rline)):
        if "Total :" in rline[m]:
            end = m-1
            break
    atom_properties=[]
    for line in rline[start:end]:
        ap=line.split()
        ap.remove(':')
        atom_properties.append(ap)

    return atom_properties

# filter
def elec_filter(list,fukui_list):
    filtrate=[]
    for atom in sorted(list,key=(lambda x:float(x[2])),reverse=True):
        if float(atom[2])>0:
            if float(atom[2])+float(atom[3])>0:
                filtrate.append(atom)
            else:
                continue
        else:
            continue
    if len(filtrate) == 0:
        for atom in sorted(fukui_list,key=(lambda x:float(x[5])),reverse=True):
            if float(atom[-1]) > 0:
                filtrate.append([atom[0],atom[1],atom[-1],0])
    #print (filtrate)
    return filtrate

def nuc_filter(list,fukui_list):
    filtrate=[]
    for atom in sorted(list,key=(lambda x:float(x[3]))):
        if float(atom[3])<0:
            if float(atom[2])+float(atom[3])<=0:
                filtrate.append(atom)
            else:
                continue
        else:
            continue

    if len(filtrate) == 0:
        for atom in sorted(fukui_list,key=(lambda x:float(x[5]))):
            if float(atom[-1]) < 0:
                filtrate.append([atom[0],atom[1],0,atom[-1]])
    #print (filtrate)    
    return filtrate

def keep_top(list, number, H_indicator):
    top=[]
    get_number=float(number)
    if get_number == 0:
        return list
    elif get_number> len(list):
        print('* top number in setting > the number of candidates in this filter, return all the candidates')
        return list
    elif get_number>=1 and get_number<=len(list):
        true_number=int(get_number)
        top=list[:true_number]
        if H_indicator == 1:
            H_number=sum(row.count('H') for row in top)
            #print(H_number)
            H_list=[]
            for atom in list:
                if "H" in atom:
                    H_list.append(atom)
            top.extend(H_list[H_number:3])
        return top
    elif 0<get_number<1:
        presentage=math.ceil(get_number * len(list))
        top=list[:presentage]
        if H_indicator == 1:
            H_number=sum(row.count('H') for row in top)
            #print(H_number)
            H_list=[]
            for atom in list:
                if "H" in atom:
                    H_list.append(atom)
            top.extend(H_list[H_number:3])
        return top

yes_choices = ['yes', 'y']
no_choices = ['no', 'n']

while True:
    weight_H=input("Would you like to weight \"H\" more? Y or N:")
    if weight_H.lower() in yes_choices:
        print('user typed yes, top 3 "H"s in each filter will be considered')
        indicator_H=1
        break
    elif weight_H.lower() in no_choices:
        print('user typed no')
        indicator_H=0
        break
    else:
        print("please type your answer again, Y or N:")
        continue


# read setting file
setting_file=open('reactants_docking_setting.txt',"r")
sline=setting_file.readlines()
#print(sline)

# files
xyz_file1=str(sline[1].rstrip('\n'))
xyz_file2=str(sline[3].rstrip('\n'))
out_file1=str(sline[5].rstrip('\n'))
out_file2=str(sline[7].rstrip('\n'))

# set critera
elec_cutoff=float(sline[10])
nucl_cutoff=float(sline[12])

# output directory
directory = str(sline[14].rstrip('\n'))
isExist = os.path.exists(directory)
if not isExist:

   # Create a new directory because it does not exist
   os.makedirs(directory)

# get information from .out files
info_out1=molecule_info(out_file1)
info_out2=molecule_info(out_file2)
fukui_out1=fukui_info(out_file1)
fukui_out2=fukui_info(out_file2)

# sort electro or nucleo
# print each sort results
print('---Candidates---')
# 1
print(str(out_file1.split('/')[-1].split('.')[0])+'_eletro_high2low: ')
out1_eletro_high2low=keep_top(elec_filter(info_out1,fukui_out1),sline[18],indicator_H)
print([i[0]+i[1] for i in out1_eletro_high2low])
# 2
print(str(out_file1.split('/')[-1].split('.')[0])+'_nucleo_low2high: ')
out1_nucleo_low2high=keep_top(nuc_filter(info_out1,fukui_out1),sline[18],indicator_H)
print([i[0]+i[1] for i in out1_nucleo_low2high])
# 3
print(str(out_file2.split('/')[-1].split('.')[0])+'_eletro_high2low: ')
out2_eletro_high2low=keep_top(elec_filter(info_out2,fukui_out2),sline[18],indicator_H)
print([i[0]+i[1] for i in out2_eletro_high2low])
# 4
print(str(out_file2.split('/')[-1].split('.')[0])+'_nucleo_low2high: ')
out2_nucleo_low2high=keep_top(nuc_filter(info_out2,fukui_out2),sline[18],indicator_H)
print([i[0]+i[1] for i in out2_nucleo_low2high])

# about output files
print('---Acquire elec/nuc cutoffs---')
print(str(elec_cutoff)+' '+str(nucl_cutoff))
print('---Output files---')

for atom in out1_eletro_high2low:
    hot1=int(atom[0])
    if float(atom[2])>= float(out1_eletro_high2low[0][2])-elec_cutoff:
        for atom2 in out2_nucleo_low2high:
            hot2=int(atom2[0])
            if float(atom2[3]) <= float(out1_nucleo_low2high[0][3])+nucl_cutoff:
                produce_script=('python molecule_rearrangement.py '+str(xyz_file1)+' '+str(xyz_file2)+' '+str(hot1)+' '+str(hot2)+' '+str(directory)+'/'+ str(xyz_file1.split('/')[-1].split('.')[0])+'_'+str(xyz_file2.split('/')[-1].split('.')[0])+'_ElecNu_'+str(hot1)+str(atom[1])+'_'+str(hot2)+str(atom2[1])+'.xyz')
                print('-------------------')
                print(produce_script.split()[-1])
                produce_file=os.system(produce_script)
            else:
                continue
    else:
        continue

print('---Exchange---')

for atom in out2_eletro_high2low:
    hot1=int(atom[0])
    if float(atom[2])>= float(out2_eletro_high2low[0][2])-elec_cutoff:
        for atom2 in out1_nucleo_low2high:
            hot2=int(atom2[0])
            if float(atom2[3]) <= float(out1_nucleo_low2high[0][3])+nucl_cutoff:
                produce_script=('python molecule_rearrangement.py '+str(xyz_file2)+' '+str(xyz_file1)+' '+str(hot1)+' '+str(hot2)+' '+str(directory)+'/'+str(xyz_file2.split('/')[-1].split('.')[0])+'_'+str(xyz_file1.split('/')[-1].split('.')[0])+'_ElecNu_'+str(hot1)+str(atom[1])+'_'+str(hot2)+str(atom2[1])+'.xyz')
                print('-------------------')
                print(produce_script.split()[-1])
                produce_file=os.system(produce_script)
            else:
                continue
    else:
        continue
