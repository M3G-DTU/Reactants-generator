import numpy as np
import geatpy as ea
from molrearr.read_xyz2list import *
from molrearr.drone_move import *
from molrearr.quality_check import *

setting_file=open('reactants_docking_setting.txt',"r")
sline=setting_file.readlines()

strict_filter_level=float(sline[21])
other_dist=float(sline[23])

def angle2dist(mol1,mol2,hot1,hot2,anglelist):
    molecule1=read_xyz2list(mol1)

    new2=drone_move(mol1,hot1,mol2,hot2,anglelist)
    molecule_2finallist=read_xyz2list(new2)
    
    min_dist=other_atom_in_mol2_to_hot1(molecule1,molecule_2finallist,hot2)

    return min_dist

def block_due2anlge(mol1,mol2,hot1,hot2,anglelist,factor):
    new2=drone_move(mol1,hot1,mol2,hot2,anglelist)
    molecule1=read_xyz2list(mol1)
    molecule2=read_fakexyz2list(new2)
    block_number=point_between_hotspots(molecule1,molecule2,hot1,hot2,factor)
    return int(block_number[-1])
    



def geat_optimization(mol1,mol2,hot1,hot2):

    # 构建问题
    #r = 1  # 目标函数需要用到的额外数据
    @ea.Problem.single
    def evalVars(Vars):  # 定义目标函数（含约束）
        #f = np.sum((Vars - r) ** 2)  # 计算目标函数值
        x1 = Vars[0]
        x2 = Vars[1]
        x3 = Vars[2]
        x4 = Vars[3]
        x5 = Vars[4]
        anglelist=[x1,x2,x3,x4,x5]
        f= angle2dist(mol1,mol2,hot1,hot2,anglelist)
        CV = np.array([other_dist-angle2dist(mol1,mol2,hot1,hot2,anglelist),
                        block_due2anlge(mol1,mol2,hot1,hot2,anglelist,strict_filter_level)-1])  # 计算违反约束程度

        return f, CV

    problem = ea.Problem(name='rotate the right angle',
                            M=1,  # 目标维数
                            maxormins=[-1],  # 目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标
                            Dim=5,  # 决策变量维数
                            varTypes=[1,1,1,1,1],  # 决策变量的类型列表，0：实数；1：整数
                            lb=[-60,-60,-60,-60,-60],  # 决策变量下界
                            ub=[60,60,60,60,60],  # 决策变量上界
                            evalVars=evalVars)
    # 构建算法
    algorithm = ea.soea_SEGA_templet(problem,
                                        ea.Population(Encoding='RI', NIND=20),
                                        MAXGEN=50,  # 最大进化代数。
                                        logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                        trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                        maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, seed=1, verbose=False, drawing=0, outputMsg=False, drawLog=False, saveFlag=False, dirName='result')
    
    return res['Vars'][0]

def optimal_position(mol1,mol2,hot1,hot2):
    anglelist=geat_optimization(mol1,mol2,hot1,hot2)
    optimal=drone_move(mol1,hot1,mol2,hot2,anglelist)

    return optimal


if __name__ == '__main__':
    
    mol1='benzo.xyz'
    mol2='new1.xyz'
    hot1=2
    hot2=4
    final= geat_optimization(mol1,mol2,hot1,hot2)
    print(final)
    optimal=optimal_position(mol1,mol2,hot1,hot2)
    print(optimal)