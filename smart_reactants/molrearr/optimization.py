import numpy as np
import geatpy as ea
from molrearr.read_xyz2list import *
from molrearr.drone_move import *
from molrearr.quality_check import *

def angle2dist(mol1: str, mol2: str, hot1: int, hot2: int, anglelist: list) -> float:
    """Calculate the distance between the hotspot of molecule 1 and the hotspot of molecule 2 after applying the specified angles to molecule 2."""
    molecule1=read_xyz2list(mol1)

    new2=drone_move(mol1,hot1,mol2,hot2,anglelist)
    molecule_2finallist=read_xyz2list(new2)
    
    min_dist=other_atom_in_mol2_to_hot1(molecule1,molecule_2finallist,hot2)

    return min_dist

def block_due2angle(mol1: str, mol2: str, hot1: int, hot2: int, anglelist: list, strict_filter_level: int) -> int:
    """Calculate the number of blocks due to the specified angles between the hotspots of two molecules."""
    new2=drone_move(mol1,hot1,mol2,hot2,anglelist)
    molecule1=read_xyz2list(mol1)
    molecule2=read_fakexyz2list(new2)
    block_number=point_between_hotspots(molecule1,molecule2,hot1,hot2,strict_filter_level)
    return int(block_number[-1])
    



def geat_optimization(mol1: str, mol2: str, hot1: int, hot2: int, minimal_distance: float, strict_filter_level: int) -> list:
    """Perform optimization using the GEAT algorithm to find the optimal angles for molecule 2 to minimize the distance between the hotspots of two molecules while considering blocking constraints."""
    # 构建问题 -> Build issues
    #r = 1  # 目标函数需要用到的额外数据 -> Extra data needed for the objective function
    @ea.Problem.single
    def evalVars(Vars):  # 定义目标函数（含约束）-> Define objective function (with constraints)
        #f = np.sum((Vars - r) ** 2)  # 计算目标函数值 -> Calculate the objective function value
        x1 = Vars[0]
        x2 = Vars[1]
        x3 = Vars[2]
        x4 = Vars[3]
        x5 = Vars[4]
        anglelist=[x1,x2,x3,x4,x5]
        f= angle2dist(mol1,mol2,hot1,hot2,anglelist)
        CV = np.array([minimal_distance-angle2dist(mol1,mol2,hot1,hot2,anglelist),
                        block_due2angle(mol1,mol2,hot1,hot2,anglelist,strict_filter_level)-1])  # 计算违反约束程度 -> Calculate the degree of constraint violation

        return f, CV

    problem = ea.Problem(name='rotate the right angle',
                            M=1,  # 目标维数 -> Objective dimension
                            maxormins=[-1],  # 目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标 -> List of objective minimization/maximization flags, 1: minimize this objective; -1: maximize this objective
                            Dim=5,  # 决策变量维数 -> Decision variable dimension
                            varTypes=[1,1,1,1,1],  # 决策变量的类型列表，0：实数；1：整数 -> List of decision variable types, 0: real; 1: integer
                            lb=[-60,-60,-60,-60,-60],  # 决策变量下界 -> Decision variable lower bounds
                            ub=[60,60,60,60,60],  # 决策变量上界 -> Decision variable upper bounds
                            evalVars=evalVars)
    # 构建算法 -> Build algorithm
    algorithm = ea.soea_SEGA_templet(problem,
                                        ea.Population(Encoding='RI', NIND=20),
                                        MAXGEN=50,  # 最大进化代数 -> Maximum number of generations
                                        logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录 -> Log interval, 0 means no logging
                                        trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。 -> Stagnation threshold for single-objective optimization.
                                        maxTrappedCount=10)  # 进化停滞计数器最大上限值。-> Maximum limit for the evolution stagnation counter.
    # 求解 -> Solve
    res = ea.optimize(algorithm, seed=1, verbose=False, drawing=0, outputMsg=False, drawLog=False, saveFlag=False, dirName='result')
    
    return res['Vars'][0]

def optimal_position(mol1: str, mol2: str, hot1: int, hot2: int, minimal_distance: float, strict_filter_level: int) -> str:
    """Find the optimal position for molecule 2 relative to molecule 1 using the GEAT optimization results."""
    anglelist = geat_optimization(mol1, mol2, hot1, hot2, minimal_distance, strict_filter_level)
    optimal = drone_move(mol1, hot1, mol2, hot2, anglelist)

    return optimal


if __name__ == '__main__':
    
    mol1='benzo.xyz'
    mol2='new1.xyz'
    hot1=2
    hot2=4
    minimal_distance=1.6
    strict_filter_level=3
    final= geat_optimization(mol1,mol2,hot1,hot2,minimal_distance, strict_filter_level)
    print(final)
    optimal=optimal_position(mol1,mol2,hot1,hot2,minimal_distance, strict_filter_level)
    print(optimal)