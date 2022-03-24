from robot_trans import *
from configs.dimension_params import *
import numpy as np
import math

def M_matrix_reduced(q_list):
    
    """mass_matrix.py
    Reduced M matrix, omitting revolution. Only the first term in the formula.
    All terms can be determined from q_list
    """

    M = np.zeros((12,12))

    m = np.array([M1,M2,M3,M4])

    for i in range(12):
        for j in range(12):
            s = math.floor(np.max((i,j)) / 3) # which section
            si = math.floor(i/3) # the section number of si
            sj = math.floor(j/3) # the sectino number of sj
            for ii in range(s,4):
                if sj > s or si > s:
                    pass # because at least one of the derivative will be 0
                else:
                    dHdi = 
                    dHdj
                    M[i,j] += m[s] * 

    return M

if __name__ == '__main__':

    q = np.array([50e-3,60e-3,40e-3, \
                53e-3,50e-3,55e-3, \
                49e-3,52e-3,50e-3, \
                47e-3,50e-3,58e-3],dtype=np.float64)

    # print(math.floor(np.max((0,4))/3))
    M = M_matrix_reduced(q)
    print(M)


