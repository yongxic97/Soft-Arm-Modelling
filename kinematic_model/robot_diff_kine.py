# First written for Graduation Project 'Research on dynamic modelling of a soft robot arm with peumatic bellows actuators and a
# damping support structure'. Functions in this module deals with differential kinematics of the robot.
#
# Author: Yongxi Cao  (Affiliation when doing the project: Dept. Mechincal & Energy Engineering, SUSTech)
# Email: yongxi.cao.2000@gmail.com
# v1.0: 2022-Mar-12
#
# function list:
#   1. Gamma(aster)
#   2. dHdqLocal(H,q,i)
#   3. dHdqWorld(H_base,H_top,q,way,section)
#   4. dHdq2ndLocal()
#   3. dGammad(aster,dWhat)
import numpy as np
import numpy.matlib
import robot_trans
from configs.numerical_settings import *
from configs.dimension_params import *
import copy

def Gamma(aster):

    """kinematic_model\\robot_diff_kine.py
    A linear operator used for convinently calculating ratation rates.
    
    Input: Some 3*3 matrix.
    
    Output: Corresponding rasult after gamma operator.
    """

    preMut = np.matlib.zeros((3,3))
    postMut = np.matlib.zeros((3,3))
    preMut[1,0],preMut[0,2],preMut[2,1] = 1,1,1
    postMut[1,0],postMut[0,2],postMut[2,1] = 1,1,1

    mutRes = preMut@aster@postMut
    diag = np.diagonal(mutRes)

    return diag


def dHdqLocal(q,way):
    """kinematic model\\robot_diff_kine.py
    First order derivative of H transformation of a single section. The meaning of this function: When the length of 
    a certain muscle changes, the H matrix's rate of changing is indicated by the dH matrix.
    
    Presumption: this muscle is inside the section of the H matrix, not in other sections.

    Input:  q: current length of three muscles, a 1*3 veecor.
            way: which way of muscle to move, i from {1,2,3}, not {0,1,2}!
           
    Output: the partial derivative matrix dHdq.
    """
    # 暂时假设三者效应不能线性分割，否则奇异位置计算会遇到问题，且暂时无法解决。
    # 按phi kappa length三者效应叠加处理。
    stepLength = INTERSECTION1 # 0.001mm

    # 这里是为数值计算考量，除以的倍数10000放大和缩小量级貌似影响很小
    # if q[0] == q[1] and q[0] == q[2]:# 为避免奇异性，在三路等长时给一个待改变方向的微小弧长，长度相对于数值计算的步长也足够小
    #     q[way-1] += INTERSECTION1 / 100000
    # 经过简单测试，当输入的三路肌肉长度显著不一样长时（至少其中两路差大于1mm）以上的微调对保留的结果没有影响
    # 但是这种处理在奇异位置的有效性，有待进一步确认。如果以后发现计算有问题，再回来检查这里。
    H1 = robot_trans.actuator2homoTrans(q)
    k1 = robot_trans.actuator2config(q)
    q_new = q.copy()
    q_new[way-1] = q[way-1] + stepLength
    # calculate new k and H caused by change of three muscles
    H2 = robot_trans.actuator2homoTrans(q_new)
    k2 = robot_trans.actuator2config(q_new)

    # k_temp1 = np.array([k2[0],k1[1],k1[2]]) # new k by change of way1
    # k_temp2 = np.array([k1[0],k2[1],k1[2]]) # new k by change of way2
    # k_temp3 = np.array([k1[0],k1[1],k2[2]]) # new k by change of way3

    # H2_1 = robot_trans.config2homoTransForDiffCal(k_temp1) # new k assuming change of only phi
    # H2_2 = robot_trans.config2homoTransForDiffCal(k_temp2) # new k assuming change of only kappa
    # H2_3 = robot_trans.config2homoTransForDiffCal(k_temp3) # new k assuming change of only length
    # diffk = k2 - k1
    # the fourth row are just all zeros.
    H2 = robot_trans.actuator2homoTrans(q_new)
    # dH = (H2_1+H2_2+H2_3-3*H1) / stepLength
    dH = (H2-H1) / stepLength
            
    return dH


def dHdqWorld(H_base,H_top,q,way,section):
    """kinematic model\\robot_diff_kine.py
    Derivative of H transformation in the world frame sense. The function takes the two head's H matrices directly.
    
    Presumption: this muscle is inside the section of the H matrix, not in other sections.

    Input:  H_base and H_top: the H matrices from bottom to the base of this section
            q: current length of three muscles, a 1*3 veecor.
            way: which way of muscle to move, i from {1,2,3}, not {0,1,2}!
            section: which section is the change happening in (not quite useful for now, but will indicate something).
    Output: the partial derivative matrix dHdq.
    """
    sectional_dH = dHdqLocal(q, way)

    dH = H_base@sectional_dH@H_top

    return dH


def dHdq2ndLocal(q,way1,way2):
    """kinematic model\\robot_diff_kine.py
    Second order derivative of H transformation of a single section. The meaning of this function: When the length of 
    two certain muscles change, the H matrix's second order rate of changing is indicated by the ddH matrix.
    
    Presumption: this muscle is inside the section of the H matrix, not in other sections.

    Input:  q: current length of three muscles, a 1*3 veecor.
            way1,way2: which ways of muscle to move, i from {1,2,3}, not {0,1,2}!
           
    Output: the partial derivative matrix dHdq.
    """

    stepLength = INTERSECTION1 # 0.001mm

    q_new1 = q.copy()
    q_new2 = q.copy()
    q_new1[way1-1] = q[way1-1] + float(stepLength)
    q_new2[way2-1] = q[way2-1] + float(stepLength)

    # calculate k1 and H1
    H1 = robot_trans.actuator2homoTrans(q)
    k1 = robot_trans.actuator2config(q)
   
    return ddH


def dHdt(q_list,dq_list):

    """kinematic model\\robot_diff_kine.py
    First order time derivative of world frame H matrix.

    Input: q_list: all muscle lengths. Length of the vector: 3*numOfSection
           dq_list: all muscle velocities. Length of the vector also 3*numOfSection. Unit: meter/second
           (abolish) numOfSection: How many sections are involved. Now suppose fixed at 4 sections.

    Output: dH: time derivative of world frame H matrix.
    """

    H = np.zeros((4,4,4))

    for i in range(4):
        H[:,:,i] = robot_trans.actuator2homoTrans(q[(i*3):(i*3+2)]) # calculate local H matrices.

    dH = np.zeros((4,4,4,3))
    dHdt = np.zeros((4,4))
    dHworld = np.zeros((4,4,4,3))

    for i in range(4):
        for j in range(3):
            
            dH[:,:,i,j] = dHdqLocal(q[(3*i):(3*i+2)], j+1)         # calculate local dH matrices.
            frontH = np.identity(4)
            endH = np.identity(4)

            for s in range(i):
                frontH = frontH @ H[:,:,s]
            for s in range(i+1,4):
                endH = endH @ H[:,:,s]

            dHworld[:,:,i,j] = dHdqWorld(frontH, endH, q[(3*i):(3*i+2)], j+1, i)
            dHdt += dHworld[:,:,i,j] * dq_list[3*i+j]

    return dHdt


def drdq(r,q,section,way):
    
    """kinematic_model\\robot_diff_kine.py
    How a position changes with the change of one muscle.
    Note that r must be extended: [x y z 1] as a 1*4 vector.
    Input:  r: local position in i* frame.
            q: state (length) of all muscles. (typically here its length: 3*4 = 12, from bottom to head) 
            section: the section the change of muscle length happens in. 
            way: the way of muscle that changes.
            # H_base: the H trans from world frame to the local frame that r lies in.
    Output: derivative of r relative to q_ik in world frame.
    """
    H_base = np.identity(4)
    H_up = np.identity(4)

    for i in range(section):
        H_base = H_base @ robot_trans.actuator2homoTrans(q[(3*i):(3*i+2)])
    
    for i in range(section+1,4):
        H_up =  H_up @ robot_trans.actuator2homoTrans(q[(3*i):(3*i+2)])

    dH = dHdqWorld(H_base, H_up, q[(3*section):(3*section+2)], way, section)
    dr = dH @ r # 最简单的case，计算某个移动坐标系的原点的位置，这里r就是[0 0 0 1]的形式
                # 对应的dH就是从底座到这个移动坐标系的H矩阵对下面某个肌肉变化的求导

    return dr


def drdt(r_local,q_list,dq_list):

    """
    Time derivative of a point.
    r_local: position vector, need to be extend to [x y z 1]
    """

    dH = dHdt(q_list, dq_list)
    dr = dH @ r_local

    return dr


def dGammad(aster,dWhat):
    return ddiag

if __name__ == '__main__':

    # Following are tests to the functions above.
    
    # diffk = np.array([0.001,0.002,0.003])
    # H1 = np.matlib.zeros((4,4))
    # H2 = np.matlib.zeros((4,4))
    # H1[0,0] = 5
    # H2[0,0] = 5.001

    # res1 = np.sum((H2[0,0]-H1[0,0])/diffk)
    # res2 = (H2[0,0]-H1[0,0])/diffk[0] + (H2[0,0]-H1[0,0])/diffk[1] + (H2[0,0]-H1[0,0])/diffk[2]
    # print(res1,res2) # this way of using numpy.sum is correct.
    dh = dHdqLocal(np.array([50e-3,50e-3,50e-3],dtype=np.float64),way=2) # 输入以米为单位,64位精度很重要，否则数值上误差会很大！
    print(dh)

    # H = np.ones((4,4))
    # H = H * 3
    # print(H)

    H = np.zeros((4,4))
    H[2,:] = 4,4,4,4
    h = np.array([1,2,3,4])
    print(H,h)
    H = H@h
    print(H)