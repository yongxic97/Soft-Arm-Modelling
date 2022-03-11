# First written for Graduation Project 'Dynamic modelling of a soft robot arm with peumatic bellows and a
# damping support structure'. Functions in this file deals with robot-specefic transformations of the robot.
# Arthor: Yongxi Cao 
# Email: yongxi.cao.2000@gmail.com
# v1.0: 2022-Mar-12
# function list:
#   1. actuator2config()
import numpy as np

def actuator2config(q,r_b):
    """kinetic model\\robot_specific_trans.py

        Performs robot-specific transformation to constant curvature continuum robot arms with 3 ways of muscles.

        Adapted from FESTO paper [2015 T-RO], can be applied to our prototype ESNAKE.

        Input: q = [q1 q2 q3] is the length vector. representing the length of 3 muscles of one section. Unit: mm.
               r_b is the distance of the robot center and one bellow's center. Unit: mm.
               
        Output: k = [phi kappa length] is the constant curvature configuration set. Unit: mm.
    """
    q1,q2,q3 = q[0],q[1],q[2]

    phi   = np.arctan2(np.sqrt(3)*(q3-q2),q2+q3-2*q1)
    kappa = 2 * sqrt(q1**2+q2**2+q3**2-q1*q2-q2*q3-q3*q1) \
            / ((q1+q2+q3)*r_b)
    leng  = (q1+q2+q3) / 3

    k = np.array([phi,kappa,leng])

    return k
    
def config2homoTrans(k):
    """kinetic model\\robot_specific_trans.py

        Transforms constant curvature configuration set to corresponding homogeneous transformation matrix. 

        This is the well-recognized way to perform transformation. (5 D-H transforms)

        Input: k = [phi kappa length] is the constant curvature configuration set. Unit: mm.

    """

    return H

if __name__ == '__main__':
    actuator2config(q)