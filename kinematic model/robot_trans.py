# First written for Graduation Project 'Research on dynamic modelling of a soft robot arm with peumatic bellows actuators and a
# damping support structure'. Functions in this module deals with robot-specefic and robot-independent transformations of the robot.
#
# Author: Yongxi Cao  (Affiliation when doing the project: Dept. Mechincal & Energy Engineering, SUSTech)
# Email: yongxi.cao.2000@gmail.com
# v1.0: 2022-Mar-12
#
# function list:
#   1. actuator2config()
#   2. config2homoTrans()
#   3. config2rotMat()
#   4. config2tipPosition()
#   5. actuator2homoTrans()
#   6. actuator2rotMat()
#   7. actuator2tipPosition()
import numpy as np


def actuator2config(q,r_b):

    """kinematic model\\robot_trans.py

        Performs robot-specific transformation to constant curvature continuum robot arms with 3 ways of muscles.

        Adapted from FESTO paper [2015 T-RO], can be applied to our prototype ESNAKE.

        Input: q = [q1 q2 q3] is the length vector. representing the length of 3 muscles of one section. Unit: mm.
               r_b is the distance of the robot center and one bellow's center. Should be a fixed parameter. Unit: mm.
               
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

    """kinematic model\\robot_trans.py

        Transforms constant curvature configuration set to corresponding homogeneous transformation matrix. 

        This is the well-recognized way to perform transformation. (5 D-H transforms)

        Input: k = [phi kappa length] is the constant curvature configuration set. Unit: mm.

        Output: H is the homogeneous transformation matrix of the end tip from the bottom.
    """

    phi,kappa,leng = k[0],k[1],k[2]
    kl = kappa * leng

    H = np.matlib.zeros((4,4))

    if abs(kappa) < 1e-6: # kappa is sufficiently small that can be recognized as zero
        
        H[0,0] = (np.cos(phi) ** 2) * (np.cos(kl) - 1) + 1
        H[0,1] = np.sin(phi) * np.cos(phi) * (np.cos(kl) - 1)
        H[0,2] = np.cos(phi) * np.sin(kl)
        H[0,3] = np.cos(phi) * (1-np.cos(kl)) / kappa

        H[1,0] = np.sin(phi) * np.cos(phi) *(np.cos(kl) - 1)
        H[1,1] = (np.cos(phi) ** 2) * (1 - np.cos(kl)) + np.cos(kl)
        H[1,2] = np.sin(phi) * np.sin(kl)
        H[1,3] = np.sin(phi) * (1-np.cos(kl)) / kappa

        H[2,0] = -np.cos(phi) * np.sin(kl)
        H[2,1] = -np.sin(phi) * np.sin(kl)
        H[2,2] = np.cos(kl)
        H[2,3] = np.sin(kl) / kappa

        H[3,3] = 1
    
    else:
        H[0,0],H[1,1],H[2,2],H[3,3],H[2,3] = 1,1,1,1,leng

    return H


def config2rotMat(k):

    """kinematic model\\robot_trans.py

        Transforms constant curvature configuration set to corresponding homogeneous transformation matrix. And take
        the rotation matrix part of it. Representing the orientation of the tip. Can be used for calculating angular
        velocity. 

        This is the well-recognized way to perform transformation. (5 D-H transforms)

        Input: k = [phi kappa length] is the constant curvature configuration set. Unit: mm.

        Output: R is the rotational matrix of the end tip from the bottom.
    """
    
    H = config2homoTrans(k)

    R = np.split(H,3,axis=0)[0]
    R = np.split(R,3,axis=1)[0]

    return R


def config2tipPosition(k):

    """kinematic model\\robot_trans.py

        Transforms constant curvature configuration set to corresponding homogeneous transformation matrix. And take
        the translational part of it. Representing the position of the tip. Can be used for calculating translational
        velocity. 

        This is the well-recognized way to perform transformation. (5 D-H transforms)

        Input: k = [phi kappa length] is the constant curvature configuration set. Unit: mm.

        Output: r is the position vector of the end tip from the bottom.
    """

    H = config2homoTrans(k)
    r = H[0:2,3]

    return r


def actuator2homoTrans(q,r_b):
    """kinematic model\\robot_trans.py

        Combined robot-independent and robot-specific mapping from muscle length to H matrix.

        Input: q = [q1 q2 q3] is the length vector. representing the length of 3 muscles of one section. Unit: mm.
                r_b is the distance of the robot center and one bellow's center. Should be a fixed parameter. Unit: mm.
        
        Output: H is the homogeneous transformation matrix of the end tip from the bottom.
    """
    k = actuator2config(q, r_b)
    H = config2homoTrans(k)

    return H


def actuator2rotMat(q,r_b):
    """kinematic model\\robot_trans.py

        Combined robot-independent and robot-specific mapping from muscle length to R matrix.

        Input: q = [q1 q2 q3] is the length vector. representing the length of 3 muscles of one section. Unit: mm.
                r_b is the distance of the robot center and one bellow's center. Should be a fixed parameter. Unit: mm.
    
        Output: R is the rotational matrix of the end tip from the bottom.
    """
    k = actuator2config(q, r_b)
    R = config2rotMat(k)

    return R


def actuator2tipPosition(q,r_b):
    """kinematic model\\robot_trans.py
    
        Combined robot-independent and robot-specific mapping from muscle length to tip position.

        Input: q = [q1 q2 q3] is the length vector. representing the length of 3 muscles of one section. Unit: mm.
                r_b is the distance of the robot center and one bellow's center. Should be a fixed parameter. Unit: mm.

        Output: r is the position vector of the end tip from the bottom.
    """
    k = actuator2config(q, r_b)
    r = config2tipPosition(k)

    return r

if __name__ == '__main__':

    # TODO: Test the functions above.
    pass