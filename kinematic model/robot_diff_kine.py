# First written for Graduation Project 'Research on dynamic modelling of a soft robot arm with peumatic bellows actuators and a
# damping support structure'. Functions in this module deals with differential kinematics of the robot.
#
# Author: Yongxi Cao  (Affiliation when doing the project: Dept. Mechincal & Energy Engineering, SUSTech)
# Email: yongxi.cao.2000@gmail.com
# v1.0: 2022-Mar-12
#
# function list:
#   1. Gamma(aster)
#   2. dGammad(aster,dWhat)
#   3. dHdq(H,q,i)


def Gamma(aster):

    """kinematic_model\\robot_diff_kine.py
    A linear operator used for convinently calculating ratation rates.
    
    Input: Some 3*3 matrix.
    
    Output: Corresponding rasult after gamma operator.
    """

    preMut = np.matlib.zeros((4,4))

    return diag