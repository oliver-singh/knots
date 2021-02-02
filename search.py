import sys
import os
from sympy.abc import t

sys.path.append(os.path.dirname(__file__))
from braid import *

def DV(n):
    knot = ( full_twist(n) *
            (empty(2) + full_twist(-n+2)) *
            full_twist(-3) * Braid([-1,-1]) * (Braid([3,4])**-3) * Braid([-5,-5]) *
             -braid_range(n))
    knot = knot.simplify()
    return knot
def figeight(n):
    knot = Braid([3,2,1,4,3,2,2,3,1,2,-1,-1,-4,-3,-5,-5]) * (-braid_range(n))
    return knot

def tref(n):
    knot = Braid([4,3,2,5,4,3,3,4,2,3,1,-2,1,-5,-4,-6,-6]) * (-braid_range(2,n))
    return knot

def twisted_braid(twists, displacements, knot_n_function, n=11):
    knot = Braid()
    for twist, displacement in zip(twists,displacements):
        assert abs(twist) + displacement <= n, "too many strands!"
        knot = knot * (empty(displacement) + full_twist(twist))
    knot = knot * knot_n_function(n)
    knot = knot.super_simplify()
    return knot






if __name__ == "__main__":
    for i in range(5,11,2):
        for i_dis in range(11-i+1):
            for j in range(5, 11,2):
                for j_dis in range(11 - j+1):
                    fr = 11**2 - 9 **2 - i**2 + j**2
                    if fr > 0:
                        print("********")
                        knot = twisted_braid([11,-9,-i,j], [0,2,i_dis,j_dis], figeight)
                        
                        
                        print("twists: ",[11,-9, -i, j], "displaced by ",[0,2,i_dis, j_dis])
                        obs = 4+(1/4)*fr -4
                        print("obstruction:", 4+(1/4)*fr -4)
                        pol = knot.alexander_poly()
                        degree = pol.degree(t)
                        print("degree", degree)
                        if degree < obs:
                            print("*****************************************success!!!!!*********************************************************")