from math import pi
import math
import numpy as np
from numpy.linalg import inv

sin = math.sin
cos = math.cos
tan = math.tan

#scalar knowns
R1 = 4.8 #inches - pg 96
R2 = 2 #inches - pg 96
R6 = 3.65 #inches - pg 96
theta_1 = -(pi)# pg 96
theta_5 = -(pi)/2# pg 96
theta_6 = 0# pg 96

#Initial guess values
theta_3 = 5.236 #inches
R3 = 3 #inches
theta_4 = 5.236 #inches
R4 = 11 #inches
R5 =  4 #inches
x = np.array([theta_3,R3,theta_4,R4,R5], dtype=np.float)
#Input Angle
theta_2 = 0 #radians



while theta_2 < (2*pi):
    #finding the sines and cosines of all the angles
    ct2 = cos(theta_2)
    st2 = sin(theta_2)
    ct3 = cos(theta_3)
    st3 = sin(theta_3)
    ct4 = cos(theta_4)
    st4 = sin(theta_4)
    ct5 = cos(theta_5)
    st5 = sin(theta_5)
    ct6 = cos(theta_6)
    st6 = sin(theta_6)
    #Loop Counter
    i = 0

    #Newton's Method Loop
    while i<100:


        #find the values of the VLEs provided on page 96
        f1 = R2*ct2-R3*ct3+R1
        f2 = R2*st2-R3*st3
        f3 = R6-R4*ct4+R1
        f4 = -R5-R4*st4
        f5 = theta_4-theta_3
        f = [f1, f2, f3, f4, f5]
        fa = np.array(f,dtype=np.float)
        #finding the derivatives
        dfdt3 = np.array([[R3*st3], [-R3*ct3],[0], [0], [-1]], dtype=np.float)
        dfdr3 = np.array([[-ct3], [-st3], [0], [0], [0]], dtype=np.float)
        dfdt4 = np.array([[0], [0], [R4*st4], [-R4*ct4], [1]], dtype=np.float)
        dfdr4 = np.array([[0], [0], [-ct4], [-st4], [0]], dtype=np.float)
        dfdr5 = np.array([[0], [0], [0], [-1], [0]], dtype=np.float)
        #Making 5x5 array of derivatives
        A = np.hstack((dfdt3,dfdr3,dfdt4,dfdr4,dfdr5))
        #Takes the inverse of Matrix A
        ainv = inv(A)
        #Newton's Method Applied
        x = x-(ainv*fa)
        #Extracts values
        theta_3 = x[0][0]
        R3 = x[1][0]
        theta_4 = x[2][0]
        R4 = x[3][0]
        R5 = x[4][0]
        i+=1
    theta_2 +=.1
    print(theta_3, R3, theta_4, R4, R5)
