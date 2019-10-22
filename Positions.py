from math import math.pi as pi
import math

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
R3 = 3 #inches
R4 = 11 #inches
R5 =  4 #inches
theta_3 = 5.236 #inches
theta_4 = 5.236 #inches

#Input Angle
theta_2 = 4.94 #inches

#Loop Counter
i = 0

#Main Loop
while i<6:
    #finding the sines and cosines of all the angles
    ct2 = cos(theta_2)
    st2 = sint(theta_2)
    ct3 = cos(theta_3)
    st3 = sint(theta_3)
    ct4 = cos(theta_4)
    st4 = sint(theta_4)
    ct5 = cos(theta_5)
    st5 = sint(theta_5)
    ct6 = cos(theta_6)
    st6 = sint(theta_6)

    #find the values of the VLEs provided on page 96
    f1 = R2*ct2-R3*ct3+R1
    f2 = R2*st2-R3*st3
    f3 = R6-R4*ct4+R1
    f4 = -R5-R4*st4
    f5 = theta_4-theta_3
    f = [f1, f2, f3, f4, f5]
