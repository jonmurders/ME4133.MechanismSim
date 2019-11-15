from math import pi
import math
import numpy as np
from numpy.linalg import inv
import pandas as pd
from matplotlib import pyplot as plt

sin = math.sin
cos = math.cos
tan = math.tan
#Bring In Input CSV
input = pd.read_csv('Input.csv')

#loading scalar knowns from input sheet
R1 = input.loc[0,'Value']  #inches - pg 96
R2 = input.loc[1,'Value'] #inches - pg 96
R6 = input.loc[5,'Value'] #inches - pg 96
theta_1 = input.loc[6,'Value']# pg 96
theta_5 = input.loc[10,'Value']# pg 96
theta_6 = input.loc[11,'Value']# pg 96
m2 = input.loc[12,'Value']
m3 = input.loc[13,'Value']
m4 = input.loc[14,'Value']
m5 = input.loc[15,'Value']
i2 = input.loc[16,'Value']
i3 = input.loc[17,'Value']
i4 = input.loc[18,'Value']
i5 = input.loc[19,'Value']
#Initial guess values from the input sheet
theta_3 = input.loc[8,'Value'] #radians
R3 = input.loc[2,'Value'] #meter
theta_4 = input.loc[9,'Value'] #radians
R4 = input.loc[3,'Value'] #meters
R5 =  input.loc[4,'Value'] #meters
x = np.array([theta_3,R3,theta_4,R4,R5], dtype=np.float)
#Input Setup from the input sheet
a2 = input.loc[20,'Value'] #radians per second per second
theta_2 = input.loc[7,'Value'] #radians
w2 = input.loc[22,'Value'] #radians per second
#Data collection table
positions = pd.DataFrame(columns=['theta_2','theta_3','R3','theta_4','R4','R5','h3','f3','h4','f4','f5','h3p','f3p','h4p','f4p','f5p','T','F12','F14','F23','F34','F45'])

r=0 #row Counter

while  theta_2 < input.loc[21,'Value']:
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
        f1 = R2*ct2-R3*ct3+R1 #R1x
        f2 = R2*st2-R3*st3 #R1y
        f3 = R5 - R4*ct4-R1 #R2x
        f4 = R6 -R4*st4 #R2y
        f5 = theta_4-theta_3 #G1
        f = [f1, f2, f3, f4, f5]
        fa = np.array(f,dtype=np.float)
        #finding the derivatives
        dfdh3 = np.array([[R3*st3], [-R3*ct3],[0], [0], [-1]], dtype=np.float)
        dfdr3 = np.array([[-ct3], [-st3], [0], [0], [0]], dtype=np.float)
        dfdt4 = np.array([[0], [0], [R4*st4], [-R4*ct4], [1]], dtype=np.float)
        dfdr4 = np.array([[0], [0], [-ct4], [-st4], [0]], dtype=np.float)
        dfdr5 = np.array([[0], [0], [0], [-1], [0]], dtype=np.float)
        #Making 5x5 array of derivatives
        A = np.hstack((dfdh3,dfdr3,dfdt4,dfdr4,dfdr5))
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
    #Logging Data into the table
    positions.loc[r,'theta_2'] = theta_2
    positions.loc[r,'theta_3'] = theta_3
    positions.loc[r,'R3'] = R3
    positions.loc[r,'theta_4'] = theta_4
    positions.loc[r,'R4'] = R4
    positions.loc[r,'R5'] = R5

    # velocity coefficents matricies
    B = np.array([[R3*st3, 0 , -ct3, 0, 0],
                 [-R3*ct3, 0, -st3, 0, 0],
                 [0, R4*ct4, 0, -ct4, ct5],
                 [0, -R4*st4, 0, -st4, st5],
                 [-1, 1, 0, 0, 0,]], dtype=np.float)

    C = np.array([[R2*st2],[-R2*ct2],[0],[0],[0]], dtype=np.float)

    binv = inv(B)
    y = binv*C
    #storing kinematic coefficent values
    positions.loc[r,'h3'] = y[0][0]
    positions.loc[r,'f3'] = y[1][0]
    positions.loc[r,'h4'] = y[2][0]
    positions.loc[r,'f4'] = y[3][0]
    positions.loc[r,'f5'] = y[4][0]
    h3 = y[0][0]
    f3 = y[1][0]
    h4 = y[2][0]
    f4 = y[3][0]
    f5 = y[4][0]
    #solving for second kinematic coefficents
    B = np.array([[R3*st3, 0 , -ct3, 0, 0],
                 [-R3*ct3, 0, -st3, 0, 0],
                 [0, R4*st4, 0, -ct4, ct5],
                 [0, -R4*st4, 0, -st4, st5],
                 [-1, 1, 0, 0, 0,]], dtype=np.float)

    C = np.array([[R2*ct2 - 2*f3*h3*ct3 - R3*st3*h3**2],[R2*st2 - 2*f3*h3*ct3 - R3*st3*h3**2],[-R4*ct4*h4**2],[2*f4*h4*ct4 - R4*st4*h4**2],[0]], dtype=np.float)

    binv = inv(B)
    y = binv*C
    #storing second kinematic coefficents
    positions.loc[r,'h3p'] = y[0][0]
    positions.loc[r,'f3p'] = y[1][0]
    positions.loc[r,'h4p'] = y[2][0]
    positions.loc[r,'f4p'] = y[3][0]
    positions.loc[r,'f5p'] = y[4][0]

    h3p = y[0][0]
    f3p = y[1][0]
    h4p = y[2][0]
    f4p = y[3][0]
    f5p = y[4][0]
    #IDP Analysis Start

    # kinematic coefficents of CGs
    fg4x = -(R4*.5)*st3*h3
    fg4y = (R4*.5)*ct3*h3
    fg3x = -R2*st2
    fg3y = R2*ct2
    fg2x = -(R2*.5)*st2
    fg2y = (R2*.5)*ct2
    fg5x = f5
    fg5y = 0

    #Second Order Kinematic Coefficients
    fg4xp = (-(R4*.5)*ct3*h3**2)+(-(R4*.5)*st3*h3p)
    fg4yp = (-(R4*.5)*st3*h3**2)+((R4*.5)*ct3*h3p)
    fg3xp = -R2*ct2
    fg3yp = -R2*st2
    fg2xp = -(R2*.5)*ct2
    fg2yp = -(R2*.5)*st2
    fg5xp = f5p
    fg5yp = 0

    #Accelerations
    ag4x = fg4x*a2+fg4xp*w2**2
    ag4y = fg4y*a2+fg4yp*w2**2
    ag3x = fg3x*a2+fg3xp*w2**2
    ag3y = fg3y*a2+fg3yp*w2**2
    ag2x = fg2x*a2+fg2xp*w2**2
    ag2y = fg2y*a2+fg2yp*w2**2
    ag5x = fg5x*a2+fg5xp*w2**2
    ag5y = fg5y*a2+fg5yp*w2**2
    a3 = h3*w2+h3p*w2**2
    ag2 = ((ag2x**2)+ag2y**2)**.5
    ag4 = ((ag4x**2)+ag4y**2)**.5
    #Solving forces
    F45x = m5*ag5x
    F45y = m5*ag5y - m5*9.81
    F45 = ((F45x**2)+F45y**2)**.5
    F34 = ((i4*a3)+(.5*R4*m4*ag4)-(.5*R4*m4*9.81*st3)+(.5*R4*m4*9.81*ct3)-(F45*R4))/R3
    F34x = F34*cos(180-theta_3)
    F34y = F34*sin(180-theta_3)
    F23x = m3*ag3x - F34x
    F23y = m3*ag3x - F34y - m3*9.81
    F23 = ((F23x**2)+F23y**2)**.5
    F14x = m4*ag4x - F45x - F34x
    F14y = m4*ag4y - F45y - F34y - m4*9.81
    F14 = ((F14x**2)+F14y**2)**.5
    F12x = m2*ag2x - F23x
    F12y = m2*ag2x - F23y - m3*9.81
    F12 = ((F12x**2)+F12y**2)**.5
    T = i2*a2+(.5*R2*m2*ag2)-(.5*R2*m2*9.81)+R2*F23

    #Storing Values at theta 2
    positions.loc[r,'theta_2'] = theta_2
    positions.loc[r,'T'] = T
    positions.loc[r,'F12'] = F12
    positions.loc[r,'F14'] = F14
    positions.loc[r,'F23'] = F23
    positions.loc[r,'F34'] = F34
    positions.loc[r,'F45'] = F45

    theta_2 +=.01
    r += 1

#Generating t2 vs h3 plot
plt.figure(1)
plt.plot(positions.theta_2,positions.theta_3)
titlet3 = 'Theta 2 vs Theta3'
plt.title(titlet3)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Theta 3 (radians)')
plt.savefig(titlet3)

#Generating t2 vs R3 plot
plt.figure(2)
plt.plot(positions.theta_2,positions.R3)
titler3 = 'Theta 2 vs Vector R3'
plt.title('Theta 2 vs Vector R3')
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Vector R3 (m)')
plt.savefig(titler3)

#Generating t2 vs t4 plot
plt.figure(3)
plt.plot(positions.theta_2,positions.theta_4)
titlet4 = 'Theta 2 vs Theta 4'
plt.title(titlet4)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Theta 4 (radians)')
plt.savefig(titlet4)

#Generating t2 vs R4 plot
plt.figure(4)
plt.plot(positions.theta_2,positions.R4)
titler4 = 'Theta 2 vs Vector R4'
plt.title(titler4)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Vector R4 (m)')
plt.savefig(titler4)

#Generating t2 vs R5 plot
plt.figure(5)
plt.plot(positions.theta_2,positions.R5)
titler5 = 'Theta 2 vs Vector R5'
plt.title(titler5)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Vector R5 (m)')
plt.savefig(titler5)

#Generating t2 vs h3 plot
plt.figure(6)
plt.plot(positions.theta_2,positions.h3)
titleh3 = 'Theta 2 vs h3'
plt.title(titleh3)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('h3 (radians/s)')
plt.savefig(titleh3)

#Generating t2 vs R3 plot
plt.figure(7)
plt.plot(positions.theta_2,positions.f3)
titler7 = 'Theta 2 vs f3'
plt.title('Theta 2 vs f3')
plt.xlabel('Theta 2 (radians)')
plt.ylabel('r3 (m/s)')
plt.savefig(titler7)

#Generating t2 vs t4 plot
plt.figure(8)
plt.plot(positions.theta_2,positions.h4)
titlet8 = 'Theta 2 vs h4'
plt.title(titlet4)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('h4 (radians/s)')
plt.savefig(titlet8)

#Generating t2 vs R4 plot
plt.figure(9)
plt.plot(positions.theta_2,positions.f4)
titler9 = 'Theta 2 vs f4'
plt.title(titler4)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('f4 (m/s)')
plt.savefig(titler9)

#Generating t2 vs R5 plot
plt.figure(10)
plt.plot(positions.theta_2,positions.f5)
titler10 = 'Theta 2 vs f5'
plt.title(titler10)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('f5 (m/s)')
plt.savefig(titler10)

#Generating t2 vs h3' plot
plt.figure(11)
plt.plot(positions.theta_2,positions.h3p)
titleh3p = 'Theta 2 vs h3\''
plt.title(titleh3p)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('h3\' (radians/s)')
plt.savefig(titleh3p)

#Generating t2 vs R3' plot
plt.figure(12)
plt.plot(positions.theta_2,positions.f3p)
titler3p = 'Theta 2 vs f3\''
plt.title(titler3p)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('r3\' (m/s)')
plt.savefig(titler3p)

#Generating t2 vs t4 plot
plt.figure(13)
plt.plot(positions.theta_2,positions.h4p)
titlet4p = 'Theta 2 vs h4\''
plt.title(titlet4p)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('h4\' (radians/s)')
plt.savefig(titlet4p)

#Generating t2 vs R4 plot
plt.figure(14)
plt.plot(positions.theta_2,positions.f4p)
titler4p = 'Theta 2 vs f4\''
plt.title(titler4p)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('f4\' (m/s)')
plt.savefig(titler4p)

#Generating t2 vs R5 plot
plt.figure(15)
plt.plot(positions.theta_2,positions.f5p)
titler5p = 'Theta 2 vs f5\''
plt.title(titler5p)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('f5\' (m/s)')
plt.savefig(titler5p)

#Generating t2 vs T plot
plt.figure(16)
plt.plot(positions.theta_2,positions['T'])
title16 = 'Theta 2 vs T'
plt.title(title16)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Torque (kNm)')
plt.savefig(title16)

#Generating t2 vs F12 plot
plt.figure(17)
plt.plot(positions.theta_2,positions.F12)
title17 = 'Theta 2 vs F12'
plt.title(title17)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Force (kN)')
plt.savefig(title17)

#Generating t2 vs F14 plot
plt.figure(18)
plt.plot(positions.theta_2,positions.F14)
title18 = 'Theta 2 vs F14'
plt.title(title18)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Force (kN)')
plt.savefig(title18)

#Generating t2 vs F23 plot
plt.figure(19)
plt.plot(positions.theta_2,positions.F23)
title19 = 'Theta 2 vs F23'
plt.title(title19)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Force (kN)')
plt.savefig(title19)

#Generating t2 vs F34 plot
plt.figure(20)
plt.plot(positions.theta_2,positions.F34)
title20 = 'Theta 2 vs F34'
plt.title(title20)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Force (kN)')
plt.savefig(title20)

#Generating t2 vs F45 plot
plt.figure(21)
plt.plot(positions.theta_2,positions.F45)
title21 = 'Theta 2 vs F45'
plt.title(title21)
plt.xlabel('Theta 2 (radians)')
plt.ylabel('Force (kN)')
plt.savefig(title21)
