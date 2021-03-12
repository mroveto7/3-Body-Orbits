import string
import sys
import numpy

np = numpy

def gravder(x, m, par, deriv):
    G = par[0]
    r = numpy.zeros(3)

    r[0] = np.sqrt((x[2,0] - x[1,0])*(x[2,0] - x[1,0]) + (x[2,1] - x[1,1])*(x[2,1] - x[1,1]))
    r[1] = np.sqrt((x[0,0] - x[2,0])*(x[0,0] - x[2,0]) + (x[0,1] - x[2,1])*(x[0,1] - x[2,1]))
    r[2] = np.sqrt((x[1,0] - x[0,0])*(x[1,0] - x[0,0]) + (x[1,1] - x[0,1])*(x[1,1] - x[0,1]))

    deriv[0,0] = x[0,2]
    deriv[0,1] = x[0,3]
    deriv[0,2] = G*m[1]*((x[1,0] - x[0,0])/np.power(r[2],3)) + G*m[2]*((x[2,0] - x[0,0])/np.power(r[1],3))
    deriv[0,3] = G*m[1]*((x[1,1] - x[0,1])/np.power(r[2],3)) + G*m[2]*((x[2,1] - x[0,1])/np.power(r[1],3))
    deriv[1,0] = x[1,2]
    deriv[1,1] = x[1,3]
    deriv[1,2] = G*m[2]*((x[2,0] - x[1,0])/np.power(r[0],3)) + G*m[0]*((x[0,0] - x[1,0])/np.power(r[2],3))
    deriv[1,3] = G*m[2]*((x[2,1] - x[1,1])/np.power(r[0],3)) + G*m[0]*((x[0,1] - x[1,1])/np.power(r[2],3))
    deriv[2,0] = x[2,2]
    deriv[2,1] = x[2,3]
    deriv[2,2] = G*m[0]*((x[0,0] - x[2,0])/np.power(r[1],3)) + G*m[1]*((x[1,0] - x[2,0])/np.power(r[0],3))
    deriv[2,3] = G*m[0]*((x[0,1] - x[2,1])/np.power(r[1],3)) + G*m[1]*((x[1,1] - x[2,1])/np.power(r[0],3))

    return 0


def rk4(x, dt, gravder, m, par):

    k1 = np.zeros((3,4))
    k2 = np.zeros((3,4))
    k3 = np.zeros((3,4))
    k4 = np.zeros((3,4))
    qtemp = np.zeros((3,4))

    for i in range(0,3):
        for j in range(0,4):
            qtemp[i,j] = x[i,j]

    c1 = gravder(qtemp, m, par, k1)
    for i in range(0,3):
        for j in range(0,4):
            k1[i,j] *= dt
            qtemp[i,j] = x[i,j] + k1[i,j]/2

    c2 = gravder(qtemp, m, par, k2)
    for i in range(0,3):
        for j in range(0,4):
            k2[i,j] *= dt
            qtemp[i,j] = x[i,j] + k2[i,j]/2

    c3 = gravder(qtemp, m, par, k3)
    for i in range(0,3):
        for j in range(0,4):
            k3[i,j] *= dt
            qtemp[i,j] = x[i,j] + k3[i,j]

    c4 = gravder(qtemp, m, par, k4)
    for i in range(0,3):
        for j in range(0,4):
            k4[i,j] *= dt

    for i in range(0,3):
        for j in range(0,4):
            x[i,j] += (k1[i,j] + 2*k2[i,j] + 2*k3[i,j] + k4[i,j])/6

    return c1 + c2 + c3 + c4

####functions defined, new stuff follows####


q = np.zeros((3,4))
q0 = np.zeros((3,4))
m = np.zeros(3)
t = 0.0
dt = 0.0001
par = np.zeros(1)
ns = 80000

par[0] = 1

file = open('three_body.dat','w')
####################################
#Sun-Earth-Moon
#par[0] = 4*np.pi*np.pi
#  
#m[0] = 1
#m[1] = 3.003e-6
#m[2] = 3.694e-8
#  
#  
#q0[0,0] = 0
#q0[0,1] = 0
#q0[0,2] = 0
#q0[0,3] = 0
#q0[1,0] = 1
#q0[1,1] = 0
#q0[1,2] = 0
#q0[1,3] = 6.28
#q0[2,0] = 1
#q0[2,1] = 0.00257
#q0[2,2] = -0.1937
#q0[2,3] = 6.28
#####################################

#####################################
#Equal Mass - Lagrange Solution
#m[0] = 1
#m[1] = 1
#m[2] = 1
#
#q0[0,0] = -0.5
#q0[0,1] = np.sqrt(3)/6
#q0[0,2] = -np.sqrt(3.)/6
#q0[0,3] = -0.5
#q0[1,0] = 0.5
#q0[1,1] = np.sqrt(3)/6
#q0[1,2] = -np.sqrt(3)/6
#q0[1,3] = 0.5
#q0[2,0] = 0
#q0[2,1] = -np.sqrt(3)/3
#q0[2,2] = np.sqrt(3)/3
#q0[2,3] = 0

#####################################

#####################################
#Unequal Mass - Lagrange Solution
#m[0] = 1
#m[1] = 2
#m[2] = 3
#
#q0[0,0] = -7/12
#q0[0,1] = 0.25*np.sqrt(3)
#q0[0,2] = -0.25*np.sqrt(3)
#q0[0,3] = -7/12
#q0[1,0] = 5/12
#q0[1,1] = 0.25*np.sqrt(3)
#q0[1,2] = -0.25*np.sqrt(3)
#q0[1,3] = 5/12
#q0[2,0] = -1/12
#q0[2,1] = -0.25*np.sqrt(3)
#q0[2,2] = 0.25*np.sqrt(3)
#q0[2,3] = -1/12

#####################################

#####################################
#Equal Mass - Euler Solution
#m[0] = 1
#m[1] = 1
#m[2] = 1
#
#q0[0,0] = -1      #x1
#q0[0,1] = 0       #y1
#q0[0,2] = 0       #vx1
#q0[0,3] = -6.5    #vy1
#q0[1,0] = 0       #x2
#q0[1,1] = 0       #y2
#q0[1,2] = 0       #vx2
#q0[1,3] = 0       #vy2
#q0[2,0] = 1       #x3
#q0[2,1] = 0       #y3
#q0[2,2] = 0       #vx3
#q0[2,3] = 6.5     #vy3

#####################################

#####################################
#Figure Eight
#m[0] = 1.;
#m[1] = 1.;
#m[2] = 1.;
#
#q0[0,0] = 0.97000436
#q0[0,1] = -0.24308753
#q0[0,2] = 0.93240737/2
#q0[0,3] = 0.86473146/2
#q0[1,0] = -0.97000436
#q0[1,1] = 0.24308753
#q0[1,2] = 0.93240737/2
#q0[1,3] = 0.86473146/2
#q0[2,0] = 0
#q0[2,1] = 0
#q0[2,2] = -0.93240737
#q0[2,3] = -0.86473146

#####################################

#####################################
#Special Solutions
#####################################

#####################################
#MOTH I - Equal Mass
#p1 = 0.46445
#p2 = 0.396060

#####################################

#####################################
#Yin-Yang 2b - Equal Mass
#p1 = 0.417343
#p2 = 0.313100

#####################################

#####################################
#Butterfly II
p1 = 0.392955
p2 = 0.097579

#####################################

m[0] = 1
m[1] = 1
m[2] = 1

q0[0,0] = -1
q0[0,1] = 0
q0[0,2] = p1
q0[0,3] = p2
q0[1,0] = 1
q0[1,1] = 0
q0[1,2] = p1
q0[1,3] = p2
q0[2,0] = 0
q0[2,1] = 0
q0[2,2] = -2*p1
q0[2,3] = -2*p2

#####################################

 
for i in range(0,3):
    for j in range(0,4):
        q[i,j] = q0[i,j]

for n in range(0,ns):
   file.write(str(t) + '\t')
   for i in range(0,3):
       file.write(str(q[i,0]) + '\t' + str(q[i,1]) + '\t')
   file.write('\n')

   rk4(q, dt, gravder, m, par)
   t += dt

file.close()
