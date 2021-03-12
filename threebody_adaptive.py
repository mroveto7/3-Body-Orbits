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


def rk4(x, dt, gravder, m, par, Ravg, delta, test1):

    k1 = np.zeros((3,4))
    k2 = np.zeros((3,4))
    k3 = np.zeros((3,4))
    k4 = np.zeros((3,4))
    k5 = np.zeros((3,4))
    k6 = np.zeros((3,4))
    qtemp = np.zeros((3,4))
    test2 = np.zeros((3,4))
    R = np.zeros((3,2))
    tol = 0.0001

    for i in range(0,3):
        for j in range(0,4):
            qtemp[i,j] = x[i,j]

    c1 = gravder(qtemp, m, par, k1)
    for i in range(0,3):
        for j in range(0,4):
            k1[i,j] *= dt
            qtemp[i,j] = x[i,j] + k1[i,j]/4

    c2 = gravder(qtemp, m, par, k2)
    for i in range(0,3):
        for j in range(0,4):
            k2[i,j] *= dt
            qtemp[i,j] = x[i,j] + (3/32)*k1[i,j] + (9/32)*k2[i,j]

    c3 = gravder(qtemp, m, par, k3)
    for i in range(0,3):
        for j in range(0,4):
            k3[i,j] *= dt
            qtemp[i,j] = x[i,j] + (1932/2197)*k1[i,j] + (7200/2197)*k2[i,j] + (7296/2197)*k3[i,j]

    c4 = gravder(qtemp, m, par, k4)
    for i in range(0,3):
        for j in range(0,4):
            k4[i,j] *= dt
            qtemp[i,j] = x[i,j] + (439/216)*k1[i,j] - 8*k2[i,j] + (3680/513)*k3[i,j] - (845/4104)*k4[i,j]

    c5 = gravder(qtemp, m, par, k5)
    for i in range(0,3):
        for j in range(0,4):
            k5[i,j] *= dt
            qtemp[i,j] = x[i,j] - (8/27)*k1[i,j] + 2*k2[i,j] - (3544/2565)*k3[i,j] + (1859/4104)*k4[i,j] - (11/40)*k5[i,j]

    c6 = gravder(qtemp, m, par, k6)
    for i in range(0,3):
        for j in range(0,4):
            k6[i,j] *= dt

    for i in range(0,3):
        for j in range(0,4):
            test1[i,j] = x[i,j] + (25/216)*k1[i,j] + (1408/2565)*k3[i,j] + (2197/4104)*k4[i,j] - k5[i,j]/5
            test2[i,j] = x[i,j] + (16/135)*k1[i,j] + (6656/12825)*k3[i,j] + (28561/56430)*k4[i,j] - (9/50)*k5[i,j] + (2/55)*k6[i,j]

    for i in range(0,3):
        R[i,0] = (1/dt)*np.absolute(test2[i,0] - test1[i,0])
        R[i,1] = (1/dt)*np.absolute(test2[i,1] - test1[i,1])
    
    Ravg = (R[0,0] + R[0,1] + R[1,0] + R[1,1] + R[2,0] + R[2,1])/2
    delta = 0.84*np.power((tol/Ravg),(1/4))

    return Ravg, delta

####functions defined, new stuff follows####


q = np.zeros((3,4))
q0 = np.zeros((3,4))
m = np.zeros(3)
t = 0.0
dt = 1./10000.
par = np.zeros(1)
ns = 100000000
tol = 0.000055

par[0] = 1

test1 = np.zeros((3,4))
Ravg = 0.0
delta = 0.0

file = open('three_body.dat','w')

#####################################
#Special Solutions
#####################################
#Butterfly II
#p1 = 0.392955
#p2 = 0.097579
#####################################
#Yin-Yang 2b
#p1 = 0.417343
#p2 = 0.313100
#####################################
#Dragonfly
#p1 = 0.080584
#p2 = 0.588836
#####################################
#Butterfly IV
p1 = 0.350112
p2 = 0.079339

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
        
n = 0
while n < ns:
   for i in range(0,3):
       file.write(str(q[i,0]) + '\t' + str(q[i,1]) + '\t')
   file.write('\n')

   Ravg, delta = rk4(q, dt, gravder, m, par, Ravg, delta, test1)
   #print(Ravg)
   #print(delta)
   if Ravg <= tol:
       for i in range(0,3):
           for j in range(0,4):
               q[i,j] = test1[i,j]
           dt *= delta
           n += 1
   else:
       dt *= delta
       n = n
       #print('it happened')
           

file.close()
