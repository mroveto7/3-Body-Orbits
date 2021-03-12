import string ##needed to do string stuff
import sys    ##needed for input argument commands

fname = sys.argv[1] ##filename is the command line argument right after the program name

f = open(fname,'r')##opens the file with read permission only
w = f.readlines()##makes a list, w, that is just all the lines
f.close

numpts = len(w)##sets numpts equal to the first element of list w
z = []##need to create an array of length zero to store each element of list w

for i in range(0,len(w)):##loop over elements of list w, append the array of z
    z.append(w[i])##each element of z now consists of a single LINE of original data file
    
t=[]
x1=[]
y1=[]
x2=[]
y2=[]
x3=[]
y3=[]
    
    
for i in range(0,len(w)):##loop over our array of lines z
    tmp = str(z[i])##casts each element of z as a string so we can split it
    tmp2 = tmp.split()##creates a temporary array of length equal to number of elements in each line
    t.append(float(tmp2[0]))##appends x by the first element of tmp2, which is first column of data file
    x1.append(float(tmp2[1]))##appends y by the second element of tmp2, which is second column of data file
    y1.append(float(tmp2[2]))
    x2.append(float(tmp2[3]))
    y2.append(float(tmp2[4]))
    x3.append(float(tmp2[5]))
    y3.append(float(tmp2[6]))

fbody1 = open('body1.dat','w')
    
for i in range(0,len(w)):
    fbody1.write(str(x1[i]) + '\t' + str(y1[i]) + '\n')

fbody1.close()

fbody2 = open('body2.dat','w')
    
for i in range(0,len(w)):
    fbody2.write(str(x2[i]) + '\t' + str(y2[i]) + '\n')

fbody2.close()

fbody3 = open('body3.dat','w')
    
for i in range(0,len(w)):
    fbody3.write(str(x3[i]) + '\t' + str(y3[i]) + '\n')

fbody3.close()
