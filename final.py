"""
Assignment 2 EE2703: Applied Programming Lab
Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN2_EE19B077.py ,filename..netlist
Output: print the node voltages as well as the current through voltaage sources
"""
import numpy as np
import math
import sys


#Class that contains the information about all components present in the circuit
class Components:
    
    name=""
    Node_1=0
    Node_2=0
    Value=0
    
    def __init__(self,i_name,i_Node1,i_Node2,i_value):
        self.name = i_name
        self.Node_1 = i_Node1
        self.Node_2 = i_Node2
        self.Value = i_value

#Function to solve the nodal equations
def Solve():
    
    KCL = []    #Array of all nodal equations
    volt = []   #Voltage value of the voltage source
    current = []#Current Value of the Current source
    i=0         #i everywhere used as Control variable
   
    while(i < (max_node+2)):
        KCL.append({'0':0}) #Creates an array of dictionaries for each node referenced as 0 to last node
        k=0                 #k everywhere used as Control variable
        while(k <= max_node):
            KCL[i][str(k)] = 0
            k=k+1
        i=i+1
    
    i=0
    while(i <= max_node):   
        current.append(0)
        i=i+1

    count=1     #Count of no. of voltage sources
    i=0
    while(i < len(Component)):  #Nodal equations are formed in this
        
        if("R" in Component[i].name or "C" in Component[i].name or "L" in Component[i].name):
            KCL[int(Component[i].Node_1)][Component[i].Node_1] = (1/(Component[i].Value))  + KCL[int(Component[i].Node_1)][Component[i].Node_1]
            KCL[int(Component[i].Node_1)][Component[i].Node_2] = KCL[int(Component[i].Node_1)][Component[i].Node_2] - (1/(Component[i].Value))
            KCL[int(Component[i].Node_2)][Component[i].Node_2] = (1/(Component[i].Value))  + KCL[int(Component[i].Node_2)][Component[i].Node_2]
            KCL[int(Component[i].Node_2)][Component[i].Node_1] = KCL[int(Component[i].Node_2)][Component[i].Node_1] - (1/(Component[i].Value))
        
        if("I" in Component[i].name):
            current[int(Component[i].Node_1)] = float(Component[i].Value)
            current[int(Component[i].Node_2)] = -float(Component[i].Value)
        
        if("V" in Component[i].name):
            KCL[int(Component[i].Node_1)][str(max_node+count)] = 1
            KCL[int(Component[i].Node_2)][str(max_node+count)] = -1
            KCL[(max_node + count)][Component[i].Node_2] = 1
            KCL[(max_node + count)][Component[i].Node_1] = -1
            volt.append(Component[i].Value)
            count=count+1
        
        i=i+1

    i=0
    
    #Set ground as zero
    while(i < (max_node+2)):
        KCL[i]['0'] = 0
        i=i+1
    KCL[0] = 0

    #   [M_ac]*[X] = b OR [M_dc]*[X] = b
    b = np.array([(max_node+count)])
    M_ac = np.empty([(max_node+count-1),(max_node+count-1)],dtype=complex)
    M_dc = np.empty([(max_node+count-1),(max_node+count-1)])
    temp = []
    i=0

    while(i < max_node):    #To get b
        temp.append(current[i])
        i=i+1
    while(i < (max_node+count-1)):
        temp.append(volt[(i - max_node)])
        i = i+1
    b=temp

    i=1

    while(i < (max_node+count)):    #To get M_ac and M_dc
        k=0
        while(k < (max_node+count-1)):
            if(ac == 1):
                try:
                    M_ac[i-1][k] = KCL[i][str(k+1)]
                except Exception:
                    M_ac[i-1][k] = 0
            elif(ac == 0):
                try:
                    M_dc[i-1][k] = round(KCL[i][str(k+1)],3)
                except Exception:
                    M_dc[i-1][k] = 0
            k=k+1
        i=i+1

    print("Voltage at GND: 0")

    if(ac == 1):    #Output M_ac
        try:
            out_ac = np.linalg.solve(M_ac,b)    #Function that solves the matrix
            i=0
            while(i < max_node):
                print("Voltage at node"+ str(i+1)+ ":"+ str(np.round(out_ac[i],10)))
                i=i+1
            while(i < max_node+count-1):
                print("Current through Voltage Source"+ str(i-(max_node+count-1))+ " :"+ str(np.round(out_ac[i],10)))
                i=i+1
        except Exception:
            print("Redundant netlist (check if integers are used to represent nodes in netlist)")

    if(ac == 0):    #Output M_dc
        try:
            out_dc = np.linalg.solve(M_dc,b)
            i=0
            while(i < max_node):
                print("Voltage at node"+ str(i+1)+ ":"+ str(np.round(out_dc[i],10)))
                i=i+1
            while(i < max_node+count-1):
                print("Current through Voltage Source" + str(i-(max_node+count-1))+ " :"+ str(np.round(out_dc[i],10)))
                i=i+1
        except Exception:
            print("Redundant netlist (check if integers are used to represent nodes in netlist)")

#Get input of file
if len(sys.argv) != 2:
    print('Wrong file name')
    sys.exit(0)
file_name = sys.argv[1]
try:
    file = open(file_name)
except Exception:
    print('Wrong file name')
    sys.exit(0)

#Reading file
data = file.read()
file.close()
line = data.split('\n')

i=0
ac=0            #Check whether circuit is ac or dc
freq = 0        #Frequence of circuit
Component = []  #Instance of class components
max_node=0      #Max No. of nodes in the circuit

while(i < len(line) and not('.ac' in line[i])):     #Check of ac or dc
    i=i+1
if(i != len(line)):
    ac=1
    word = line[i].split()
    freq = float((word[2]))

file = open(file_name)
data = file.read()
file.close()
line = data.split('\n')     #Each line of a file
i=0

start = -1        
finish = -1
len_cir = len(".circuit")
len_end = len(".end")

for k in line:
    if '.circuit' == k[:len_cir]:
        start = line.index(k)          #start index
    elif '.end' == k[:len_end]:
        finish = line.index(k)         #end index
    k = k.split('#')[0]
if start >= finish:                     #validating the code is write or not
    print('Invalid Spice Code')
    exit()

while(i < len(line) and line[i] != '.circuit'):     #Check for .circuit
    i = i+1
i = i+1

#Get the inputs from file and store in class after aniputlation of values as required
while(i < len(line) and line[i] != '.end'):
    
    word = line[i].split()      #Each word in the file
    
    if(word[1] == 'GND'):
        word[1] = '0'
    if(word[2] == 'GND'):
        word[2] = '0'

    if(int(word[1][-1]) > max_node):    #Find max no. of nodes
        max_node = int(word[1][-1])
    if(int(word[2][-1]) > max_node):
        max_node = int(word[2][-1])

    if('R' in word[0]):
        word[3] = float(word[3])

    if(ac == 0):
        if('V' in word[0]):
            word[3] = float(word[3])

    if(ac == 1):
        
        if( 'V' in word[0]):
            if( 'dc' in word[3]):
                print("ERROR!!! Cant do with a dc source on an AC circuit...Will be take as AC source with phase 0")
            elif( 'ac' in word[3]):
                try:
                    phase = (float(word[5])*math.pi)/180    #The phase of the voltage source
                except Exception:
                    phase = 0
                word[3] = complex((float(word[4])/2)*math.cos(phase),(float(word[4])/2)*math.sin(phase))

        elif( 'L' in word[0]):      #Finding impedamce
                word[3] = complex(0,(2*math.pi*freq*float(word[3])))
        elif( 'C' in word[0]):
            word[3] = complex(0,(1/(2*math.pi*freq*float(word[3]))))

    Component.append(Components(word[0],word[1][-1],word[2][-1],word[3]))   #Storing the object

    i=i+1

i=0
print()
Solve()     #Calling the function
print()


