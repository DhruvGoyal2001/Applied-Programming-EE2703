"""
Assignment 2 EE2703: Applied Programming Lab
Dhruv Goyal, EE19B077

Usage example: python3 EE2703_ASSIGN2_EE19B077.py ,filename..netlist
Output: print the node voltages as well as the current through voltaage sources
"""
from sys import argv,exit

from numpy import *

class element_info():
    ''' General class for all the one port elements '''
    def __init__(self,line):
        self.line = line
        self.words = self.line.split()
        self.name = elem_type(self.words[0])
        self.from_node = self.words[1]
        self.to_node = self.words[2]
        if len(self.words) == 5:
            self.type = 'dc'
            self.value = float(self.words[4])
        elif len(self.words) == 6:
            self.type = 'ac'
            Vm = float(self.words[4])/2
            phase = float(self.words[5])
            real = Vm*cos(phase)
            imag = Vm*sin(phase)
            self.value = complex(real,imag)
        else:
            self.type = 'RLC'
            self.value = float(self.words[3])


def elem_type(name):
    ''' Gets the element name '''
    if name[0] == 'R':
        return 'resistor'
    elif name[0] == 'L':
        return 'inductor'
    elif name[0] == 'C':
        return 'capacitor'
    elif name[0] == 'V':
        return 'voltage source'
    elif name[0] == 'I':
        return 'current source'
    elif name[0] == 'G':
        return 'voltage controlled voltage source'
    elif name[0] == 'E':
        return 'voltage controlled current source'
    elif name[0] == 'F':
        return 'current controlled current source'
    elif name[0] == 'H':
        return 'current controlled voltage source'
    else:
        print("Unknown element defined")
        exit()


def dictionary_for_node(cir_def):
    ''' Returns a dictionary of nodes from the circuit definition. By default, the GND node is assigned value 0. '''
    node_rep = {}
    nodes = [element_info(line).from_node for line in cir_def]
    nodes.extend([element_info(line).to_node for line in cir_def])
    ind = 1
    nodes = list(set(nodes))           #removes duplicate values
    for node in nodes:
        if node == 'GND' :
            node_rep[node] = 0
        else :
            node_rep[node] = ind
            ind += 1
    return node_rep

AC = '.ac'

def freq(lines):
    ''' Returns the frequency of the source '''
    f = 0
    for line in lines:
        if line[:len(AC)] == '.ac' :
            f = float(line.split()[2])
    return f

def get_key(node_ref,value):
    ''' Gets the corresponding key for a value in the dictionary '''
    for key in node_ref.keys():
        if node_ref[key] == value :
            return key

def make_dict(circuit_def,elem):
    ''' Makes a dictionary for each component of the particular type of element '''
    e = elem
    volt_dict = {}
    volt_names = [element_info(line).words[0] for line in circuit_def if element_info(line).words[0][0].lower()== e]
    for index,name in enumerate(volt_names):
        volt_dict[name] = index
    return volt_dict

def finding_node(circuit_def,node_key,node_rep):
    ''' Finds the lines and position from/to of the given node '''
    index = []
    for i in range(len(circuit_def)):
        if(node_rep[element_info(circuit_def[i]).from_node] == node_key):
            index.append((i,1))
        elif(node_rep[element_info(circuit_def[i]).to_node] == node_key):
            index.append((i,2))
    return index

def dimensions(circuit_def):
    ''' Returns a tuple : number of nodes, number of voltage sources '''
    volt_ind = [i for i in range(len(circuit_def)) if circuit_def[i].split()[0][0] == 'V']
    k = len(volt_ind)
    n = len(dictionary_for_node(circuit_def))
    return n,k

def mod_matrix(circuit_def,f,node_key,node_rep,volt_dict,ind_dict,M,b):
    ''' Updates the M and b matrix for the given node '''

    index = finding_node(circuit_def,node_key,node_rep)
    n,k = dimensions(circuit_def)

    for j in index:
         #getting all the attributes of the element using the class definition
         elem = element_info(circuit_def[j[0]])
         element_name = circuit_def[j[0]].split()[0]

         if elem.name == 'resistor':
             if j[1]== 1 :
                 #node is the from_node
                 adj_key = node_rep[elem.to_node]
                 M[node_key,node_key] += 1/(elem.value)
                 M[node_key,adj_key] -= 1/(elem.value)
             if j[1] == 2 :
                 # node is the to_node
                 adj_key = node_rep[elem.from_node]
                 M[node_key,node_key] += 1/(elem.value)
                 M[node_key,adj_key] -= 1/(elem.value)

         if elem.name == 'capacitor':
             if j[1]== 1 :
                 #node is the from_node
                 adj_key = node_rep[elem.to_node]
                 M[node_key,node_key] += complex(0, 2*pi*f*(elem.value))
                 M[node_key,adj_key] -= complex(0, 2*pi*f*(elem.value))
             if j[1] == 2 :
                 # node is the to_node
                 adj_key = node_rep[elem.from_node]
                 M[node_key,node_key] += complex(0, 2*pi*f*(elem.value))
                 M[node_key,adj_key] -= complex(0, 2*pi*f*(elem.value))

         if elem.name == 'inductor':
             try:
                 if j[1]== 1 :
                     #node is the from_node
                     adj_key = node_rep[elem.to_node]
                     M[node_key,node_key] += complex(0, 1/(2*pi*f*(elem.value)))
                     M[node_key,adj_key] -= complex(0, 1/(2*pi*f*(elem.value)))
                 if j[1] == 2 :
                     # node is the to_node
                     adj_key = node_rep[elem.from_node]
                     M[node_key,node_key] += complex(0, 1/(2*pi*f*(elem.value)))
                     M[node_key,adj_key] -= complex(0, 1/(2*pi*f*(elem.value)))

             except ZeroDivisionError:

                 index = ind_dict[element_name]

                 if j[1]== 1:
                     M[node_key,n+k+index] += 1
                     M[n+k+index,node_key] -= 1
                     b[n+k+index] = 0

                 if j[1]== 2:
                     M[node_key,n+k+index] -= 1
                     M[n+k+index,node_key] += 1
                     b[n+k+index] = 0

         if elem.name == 'voltage source':

            index = volt_dict[element_name]

            if j[1]== 1:
                adj_key = node_rep[elem.to_node]
                M[node_key,n+index] += 1
                M[n+index,node_key] -= 1
                b[n+index] = elem.value

            if j[1] == 2 :
                adj_key = node_rep[elem.from_node]
                M[node_key,n+index] -= 1
                M[n+index,node_key] +=1
                b[n+index] = elem.value

         if elem.name== 'current source' :

            if j[1]== 1:
                b[node_key] -= elem.value
            if j[1] == 2 :
                b[node_key] += elem.value

         if (element_name[0]== 'E')|(element_name[0]== 'F')|(element_name[0]== 'G')|(element_name[0]== 'H'):

            print('This Program cannot solve dependent sources')
            exit()

"""
Checking with file input is given or not

"""
if len(argv)!=2:
    print('\nUsage: %s <input file>' % argv[0])
    exit()

circuit = '.circuit'
end = '.end'
print(argv[1])
try:
    with open(argv[1]) as file:
        lines = file.readlines()


        f = freq(lines)

        len_cir = len(circuit)
        len_end = len(end)
        start = -1
        finish = -1

        for k in lines:
            if circuit == k[:len_cir]:
                start = lines.index(k)          #start index
            elif end == k[:len_end]:
                finish = lines.index(k)         #end index
            k = k.split('#')[0]
        if start >= finish:                     #validating the code is write or not
           print('Invalid Spice Code')
           exit()

        circuit_def = lines[start+1:finish]

        for k in circuit_def:
            if element_info(k).from_node == element_info(k).to_node:
                print("From and To node of an element cannot be equal. Please redefine the circuit.")

        
        node_rep = dictionary_for_node(circuit_def)
        volt_dict = make_dict(circuit_def,'v')
        ind_dict = make_dict(circuit_def,'l')
        n,k = dimensions(circuit_def)
        dim = n+k
        M = zeros((dim,dim),dtype=complex)
        b = zeros(dim,dtype=complex)
        dc_flag = False
        if f == 0:
            dc_flag = True
            M = zeros((dim+len(ind_dict),dim+len(ind_dict)),dtype=complex)
            b = zeros(dim+len(ind_dict),dtype=complex)

        for i in range(len(node_rep)):
            mod_matrix(circuit_def,f,i,node_rep,volt_dict,ind_dict,M,b)
        M[0] = 0
        M[0,0] =1
        print('The M matrix is :\n',M)
        print('The b matrix is :\n',b)
        print('The node dictionary is :',node_rep)
        try:
            x = linalg.solve(M,b)
        except Exception:
            print('The incidence matrix cannot be inverted as it is singular. Please provide a valid circuit definition')
            sys.exit()
        for i in range(n):
            print("The voltage at node {} is {}".format(get_key(node_rep,i),x[i]))
        for j in range(k):
            print('The current through source {} is {}'.format(get_key(volt_dict,j),x[n+j]))
        if dc_flag:
            for i in range(len(ind_dict)):
                print('The current through inductor {} is {}'.format(get_key(ind_dict,i),x[n+k+i]))
        print('Voltage convention : From node of the voltage source is at a lower potential')

except IOError:
   print('Invalid file. Please make sure that the circuit definition block is well defined and all component value are in scientific notation.')





