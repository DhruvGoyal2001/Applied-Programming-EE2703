''' EE2703 Apllied Programming Lab - 2021 Assignment 1
	Dhruv Goyal, EE19B077

	This program reads tokens from a netlist file and stores them in a dictionary and prints
	the tokens in reversed order without any trailing comments

	Usage: python3 assign1.py ckt1.netlist
'''

# importing functions from sys module 

from sys import argv, exit 

# defining globally so that easy to change later

Circuit = '.circuit'
End = '.end'

# checking if valid number of arguments is given

if len(argv) !=2:		
	print("Usage: %s <filename.netlist>" %argv[0])
	exit()

#defining empty dicitonary
Tokens = {} 

#reading input from given file
try:
	with open(argv[1]) as f:  							
		data = f.readlines()
		start = 0
		end = -1
		for line in data:
			if Circuit == line[:len(Circuit)]:			# obtaining start index
				start = data.index(line)
			elif End == line[:len(End)]:				# obtaining end index
				end = data.index(line)

		if start >= end:								#Checking for valid circuit definition
			print("Invalid Spice Code!")
			exit(0)

		for k in data[end-1:start:-1]:
			k = k.split('#')[0].split()			#removing any trailing comment and splitting

			Tokens.update({k[0]: k[1:]})		#updating dictionary with key as corresponding name of element

			k.reverse()							#reversing order of each token
			k = " ".join(k)						#joining with spaces in between

			print(k)							#printing out the values

	
except IOError:									#if unable to load file
	print("Invalid File")
	exit()

"""
If we had to write the same code in C, it would have been comparitively huge since we would
need to define all the fucntions which we have conviniently	imported or used here
"""
	