# Python script to convert the ID's of one of the photons in the Higgs--> photon photon decay to (eg.) a neutralino in a LHE file
# © Jamie Bamber 2017

# Use: ./LHE_photon_neutralino_v2.py "<input_filename>" "<output_filename>"

import random
import sys

# define input file and neutralino (or equiv. PDG ID) & random() seed

new_particle_PDG = '1000022'
input_file = sys.argv[1]
output_file = sys.argv[2]
random.seed(324)

with open(input_file, "r") as file1:
	raw_text = file1.read()
# split into nested lists
lists = raw_text.split("\n")
for i in range(len(lists)):
	lists[i] = lists[i].split(" ")
outlists = lists # create output lists
#

#print(lists[12])

# iterate through the rows, checking for Higgs particles
for i in range(len(lists)):
	# find photon pairs
	if len(lists[i])>7 and len(lists[i+1])>7:  
		if lists[i][7]=='22' and lists[i+1][7]=='22':
			# change PDG IDs
			R = random.choice([0, 1])
			outlists[i+R][7] = new_particle_PDG

# remake outlists into a text string
for i in range(len(outlists)):
	outlists[i] = " ".join(outlists[i])
out_text = "\n".join(outlists)

with open(output_file, 'w') as file2:
	file2.write(out_text)
	



				



	
	