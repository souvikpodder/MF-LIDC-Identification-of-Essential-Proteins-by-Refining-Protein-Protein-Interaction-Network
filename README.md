# MF-LIDC-Identification-of-Essential-Proteins-by-Refining-Protein-Protein-Interaction-Network
Source code

	run mainDriver() function:
		input: 
				none
		output:
				none
	call nodeWeight(data) function:
		input: 	
				data is the name of the data file(the ppin data file).
		output:	
				A list of node weight of PIN edges.
 
		call graphBuild(resultPath,dataPath) function:
		input:	
				1. resultPath is name of the main graph file (it will be created later.) .
				2. dataPath is name of the data file (the ppin data file).
		output: 	
				1. A dictionary of the PIN.
				2. Above dictionary written into a file with the given name written with json.

		call uniqueProtein(filePath) function:
		input:
				1. filePath is the path of the reduced edge file.
		output:
				2. a list of unique proteins.

		call nodeDelete(alpha,SD,proteinList,proteinDegList,graph,k) function 3 times with K value (1,2,3):
		input:	
				1. k is an integer for calculating threshold.
				2. alpha is the mean of edge weights.
				3. SD is the standard deviation of edge weights.
				4. proteinList is the list of unique Protein.
				5. proteinDegList is a list with protein degs as values.
		output:	
				1. Two files with reduced edges edgesAfterNodeReduction(k).txt written with json and edgesAfterNodeReduction(k).csv written with pandas.
				2. A dictionary(Graph) written into reducedNodeGraph(k).txt with json.



			call edgeWeight(k,graph) function:
			input:	
					1. k is an integer for calculating threshold.
					2. A file of reduced edge(k) written with json.
					3. graph is a dictionary of proteins as a key and connected proteins as values.
			output:
					1. A TXT file with reduced node(k).
					2. A TXT file for reduced node graph(k) written with json.

				call uniqueProteinL(filePath) function:
				input:
					1. filePath is the path of the reduced edge file.
				output:
					2. Returns a list of unique proteins.

				call nonEssential(k,data)
					input:
							1. k is an integer
							2. data is the the ppin data file

				call LIDCDriver(k,graph) function:
					input: 	none
					output:	none

					call IDC(k,graph) function:
						input: 
								1. k is an integer.
								2. A excel file of protein complex.
								3. A TXT file of reduced node(protein).
						output:
								1. A file of IDC values of each protein written with json.
								2. Returns a dictionary of each protein with their IDC value.

					call LID(k,graph) function:
						input:
								1. k is an integer.
								2. The reduced graph.
						output:
								1. A file of LID values of each protein written with json.
								2. Returns a sorted dictionary of protein and LID value pair.
					call LIDC(LIDList,IDCList,k) function:
						input:
								1. LIDList is a dictionary of proteins with their LID values.
								2. IDCList is a dictionary of proteins with their IDC values.
								3. k is an integer.
						output:
								1. A file with LIDC values.
