import numpy as np
import copy
import random

def HASTA(Q, DB, ktup):


    ### Create variables for function ###

    # a chararray for the words divided from Q
    dividedQ = np.chararray([len(Q) - ktup + 1, 1],unicode=True,itemsize=ktup)
    # an array of the nucleotides
    Nucleotides = ['A','C','G','T']
    # a chararray for the nucleotides in the hash table
    hash = np.chararray([len(Nucleotides) ** ktup, 1],unicode=True,itemsize=ktup)
    # an empty list for the indices in the hash table
    hashIndices = []
    # an array for the calculation 'Pos(Q) – Pos(DB)'
    differenceMat = np.array([],dtype=int) 



    ### Divide Q into words of length “ktup” ###

    # Run 'len(Q) - ktup' times to divide Q properly to string in 'ktup' length
    for index in range(0, len(Q) - ktup):

        # Reset string
        current_str = ''

        # Add nucleotides to new string
        for add in range(0 ,ktup):
            current_str = current_str + Q[ index + add ]

        # Assign string to list
        dividedQ[index] = current_str



    ### Create strings for Hash table ###

    # Define a counter
    counter = 0

    # Run 'ktup' times in order to get correct hash table
    for x in range(0,ktup):
        
        # define counter for the index of the nucleotides
        nucIndex = 0

        # Fill the hash table with strings of nucleotides
        for index in range(0, len(hash)):
            
            # Create and assign string to list
            hash[index] = hash[index] + Nucleotides[nucIndex]
            counter += 1

            # First round of assignments
            if x == 0:

                # Reset counter since length of 'Nucleotides' list is 4
                if counter == 4:
                    counter = 0
                    nucIndex += 1
            
            # Second round of assignments
            if x == 1:

                nucIndex += 1
                if counter == 4:
                    counter = 0
                    nucIndex = 0
                    



    ### Build Hash talbe for DB ###

    # Define the index as the legnth of the hash table
    for index in range(0, len(hash)):

        # Define variable for the condition of the while loop
        stop = 0
        # Assign dataBase to current string
        strDB = DB
        # Define list for current indices
        currentLst = []

        while stop == 0:
            
            # Find if current string appears in dataBase
            currentRes = strDB.find(hash[index,0])

            # condition when current string DOES appears in dataBase
            if currentRes != -1:
                
                # Cut the string for next iter in while loop
                strDB = strDB[currentRes + 1:]
                # Assign index to list
                currentLst.append(currentRes)

            # condition when current string does NOT appears in dataBase
            else:

                # Stops the while loop
                stop = 1   

        # Add current indices list to hash list
        hashIndices.append(currentLst)                



    ### Convert hashIndices from List to numpy array ###

    # Find the longest list in our list (of lists)
    maxLength = len(max(hashIndices))

    # 'None' padding for all lists
    for row in hashIndices:

        # When the list is shorter than longest list
        while len(row) < maxLength:

            # padding with 'None'
            row.append(None)
    
    # Assign balanced list to numpy array
    hashIndices = np.array(hashIndices)



    ### POS(Q) – POS(DB) ###

    # Run on the index of dividedQ
    for indexQ in range(0, len(dividedQ)):
        
        # Find the index in hash table, where the string fron dividedQ appears
        currentRes = dividedQ[indexQ].find(hash[:, 0])
        currentRes = np.where(currentRes == 0)

        # Find the indices in DB string (from 'hashIndices')
        DBIndices = hashIndices[currentRes].copy()

        # When the string DOES appear in DB
        if DBIndices.size > 0:

            # Goes over the indices in DB
            for indexDB in DBIndices[0]:

                # When the value of the index is NOT 'None'
                if indexDB != None:
                    
                    # Add 'Pos(Q) – Pos(DB)' to the differnceMat
                    differenceMat = np.append(differenceMat, indexQ - int(indexDB))



    ### Find indices for best alignment ###

    # Count how many times each index appears to find hotspots
    counts = np.bincount(abs(differenceMat))

    # Check if any of the incides appears more the once
    countsCheck = np.any(counts > 1)

    # When at least one index appears more than once
    if countsCheck:

        # The index for the beginning is chosen based on the most frequent index
        indexBeg = np.bincount(abs(differenceMat)).argmax()

    # When no index appears more than once 
    else:

        # The index for the beginning is chosen randomly from 'differenceMat', since there are no hotspots
        randIndex = random.randint(0,len(differenceMat)-1)
        indexBeg = abs(differenceMat[randIndex])

    # Calculate index for the end of the alignment
    indexEnd = indexBeg + len(Q) - 1

    # Print the result
    print('The query best aligned from %s to %d' %(indexBeg, indexEnd))





# Define query
query = 'CATGG'

# Define database
dataBase = 'ATTGCAGTGGTCCC'

# Define ktup
ktup = 2

# Run the function
HASTA(query, dataBase, ktup)