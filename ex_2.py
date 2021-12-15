import numpy as np
import matplotlib.pyplot as plt
import random
import copy


def seqAlignment(A, B):
    


 # This function gets two input arguments - the input is two strings

    len_A = len(A) # Find length of A
    len_B = len(B) # Find length of B
    fMatrix = np.zeros([len(A)+1, len(B)+1]) # Create F matrix
    traceBack = np.chararray([len(A), len(B)],unicode=True)
    indexMatrix = np.zeros([len(A), len(B)])
    structureMatrix = np.chararray([len(A)+2, len(B)+2],unicode=True,itemsize=2)
    newA = ''
    newB = ''
    d = 2 # Gap penalty



    ### Create structure matrix ###

    structureMatrix[:][:] = ''
    structureMatrix[1][0] = '-'
    structureMatrix[0][1] = '-'

    # Add A string
    for i in range(2,len_A+2):
        structureMatrix[i][0] = A[i-2]
    
    # Add string B
    for i in range(2,len_B+2):
        structureMatrix[0][i] = B[i-2]



    ### Create F matrix ###

    # Create gap penalty score for first row
    for i in range(0, len_A+1):
        fMatrix[i][0] = -i*d

    # Create gap penalty score for first col
    for i in range(0, len_B+1):
        fMatrix[0][i] = -i*d

    # Find scores for each match
    for i in range(1,len_A+1):
        for j in range(1,len_B+1):

            # score for match or mismatch
            if A[i-1] == B[j-1]:
                score = 2
            else:
                score = -1

            # find score for [i-1][j-1], [i][j-1] and [i-1][j]
            scoreVec = np.array([fMatrix[i-1][j-1]+score, fMatrix[i][j-1]-d, fMatrix[i-1][j]-d])
            # Assign best score to fMatrix
            fMatrix[i][j] = max(scoreVec)
            # Assign indices in list for traceBackk
            indexMatrix[i-1][j-1] = np.argmax(scoreVec)
    


    ### Create Trace back matrix ###

    i = len(A)-1
    j = len(B)-1
    stop = 0
    while stop == 0:

        if i >= 0 and j >= 0:

            if indexMatrix[i][j] == 0:
                traceBack[i][j]= 'D'
                newA = newA + A[i]
                newB = newB + B[j]
                i -= 1
                j -= 1
            elif indexMatrix[i][j] == 1:
                traceBack[i][j] = 'L'
                newA = newA + '-'
                newB = newB + B[j]
                j -= 1
            else:
                traceBack[i][j] = 'U'
                newA = newA + A[i]
                newB = newB + '-'
                i -= 1

        else:

            stop = 1


        


    newA = newA[::-1]
    newB = newB[::-1]


    ### Create output variables ###

    # find best score
    alignmentScore = fMatrix[len(A)][len(B)]

    # Create fMatrix 
    fMatrixOutput = copy.deepcopy(structureMatrix)
    fMatrixOutput[1:,1:] = copy.deepcopy(fMatrix[:][:])
    # Create traceBack output
    traceBackOutput = copy.deepcopy(structureMatrix)
    traceBackOutput[2:,2:] = copy.deepcopy(traceBack[:][:])

    # return all variables 
    return[fMatrixOutput, traceBackOutput, newA, newB, alignmentScore]


A = 'AGGCT'
B = 'AGCCAT'          
[fMatrix, traceBack, newA, newB, alignmentScore] = seqAlignment(A, B)
print(fMatrix, traceBack, newA, newB, alignmentScore)

A = 'GAATTCAGTTA'
B = 'GGATCGTTA'          
[fMatrix, traceBack, newA, newB, alignmentScore] = seqAlignment(A, B)
print(fMatrix, traceBack, newA, newB, alignmentScore)



### Find p value ###
Nucleotides = ['A', 'T', 'C', 'G']
scoreVec = np.zeros([15000])
iterNum = 15000
strings = [A, B]

# create mutated sequences and calculate their score
for iteration in range(0,iterNum):
    # create variable for new sequences
    mutatedStrings = [ [] for x in range(0,len(strings)) ]
    # Define counter 
    counter = 0

    for string in strings:
        
        # Define randomly an index for the mutation
        randomIndex = random.randint(0,len(string)-1)
        # Define randomly type of mutation
        randonMutation = random.randint(0,2)
        # Define randomly type of nucleotide
        nucNum = random.randint(0,3)
        
        # For mutation in first nucleotide
        if randomIndex == 0: 

            if randonMutation == 0: # deletion mutation
                currentString = copy.deepcopy(string[randomIndex+1:])
            elif randonMutation == 1: # Addition mutation
                currentString = copy.deepcopy(Nucleotides[nucNum] + string[randomIndex:])
            else: # Replace mutation
                currentString = copy.deepcopy(Nucleotides[nucNum] + string[randomIndex+1:])

        # For mutation in last nucleotide
        elif randomIndex == len(string): 

            if randonMutation == 0: # deletion mutation
                currentString = copy.deepcopy(string[:randomIndex])
            elif randonMutation == 1: # Addition mutation
                currentString = copy.deepcopy(string[:] + Nucleotides[nucNum])
            else: # Replace mutation
                currentString = copy.deepcopy(string[:randomIndex] + Nucleotides[nucNum])

        # For mutation in nucleotide in the middle
        else: 

            if randonMutation == 0: # deletion mutation
                currentString = copy.deepcopy(string[:randomIndex] + string[randomIndex+1:])
            elif randonMutation == 1: # Addition mutation
                currentString = copy.deepcopy(string[:randomIndex] + Nucleotides[nucNum] + string[randomIndex:])
            else: # Replace mutation
                currentString = copy.deepcopy(string[:randomIndex] + Nucleotides[nucNum] + string[randomIndex+1:])

        # Assign new mutated string to variable
        mutatedStrings[counter] = currentString
        # Update counter
        counter += 1

    # Run mutated sequences in seqAlignment function to calculate score
    funcOutput = seqAlignment(mutatedStrings[0], mutatedStrings[1])
    # Assign new score to score vector
    scoreVec[iteration] = funcOutput[4]

# Calculate and print p value
pValue = sum(scoreVec > alignmentScore)/iterNum
print(pValue)



# Plot the score of the mutated sequence in a histogram

n_bins = 15
x = scoreVec
      
plt.hist(x, n_bins)
plt.xlabel('Score')
plt.ylabel('# of sequence pair')
plt.title('Permutated sequence pairs - score histogram')


            

            
            
            



