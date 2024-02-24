import argparse
#Global Variables and Tunable Parameters
#Current Parameters seem to be the best functioning for 1b
#Finding minimizers Parameters
kmerLength_ = 35 
minimizerLength_ = 10 

#SNP and Gap Tolerance per alignment
SNPTolerance_ = 2
gapTolerance_ = 2

#Scoring metics for DP alignment
gapPenalty_ = -2
matchBonus_ = 3
mismatchPenalty_= -1

#Lower bound threshold of amount of coverage per SNP or gap
thresholdMatches_ = 4 

#reads in the fasta file and grabs the sequences from the file
def read_fasta(file_path):
    sequences = {}
    current_sequence = ""
    
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_sequence = line[1:]
                sequences[current_sequence] = ""
            else:
                sequences[current_sequence] += line
    
    return sequences
#returns a dict with names of sequences and sequences

#finds the minimizers of a given DNApiece
def findMinimizers(DNApiece,kmer = kmerLength_, minimizerLength = minimizerLength_):
    minimizerDict = {} # minimizer : [indexes]

    
    currentMinimizer = 'Z'*minimizerLength

    #slides the window of the kmer
    for kmerWindow in range(len(DNApiece)-kmer+1):
        window = DNApiece[kmerWindow:kmerWindow+kmer]
        sortedList = []

        #slides the window for the minimizers and adds all minimizers to a list
        for i in range(len(window)-minimizerLength +1):
            sortedList.append((hash(window[i:i+minimizerLength]),i))
        
        #sorts based on hash
        sortedList.sort()

        #if there is a new minimzer do the following 
        if sortedList[0][0] != currentMinimizer:
            currentMinimizer = sortedList[0][0] # update the new minimizer 
            index = sortedList[0][1] #index of the new 
            
            #if already in the dictionary add the new index found
            if currentMinimizer in minimizerDict:
                minimizerDict[currentMinimizer] += [index+kmerWindow]
            
            #else create an entry for it 
            else: 
                minimizerDict[currentMinimizer] = [index+kmerWindow]
                
    return minimizerDict

def DPSequenceAlignment(genome, read, match_score=matchBonus_, mismatch_penalty=mismatchPenalty_, gap_penalty=gapPenalty_):
    # Initialize the score matrix
    rows = len(genome) + 1
    cols = len(read) + 1

    dpMatrix = []

    for i in range(rows):
        dpMatrix.append([0] * cols)
 
    # Initialize the first row and column
    for i in range(rows):
        dpMatrix[i][0] = i * gap_penalty
    for j in range(cols):
        dpMatrix[0][j] = j * gap_penalty

    #we have to do -1 indexing instead of +1 indexing due to the direction we fill out the matrix
    #Fill in the score matrix
        
    for i in range(1, rows): #already filled out 1st column so start at 1
        for j in range(1, cols):
            #match them 
            if genome[i - 1] == read[j - 1]:
                match = dpMatrix[i - 1][j - 1] + match_score
            else:
                #mismatch them
                match = dpMatrix[i - 1][j - 1] + mismatch_penalty
            delete = dpMatrix[i - 1][j] + gap_penalty
            insert = dpMatrix[i][j - 1] + gap_penalty
            dpMatrix[i][j] = max(match, match, delete, insert)

    # Traceback to find the alignment
    genomeAlign = ""
    readAlignment = ""
    i = rows - 1
    j = cols - 1
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dpMatrix[i][j] == dpMatrix[i - 1][j - 1] + (match_score if genome[i - 1] == read[j - 1] else mismatch_penalty):
            genomeAlign = genome[i - 1] + genomeAlign
            readAlignment = read[j - 1] + readAlignment
            i -= 1
            j -= 1
        elif i > 0 and dpMatrix[i][j] == dpMatrix[i - 1][j] + gap_penalty:
            genomeAlign = genome[i - 1] + genomeAlign
            readAlignment = '-' + readAlignment
            i -= 1
        else:
            genomeAlign = '-' + genomeAlign
            readAlignment = read[j - 1] + readAlignment
            j -= 1

    return genomeAlign, readAlignment, dpMatrix[rows - 1][cols - 1]

#given a minimizer found in a read and the indexes of all of the places in the genome where that minimizer shows
#up, calcaluate their errors list and match score at each location and take the best match
def matchReadMinimizers(genome,readMinimizer,ReadIndexFound,sequence_read, genomeMinimizersIndexes):
    
    #list of errors returned
    listErrors = []

    for matchIndex in genomeMinimizersIndexes:     
        
        #initialize some variables
        mismatches = []
        SNPs = 0 
        gaps = 0

        #full sequence match 
        genomeSlice = genome[matchIndex-ReadIndexFound:len(sequence_read)+matchIndex-ReadIndexFound]
        
        #keeps track of where we are at in the genome
        trueAlignmentIndex = matchIndex-ReadIndexFound

        #alignment
        genomeAlignment, readAlignment, score = DPSequenceAlignment(genomeSlice, sequence_read)

        #check the alignment and get all of the errors in format
        for i in range(len(genomeAlignment)):

            #insertion
            if genomeAlignment[i]=='-':
                trueAlignmentIndex -=1
                gaps += 1
                mismatches.append((trueAlignmentIndex+i,'I',f' {readAlignment[i]}'))
            #deletion
            elif readAlignment[i]=='-':
                gaps += 1
                mismatches.append((trueAlignmentIndex+i,'D',f' {genomeAlignment[i]}'))
            #substitution
            elif genomeAlignment[i] != readAlignment[i]:
                mismatches.append((trueAlignmentIndex+i,'S',f' {genomeAlignment[i]} {readAlignment[i]}'))
                SNPs+=1  

        
        #if valid occurence add to the list of errors                
        if len(mismatches) > 0 and SNPs <= SNPTolerance_ and gaps <= gapTolerance_:
            listErrors.append((score,mismatches))  
        
    #sort the best errors
    listErrors.sort()

    if len(listErrors)> 0 and listErrors[-1][0]>0: 
        #return the best error
        
        return listErrors[-1][1],score
    else: 
        #print('no mismatches for this read')
        return [],0

def find_variants(reference_genome, reads):
    #reads in the file that we want
    file_path_reference = reference_genome
    file_path_reads = reads
    reference_genome = read_fasta(file_path_reference)
    reads_genomes = read_fasta(file_path_reads)

    #get the reference sequence
    for sequence_id, sequence_genome in reference_genome.items():
        #reference genome that I am comparing to 
        genome = sequence_genome 
        #all of the minimizers in the sequence
        genomeMinimizersDict = findMinimizers(sequence_genome)

    #dict that will be used to return later
    dictMistakes = {}

    #list to add things to the dict 
    mistakes = []

    #used for completion metric
    count = 0
    numberReads = len(reads_genomes.items())


    #looping through all of the reads
    for sequence_id, sequence_read in reads_genomes.items():
        #find the minimizer(s) for a given read 
        nestedDict = findMinimizers(sequence_read)
        count+=1
        
        #if minimizer exists
        if len(nestedDict)> 0:
            minimizerList = []

            #loop through all of the minimizers for that read
            for minimizer, index in nestedDict.items():
                #if the minimizer appears more than once deal with later
                if len(index) > 1:
                    print('Kmer-Minimizer appears more than once in the read, please fix me later!')
                
                #if the minimizer matches to a minimizer in the genome evaluate it
                elif minimizer in genomeMinimizersDict:

                    #index for the minimizer in th read
                    minimizerIndex = index[0]
                    
                    #match the minimizer to the location in the genome
                    errors,score = matchReadMinimizers(genome,minimizer,minimizerIndex,sequence_read, genomeMinimizersDict[minimizer])
                    minimizerList.append((score,errors))

            #find the best performing minimizer with the least read errors 
            if len(minimizerList)>0:
                minimizerList.sort()
                if len(minimizerList[-1][1])>0:
                    
                    #format the errors for submission
                    for i in minimizerList[-1][1]:
                        misString = f'>{i[1]}{i[0]}{i[2]}'
                        if misString not in mistakes:
                            dictMistakes[misString] = 1
                            mistakes.append(misString)
                        else: 
                            dictMistakes[misString]+= 1
            #metric of completion
            if count%1000 == 0:
                print('Number of reads left: ',numberReads-count)
    print('Finished')

    #bounding the amount of times an error must show up to be counted as a true error 
    finalMistakes = []
    for k in mistakes:
        if dictMistakes[k] > thresholdMatches_:
            finalMistakes.append(k)
    return finalMistakes

def main():
    parser = argparse.ArgumentParser(description='Find variants in genomic data.')
    parser.add_argument('--reference-genome', required=True, help='Path to the reference genome file')
    parser.add_argument('--reads', required=True, help='Path to the reads file')
    parser.add_argument('--output', help='Optional output file to write the mutations')

    args = parser.parse_args()

    # Call the function with the provided arguments
    mutations = find_variants(args.reference_genome, args.reads)

    if args.output:
        with open(args.output, "w") as file:
            for mutation in mutations:
                file.write(mutation + "\n")
    else:
        print("List of mutations:", mutations)

if __name__ == "__main__":
    main()