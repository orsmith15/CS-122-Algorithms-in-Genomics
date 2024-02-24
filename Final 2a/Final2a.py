import argparse
#Global Variables and Tunable Parameters
#Current Parameters seem to be the best functioning for 1b
# Window Parameters for reads 
window_ = 20

#SNP and Gap Tolerance per alignment
SNPTolerance_ = 3
gapTolerance_ = 2

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


   # depending on how many counts there are add the pieces into the kmer list that will be used 
    # to make the debruijn graph
    # if there is a higher count of a specific kmer then it is probably a repeat
def countKmers(spectrum):
    dictReadNames = {}
    kmerList = []
    # go through the spectrum
    for read_name,read_sequence in spectrum.items():

        kmerList.append(read_sequence)
        if read_sequence in dictReadNames: 
            dictReadNames[read_sequence] += [read_name]
        else: 
            dictReadNames[read_sequence] = [read_name]

    return kmerList,dictReadNames

# creates my debruijn graph
def debruijn(kmerList):
    deBruijn = {}
    for kmer in kmerList: 
        suffix = kmer[1:]
        prefix = kmer[:-1]
        if prefix not in deBruijn:
            deBruijn[prefix] = [suffix]
        else: 
            deBruijn[prefix] += [suffix]
    print('Size of deBruijn graph',len(deBruijn))
    return deBruijn 
# counts in degree and out degree to find the unbalanced nodes
def findStartingNode(deBruijn):
    outDegree = {}
    inDegree = {}
    
    for lmer in deBruijn:
        outDegree[lmer] = len(deBruijn[lmer])
        for out in deBruijn[lmer]:
           
            if out in inDegree:
                inDegree[out] += 1
            else:
                inDegree[out] = 1
    lmerList = set(list(outDegree.keys()) + list(inDegree.keys()))

    # finds all of the unbalanced nodes
    
    unbalancedNodes = []
    for lmerCheck in lmerList: 
        if lmerCheck not in inDegree:
            unbalancedNodes.append(lmerCheck)
        elif lmerCheck not in outDegree: 
            unbalancedNodes.append(lmerCheck)
        elif outDegree[lmerCheck] != inDegree[lmerCheck]:
            unbalancedNodes.append(lmerCheck)
    potentialStart = []
    if len(unbalancedNodes)> 0:
        for i in unbalancedNodes: 
            if i not in outDegree: 
                continue 
            elif i not in inDegree: 
                potentialStart.append(i)
            elif outDegree[i]>inDegree[i]:
                potentialStart.append(i)
    if len(potentialStart) != 1:
        print('WARNING: Number of potential starting nodes: ', len(potentialStart), ' List:',potentialStart)
    return potentialStart       
#returns a dict with names of sequences and hash of sequences
def hashGenome(genome,window = window_):
    hashDict = {}
    for i in range(len(genome)-window+1):
        kmer = genome[i:i+window]
        if kmer in hashDict:
            hashDict[kmer] += [i]
        else: 
            hashDict[kmer] = [i]
    return hashDict

# does a random walk through the eulerian until you hit a dead end
def randomWalk(potentialStart,deBruijn):
    path = []
    node = potentialStart
    isCycle = False
    while node in deBruijn and len(deBruijn[node]) != 0: 
        path.append(node)
        current = node
        node = deBruijn[node][0]
        deBruijn[current] = deBruijn[current][1:]
    
    path.append(node)
    if path[0] == path[-1]:
        isCycle = True
    return path, isCycle

# finds the eulerian path through the graph given a potential start and the debruijn graph
def findEulerianPath(potentialStart,deBruijn):

    if len(potentialStart) < 2:
        
        if len(potentialStart) == 0:
            startNode = list(deBruijn.keys())[0]
        else:
            startNode = potentialStart[0] 
            
        path,isCycle = randomWalk(startNode,deBruijn)

        recursiveBreak = False
        # if there are values left that means there is a cycle, need to run it back again 
        while any(deBruijn.values()) and not recursiveBreak:
            print('GOING TO SECOND RANDOM WALK')
            recursiveBreak = True
            for index,node in enumerate(path): 
                if node in deBruijn: 
                    if len(deBruijn[node]) != 0:
                        ifNotCycleWarranty = deBruijn.copy()
                        newCycle, isCycle = randomWalk(node,ifNotCycleWarranty)
                        if not isCycle: 
                            print('new path is not a cycle: ', newCycle, 'index of occurence: ',index)
                        else: 
                            deBruijn = ifNotCycleWarranty.copy()
                            path = path[:index+1] + newCycle[1:] + path[index+1:]
                            recursiveBreak = False
    if len(potentialStart)> 1 and any(deBruijn.values()): 
        print('MORE THAN 1 valid start and not all edges used up, please come back and fix me')   

    return path 

# finds the genome sequence
def findGenomeSequence(path):
    genome = path[0]
    for kmer in path[1:-1]:
        genome += kmer[-1]
    
    # adds the last letter on 
    genome += path[-1][-1]

    return genome

# compact version to call all of the functions needed to find the genome
def findGenome(spectrum):
    kmerList,dictReadNames = countKmers(spectrum)
    deBruijn = debruijn(kmerList)
    potentialStartList = findStartingNode(deBruijn)
    pathList = findEulerianPath(potentialStartList,deBruijn)
    genome = findGenomeSequence(pathList)
    return genome,pathList,dictReadNames

# finds the mismatches stringwise
def misMatches(genomeMatchIndex,sequence_read,genome):
    
    # mismatches to be returned in case of SNP/gap finding
    misMatchList = []

    # keeps track of the gap and mismatch tolerance
    gaps = 0
    mismatches = 0

    # where we are at in the sequences
    i = 0
    j = 0
    
    # while we are not at the end of the sequence compare the sequences
    while i < len(sequence_read) and j < len(genome):

        #true index of the genome
        positionRefGenome = genomeMatchIndex + j
        baseRead = sequence_read[i]
        baseGenome = genome[j]
        
        if baseRead != baseGenome:
            # if not equal first check if it is a mismatch
            if mismatches < SNPTolerance_:
                # avoids the chances of length of genome and read not being the same
                length = min(len(genome[j+1:]),len(sequence_read[i+1:]))
                misMatch,MCount = checkIndel(sequence_read[i+1:i+1+length],genome[j+1:j+1+length],mismatches,typeCheck = 'mismatch')

                # if we consider this one base change to be a mismatch and the rest is perfect match, our work is done
                if misMatch:
                  mistake = (positionRefGenome,'S',f' {baseGenome} {baseRead}')
                  misMatchList.append(mistake)
                  mismatches += 1

                  return True,misMatchList
                
            # if the one base change doesn't work then check if this is an indel
            if gaps < gapTolerance_:
                insertion = False
                deletion = False

                # check if there is an insertion
                if i+1 < len(sequence_read):
                    length =min(len(sequence_read[i+1:]),len(genome[j:]))
                    insertion,ICount = checkIndel(sequence_read[i+1:i+1+length],genome[j:j+length],mismatches)

                # check if there is a deletion 
                if j+1 <len(sequence_read):
                    length = min(len(genome[j+1:]),len(sequence_read[i:]))
                    deletion,DCount = checkIndel(sequence_read[i:i+length],genome[j+1:j+1+length],mismatches)

                # if both are true find the best one
                if deletion and insertion:
                    
                    # if neither one is better return false
                    if DCount == ICount: #can add in randomness of what I choose it to be marked as
                      return False,misMatchList # TODO SUS <-- come back to, this might mess up everything
                      # check if it is a mismatch
                    elif DCount > ICount:
                        deletion = False
                    elif ICount>DCount:
                        insertion = False

                # count as deletion and adjust pointers
                if deletion:
                    misMatchList.append((positionRefGenome,'D',f' {baseGenome}'))
                    gaps+=1
                    j+=1
                    continue
                
                # count as insertion and adjust pointers
                elif insertion:
                    misMatchList.append((positionRefGenome,'I',f' {baseRead}'))
                    gaps+=1
                    i+=1
                    continue

                #if both of them arn't true then consider it a mismatch and check again. This would be the case that there is more than 1 indel/mismatch in the sequence
                else: #if deletion and insertion both false consider that it might be a mismatch
                    # we would expect that if this isn't true then the mismatch tolerance will at some point go over and they won't match
                    if mismatches < SNPTolerance_:
                      mistake = (positionRefGenome,'S',f' {baseGenome} {baseRead}')
                      misMatchList.append(mistake)
                      mismatches += 1
                      i+=1
                      j+=1
                    else:
                      return False,misMatchList
            # if too many gaps and mismatches return
            else:
                return False,misMatchList
            
        #if equal advance the two pointers
        else:
            i+=1
            j+=1

    return True,misMatchList

#check if the base pairs match enough
def checkIndel(read,genome,mismatches,typeCheck='I/D'):
  # if indel then check that considering it an indel will be ok with the SNPtolerance
  if typeCheck == 'I/D':
    count = 0
    count = mismatches
    for i in range(len(read)):
        if read[i]!=genome[i]:
            count += 1
        if count < SNPTolerance_:
            continue
        else:
            return False, 0
    return True, count
  else:
    # if mismatch then it needs to be a perfect match 
    for i in range(len(read)):
        if read[i]!=genome[i]:
          return False, 0
    return True,0

# given a match to the split of the read check if it is a valid match spot
def NEWmatchReadMinimizers(genome,readMinimizer,ReadIndexFound,sequence_read, genomeMinimizersIndexes):
    listIndexes = []
    #for all of the places in which there is a double match in the genome calculate the spots we need to look for
    for matchIndex in genomeMinimizersIndexes:     
        
        #initialize some variables
        mismatches = []

        #full sequence match 
        genomeSlice = genome[matchIndex-ReadIndexFound:len(sequence_read)+matchIndex-ReadIndexFound]
        trueAlignmentIndex = matchIndex-ReadIndexFound

        valid,mismatches= misMatches(trueAlignmentIndex,sequence_read,genomeSlice)
        
        #if valid occurence add to the dict of Errors                
        if valid:
            listIndexes.append(matchIndex)
    return listIndexes  

def findAlignment(genome,spectrum):
    hashDict = hashGenome(genome)
    indexDict = {}
    numberReads = len(spectrum)
    count = 0
    
    for sequence_read, junk in spectrum.items():
        count+=1
        
        #nestedDict = [sequence_read[:window_],sequence_read[window_:window_*2],sequence_read[window_*2:window_*3]]
        nestedDict = [sequence_read]
        if len(nestedDict)> 0:
            #loop through all of the hashes for that read
            for index,minimizer in enumerate(nestedDict):
                # if the minimizer matches to a minimizer in the genome evaluate it
                if minimizer in hashDict:
                    
                    listValid = NEWmatchReadMinimizers(genome,minimizer,index*window_,sequence_read, hashDict[minimizer])
                    # I want this to return the number of indexes where this read matches to 
                    # if valid mark it as such in the metaGenomics hash and continue to the next genome
                    if len(listValid)>0:
                        #print('valid read with this many indicies',len(listValid))
                        indexDict[sequence_read] = listValid
                        break
                                     
        if count%1000 == 0:
            print('Genomes Left for hash:',numberReads-count)
    return indexDict



def find_variants(spectrum):
    spectrum = read_fasta(spectrum)
    genome,pathList,dictSequenceReads = findGenome(spectrum) # kmerSequence: name of reads
    
    indexDict = findAlignment(genome,dictSequenceReads)
    
    dictAnswer = {}
    for sequence_read,listValid in indexDict.items():
        if sequence_read in dictSequenceReads:
            if len(listValid) == len(dictSequenceReads[sequence_read]):
                for index,i in enumerate(listValid):
                    if i not in dictAnswer:
                        #print(str(dictSequenceReads[sequence_read]))
                        dictAnswer[i] = [dictSequenceReads[sequence_read][index]]
                    else: 
                        print('collision of indicies')
                        #dictAnswer[i] += [dictSequenceReads[sequence_read][index]]
            else: 
                print('ERROR mismatch of indexes read could be found and number if identical reads')
    answer = sorted(dictAnswer)
    answerList = []
    for i in answer: 
        if i in dictAnswer: 
            answerList.extend(dictAnswer[i])

    return answerList

def main():
    parser = argparse.ArgumentParser(description='Find variants in genomic data.')
    parser.add_argument('--spectrum', required=True, help='Path to the spectrum of the genome')
    parser.add_argument('--output', help='Optional output file to write the order of reads')

    args = parser.parse_args()

    # Call the function with the provided arguments
    mutations = find_variants(args.spectrum)

    if args.output:
        with open(args.output, "w") as file:
            for mutation in mutations:
                file.write(mutation + "\n")
    else:
        print("List of reads and genomes:", mutations)

if __name__ == "__main__":
    main()