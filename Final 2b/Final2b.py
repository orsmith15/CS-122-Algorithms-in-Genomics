import argparse
import statistics
#Global Variables and Tunable Parameters
#Current Parameters seem to be the best functioning for 1b
# Window Parameters for reads 
window_ = 12

#SNP and Gap Tolerance per alignment
SNPTolerance_ = 3
gapTolerance_ = 3

# threshold for have 2 of the same spectrum 
thresholdFor2_ = 2.5

# overlap size for disjoint eulerian paths
overlapSize_ = 10

# lower bound of frequency filter
filter_ = 5

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

def makeSpectrum(reads,window = window_):
    spectrum = {}
    for read_id,read_sequence in reads.items():
        for i in range(len(read_sequence)-window_+1):
            kmer = read_sequence[i:i+window]

            if kmer in spectrum:
                spectrum[kmer] += [read_id]
            else:
                spectrum[kmer] = [read_id]
    return spectrum # dict = kmer_sequence: read_id

# counting all of the kmers from the spectrum
def countKmers(spectrum,threshold):
    # counting all of the kmers
    countList = []
    
    for kmer_sequence, reads in spectrum.items():
        if len(reads)> threshold: 
            countList.append(len(reads))
    
    mean = statistics.mean(countList)
    maxi = max(countList)
    standard_deviation = statistics.stdev(countList)

    kmerList = []
    dictReadNames = {}
    numDuplicates = []
    for read_sequence,read_name in spectrum.items():
        #print(len(read_name))
        if len(read_name)>threshold:
            kmerList.append(read_sequence)
            if read_sequence in dictReadNames:
                dictReadNames[read_sequence] += read_name
            else:
                dictReadNames[read_sequence] = read_name
    return kmerList,dictReadNames,mean, standard_deviation

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
    print('Size of deBruijn graph, do you really want to make a copy of it?',len(deBruijn))
    return deBruijn

# counts in degree and out degree to find the unbalanced nodes
def findStartingNode(deBruijn, count = True):
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
            unbalancedNodes.append((0,lmerCheck))
        elif lmerCheck not in outDegree:
            continue
        elif outDegree[lmerCheck] != inDegree[lmerCheck] and outDegree[lmerCheck]>inDegree[lmerCheck]:
            unbalancedNodes.append((inDegree[lmerCheck],lmerCheck))

    # check this logic
    # to be a potential start node it must have more out arrows than in arrows
            # TODO potentially just choose the node that has no indegree/sort the list

    potentialStart = []
    if len(unbalancedNodes)> 0:
        unbalancedNodes.sort()
        for i in unbalancedNodes:
          potentialStart.append(i[1])


    print('Number of potential items: ', len(potentialStart))
    return potentialStart,outDegree,inDegree

def findStartingNode2(outDegree,inDegree,lmerList):
    
    unbalancedNodes = []
    for lmerCheck in lmerList:
        if lmerCheck not in inDegree:
            unbalancedNodes.append((0,lmerCheck))
        elif lmerCheck not in outDegree:
            continue
        elif outDegree[lmerCheck] != inDegree[lmerCheck] and outDegree[lmerCheck]>inDegree[lmerCheck]:
            unbalancedNodes.append((inDegree[lmerCheck],lmerCheck))
    
    potentialStart = []
    if len(unbalancedNodes)> 0:
        unbalancedNodes.sort()
        for i in unbalancedNodes:
          potentialStart.append(i[1])

    
    return potentialStart

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
        if len(deBruijn[current]) == 0: 
            deBruijn.pop(current, None)

    path.append(node)
    if path[0] == path[-1]:
        isCycle = True
    return path, isCycle

def findEulerianPath(startNode,deBruijn):

    path, isCycle=randomWalk(startNode,deBruijn)

    recursiveBreak = False
    while any(deBruijn.values()) and not recursiveBreak:
        print('GOING TO SECOND RANDOM WALK')
        recursiveBreak = True
        for index,node in enumerate(path):
            if node in deBruijn and len(deBruijn[node]) != 0:
                ifNotCycleWarranty = deBruijn.copy()
                newCycle, isCycle = randomWalk(node,ifNotCycleWarranty)
                if not isCycle:
                    continue
                else:
                    deBruijn = ifNotCycleWarranty.copy()
                    path = path[:index+1] + newCycle[1:] + path[index+1:]
                    recursiveBreak = False

    return path,deBruijn

def findGenomeSequence(path):
    genome = path[0]
    for kmer in path[1:-1]:
        genome += kmer[-1]
    genome += path[-1][-1]
    return genome

def findOverlap(seq1,seq2):
    minim = min(len(seq1),len(seq2))
    max = 0
    for i in range(1,minim+1):
        if seq1[len(seq1)-i:] == seq2[:i]:
            max = len(seq2[:i])        
    return max

def sortPairs(pair_list):
    sorted_list = []
    
    while pair_list:
        for pair in pair_list:
            if len(sorted_list) == 0:
                sorted_list.append(pair)
                pair_list.remove(pair)
                break
            elif sorted_list[-1][2] == pair[1]:
                sorted_list.append(pair)
                pair_list.remove(pair)
                break
            elif pair[2] ==sorted_list[0][1]:
                sorted_list.insert(0,pair)
                pair_list.remove(pair)
                break
        else: 
            break
    return sorted_list,pair_list

def overlapsDict(paths):
    overlapDict ={}
    overlapList = []
    print('PATHS',paths)
    for i in paths: 
        
        #overlapDict[i] = {}
        for j in paths:
            if i ==j: 
                continue
                #overlapDict[i][j] = 0
            else:
                #overlapDict[i][j] = findOverlap(i,j)
                over = findOverlap(i,j)
                if over > overlapSize_: 
                    overlapList.append((str(over),i,j))
    print('overlapList',overlapList)
    if len(overlapList) != 0:
        first = []
        last = []
        pairList = []
        overlapList.sort(reverse=True)
        for i in overlapList: 
            if i[1] not in first and i[2] not in last and sorted(i) not in pairList: 
                pairList.append(i)
        sortLists = []
        
        while any(pairList):
            sortfinalOverlap,pairList = sortPairs(pairList)
            if len(pairList)>0:
                sortLists.append(sortfinalOverlap)

        genome = ''
        for sortfinalOverlap in sortLists:
            for pair in sortfinalOverlap:
                seq = ''
                if seq == '':
                    seq = pair[1] + pair[2][int(pair[0]):]
                else:
                    seq += pair[2][int(pair[0]):]
            genome += seq
    else:
        genome = "".join(paths)

    return genome

def findGenome(reads):
    # make the spectrum 
    spectrumPaired=makeSpectrum(reads)

    kmerList,dictReadNames,mean, standard_deviation=countKmers(spectrumPaired,filter_)#5)
    deBruijn = debruijn(kmerList)

    pathList = []
    path = [1]
    starts = []
    # find starting node
    potentialStartList,outDegree,inDegree = findStartingNode(deBruijn)

    # if there are still values in DeBruijn graph
    while any(deBruijn.values()) and len(path)>0:
        if len(potentialStartList) == 0:
            startNode = list(deBruijn.keys())[0]
        else:
            startNode = potentialStartList[0]
        path,deBruijn= findEulerianPath(startNode,deBruijn)#,outDegree,inDegree)
        potentialStartList,outDegree,inDegree = findStartingNode(deBruijn)
        
        if len(path)>0:
            pathList.append(path)
            starts.append(startNode)

    genomeParts = []
    for path in pathList: 
        if len(path)> 1:
            genom = findGenomeSequence(path)
        genomeParts.append(genom)

    if len(genomeParts)>1:
        genome=overlapsDict(genomeParts)
    else:
        genome = genomeParts[0]
    
    return genome,dictReadNames

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

def misMatches(genomeMatchIndex,sequence_read,genome):

    # mismatches to be returned in case of SNP/gap finding
    misMatchList = []

    # keeps track of the gap and mismatch tolerance
    insert = 0
    delete = 0
    mismatches = 0

    # where we are at in the sequences
    i = 0
    j = 0
    
    insertN = 0
    deleteN = 0
    gapsNeeded = len(sequence_read) - 50
    if gapsNeeded == 0:
        insertN = 1
        deleteN = 1
    elif gapsNeeded <0:
        deleteN = abs(gapsNeeded)
        #print('deleteN',deleteN)
    elif gapsNeeded> 0:
        insertN = gapsNeeded
    
    while i < len(sequence_read) and j < len(genome):
        #true index of the genome
        positionRefGenome = genomeMatchIndex + j
        baseRead = sequence_read[i]
        baseGenome = genome[j]

        if baseRead != baseGenome:
            #print('sequence read mistmach',sequence_read[i:])
            #print('sequence genome mistmach',genome[j:])

            potentialMismatch = False
            # if not equal first check if it is a mismatch
            if mismatches < SNPTolerance_:
                
                another = False
                if insert<insertN or delete>deleteN:
                    another = True
                # avoids the chances of length of genome and read not being the same
                length = min(len(genome[j+1:]),len(sequence_read[i+1:]))
                misMatch,conti = checkIndel(sequence_read[i+1:i+1+length],genome[j+1:j+1+length],mismatches,another,typeCheck = 'mismatch')

                # if we consider this one base change to be a mismatch and the rest is perfect match, our work is done
                if misMatch:
                    if not conti:
                        mistake = (positionRefGenome,'S',f' {baseGenome} {baseRead}')
                        misMatchList.append(mistake)
                        mismatches += 1
                        return True,misMatchList
                    else:
                        potentialMismatch = True
            # if the one base change doesn't work then check if this is an indel
            if insert < insertN or delete < deleteN: 
                insertion = False
                deletion = False
                
                # check if there is an insertion
                if insert < insertN:
                    another = False
                    if insert+1<insertN or delete<deleteN:
                        another = True
                    
                    if i+1 < len(sequence_read):
                        length =min(len(sequence_read[i+1:]),len(genome[j:]))
                        insertion,ICount = checkIndel(sequence_read[i+1:i+1+length],genome[j:j+length],mismatches,another)

                # check if there is a deletion
                if delete < deleteN:
                    
                    another = False
                    if insert<insertN or delete+1<deleteN:
                        another = True
                    if j+1 <len(sequence_read):
                        length = min(len(genome[j+1:]),len(sequence_read[i:]))
                        deletion,DCount = checkIndel(sequence_read[i:i+length],genome[j+1:j+1+length],mismatches,another)

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
                    delete+=1
                    j+=1
                    continue

                # count as insertion and adjust pointers
                elif insertion:
                    misMatchList.append((positionRefGenome,'I',f' {baseRead}'))
                    insert+=1
                    i+=1
                    continue

                #if both of them arn't true then consider it a mismatch and check again. This would be the case that there is more than 1 indel/mismatch in the sequence
                else: #if deletion and insertion both false consider that it might be a mismatch
                    # we would expect that if this isn't true then the mismatch tolerance will at some point go over and they won't match
                    if potentialMismatch:
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

def checkIndel(read,genome,mismatches,another,typeCheck='I/D'):
  # if indel then check that considering it an indel will be ok with the SNPtolerance
  
    if typeCheck == 'I/D':
        count = 0
        count = mismatches
        for i in range(len(read)):
            if read[i]!=genome[i]:
                count += 1
            if count < SNPTolerance_ or another:
                continue
            else:
                return False, 0
        return True, count
    else:
        # if mismatch then it needs to be a perfect match
        for i in range(len(read)):
            if read[i]!=genome[i]:
                if not another:
                    return False, False
                else:
                    return True, True
        
    return True, False

def NEWmatchReadMinimizers(genome,readMinimizer,ReadIndexFound,sequence_read, genomeMinimizersIndexes):
    listIndexes = []
    
    #for all of the places in which there is a double match in the genome calculate the spots we need to look for
    for matchIndex in genomeMinimizersIndexes:
        if len(sequence_read)!= 50:
            for i in range(-abs(len(sequence_read)-50)-1,abs(len(sequence_read)-50)+1):
                #initialize some variables
                mismatches = []
                indexGenome = matchIndex + i
                #full sequence match
                genomeSlice = genome[indexGenome-ReadIndexFound:len(sequence_read)+indexGenome-ReadIndexFound]
                
                trueAlignmentIndex = indexGenome-ReadIndexFound

                valid,mismatches= misMatches(trueAlignmentIndex,sequence_read,genomeSlice)

                #if valid occurence add to the dict of Errors
                if valid:
                    listIndexes.append((len(mismatches),trueAlignmentIndex))
        else:
            #initialize some variables
            mismatches = []
            
            #full sequence match
            genomeSlice = genome[matchIndex-ReadIndexFound:len(sequence_read)+matchIndex-ReadIndexFound]
        
            trueAlignmentIndex = matchIndex-ReadIndexFound

            valid,mismatches= misMatches(trueAlignmentIndex,sequence_read,genomeSlice)

            #if valid occurence add to the dict of Errors
            if valid:
                listIndexes.append((len(mismatches),trueAlignmentIndex))

    if len(listIndexes)>0:
        listIndexes.sort()
        return listIndexes[0][0],listIndexes[0][1]
    else:
        return 1000,trueAlignmentIndex               

def findAlignment(genome,spectrum):
    hashDict = hashGenome(genome)
    indexDict = {}
    numberReads = len(spectrum)
    count = 0
    for sequence_name, sequence_read in spectrum.items():

        count+=1
        tuplesList = []
        nestedDict = [sequence_read[:window_],sequence_read[window_:window_*2],sequence_read[window_*2:window_*3],sequence_read[window_*3:window_*4]]
        
        if len(nestedDict)> 0:
            #loop through all of the hashes for that read
            for index,minimizer in enumerate(nestedDict):
                if minimizer in hashDict:
                    mismatchesNum,matchIndex = NEWmatchReadMinimizers(genome,minimizer,index*window_,sequence_read, hashDict[minimizer])
                    tuplesList.append((mismatchesNum,matchIndex))
            
            tuplesList.sort()
            for i in tuplesList: 
                if i[1]<0:
                    tuplesList.remove(i) 
             
            if len(tuplesList)>0:
                indexDict[sequence_name] = i[1]
                
            #else:
                #print('DIDNT FIND ALIGNMENT!',sequence_name)
                

        if count%1000 == 0:
            print('Genomes Left for hash:',numberReads-count)
    return indexDict

def Project2(pairedreads):
    #g=read_fasta('project2_sample2_reference_genome.fasta')
    genome,sequenceReads = findGenome(pairedreads) # kmerSequence: name of reads
    
    
    # we have the correct genome, now I need to align it properly
    indexDict = findAlignment(genome,pairedreads)
    
    
    dictAnswer = {}
    for sequence_read,index in indexDict.items():
        if index not in dictAnswer:
            dictAnswer[index] = [sequence_read]
        else:
            print('collision of indicies')
            print(sequence_read)
            dictAnswer[index] += [sequence_read]


    answer = sorted(dictAnswer)
    answerList = []
    for i in answer:
        if i in dictAnswer:
            answerList.extend(dictAnswer[i])

    return answerList, genome


def find_variants(read_file):
    pairedReads = read_fasta(read_file)
    answerList,genome = Project2(pairedReads)

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