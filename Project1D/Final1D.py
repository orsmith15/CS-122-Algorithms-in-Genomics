import argparse
#Global Variables and Tunable Parameters
#Current Parameters seem to be the best functioning for 1b
# Window Parameters for reads 
window_ = 16

# Minimizer Parameters for intial 
kmer_ = 35 # 31 <--93%(40,10,4)
minimizerLength_ = 10# 7 <-- next to test

#SNP and Gap Tolerance per alignment
SNPTolerance_ = 3
gapTolerance_ = 2

# amount of reads to initially test against
sampleSize_ = 200000
threshold_ = 2000

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

# grabs the name of the genomes to be used later
def grabName(file_path):
    sequences = {}
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_sequence = line[1:]
                sequences[current_sequence] = ""
                break
    return sequences

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

def findMinimizers(DNApiece,kmer= kmer_,minimizerLength=minimizerLength_):
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

def findReadMatchesToGenomes(metaGenomics,listValid):
    answerDict = {}
    majorityDict = {}
    listOfNoMatches = []

    #go through all of the reads in the my answer dict
    for reads in metaGenomics.keys():
        # grab all of the genomes that are valid
        genomeList = metaGenomics[reads]
        
        # if there is only 1 valid genome add it to the dict
        if len(genomeList) == 1:
            if genomeList[0] in answerDict:
                answerDict[genomeList[0]][0] += 1
                answerDict[genomeList[0]][1] += [reads]
            else: 
                answerDict[genomeList[0]] = [1,[reads]]

        # if there is more than 1 process later
        elif len(genomeList)>1: 
            majorityDict[reads] = genomeList
        
        # if there is 0 process later 
        else:
            listOfNoMatches.append(reads)

    print('Number of reads not matched at all',len(listOfNoMatches))

    answerDictSorted = dict(sorted(answerDict.items(), key=lambda item: item[1]))

    # majority list in order
    genomeMajorityList = []
    for genome in answerDictSorted.keys():
        genomeMajorityList.append(genome)

    #if there are multiple genomes per dict, choose the one that appears most frequently
    for read in majorityDict.keys():
        for k in reversed(genomeMajorityList):
            if k in majorityDict[read]:
                answerDict[k][0] += 1
                answerDict[k][1] += [read]
                break 

    return answerDict,listOfNoMatches

# given a match to the split of the read check if it is a valid match spot
def NEWmatchReadMinimizers(genome,readMinimizer,ReadIndexFound,sequence_read, genomeMinimizersIndexes):
    dictErrors = []
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
            return True
            
    return False   

def find_variants(reference_genome, reads):

    #####   PART 1: FIND GENOMES PRESENT BY TAKING SAMPLE OF READS #####
    # reads in the file that we want
    file_path_reads = reads
    reads_genomes = read_fasta(file_path_reads)

    # dictionaries to store data from reads 
    # dict -- readID: {minimizers:[indexes]}
    readDictMinimizers = {}

    # dict -- minimizers: [readIDs]
    minimizerDict = {}

    # dict -- readID: readSequence
    seqReads = {}

    
    numberReads = len(reads_genomes.items())

    # number of reads to sample
    
    count = 0 

    # random sample to grab from reads
    import random
    randomList = random.sample(range(numberReads), sampleSize_)

    # metaGenomics -- readID: [genomes matched to]
    metaGenomics = {}
    for sequence_id, sequence_read in reads_genomes.items():
        
        # if index in randomList then calculate 
        if count not in randomList:
            continue 
        count += 1
        if count >= sampleSize_:
            break
        
        # initialize
        metaGenomics[sequence_id] = []

        #find the minimizers for read
        nestedDict = findMinimizers(sequence_read)
        
        # if there are minimizers, intialize all dicts
        if len(nestedDict)> 0:

            # dict with all of the sequences
            seqReads[sequence_id] = sequence_read

            # dict with all of the sequence_read: minimizer,index
            readDictMinimizers[sequence_id] = nestedDict

            #loop through all of the minimizers for that read
            for minimizer,index in nestedDict.items():
                if minimizer not in minimizerDict: # minimizerDict minimizer:list of reads with sequence
                    minimizerDict[minimizer] = [(sequence_id,index[0])]
                else: 
                    minimizerDict[minimizer] += [(sequence_id,index[0])]

        if count%10000 == 0:
            print('Genomes Left:',numberReads-count)

    print('Finished with taking a sample of the reads')


    ######  PART 1.2  MATCH ALL OF THE GENOMES TO THIS SAMPLE #####
    folder_path =  reference_genome
    import os
    
    numberOfGenomes = 5000
    
    listValid = []
    
    # go through all of the files
    for file_name in folder_path:
        file_path = os.path.join(folder_path, file_name)
    
        # Check if the item is a file and if the file is genome
        if os.path.isfile(file_path) and 'project1d_genome_' in file_name and file_name.endswith('.fasta'):
            # read file in
            referenceGenome = read_fasta(file_path)

            # go through sequence
            for sequence_id, sequence_genome in referenceGenome.items():
                # if we have already seen it skip it 
                if sequence_id.split()[0] in listValid:
                    print('Skipped duplicate read: ',sequence_id.split()[0])
                    continue
                else: 
                    # find the minimizers
                    genomeMini = findMinimizers(sequence_genome)
                    
                    # append to seen list
                    listValid.append(sequence_id.split()[0])

                    #for all of the minimizers 
                    for minimizer,index in genomeMini.items(): 
                    # check the look up table of minimizers generated from the reads 
                        if minimizer in minimizerDict:
                            # if it is in, check all of the reads that have this minimizer
                            for read in minimizerDict[minimizer]: # read = (readName,index)
                            # if we haven't already said that this read is in the genome check it
                                if sequence_id.split()[0] not in metaGenomics[read[0]]:
                                    # if valid add read
                                    valid = NEWmatchReadMinimizers(sequence_genome,minimizer,read[1],seqReads[read[0]],index)
                                    if valid:
                                        metaGenomics[read[0]] += [sequence_id.split()[0]]
                    numberOfGenomes -=1
    if numberOfGenomes%10 == 1:
        print('Genomes Left:',numberOfGenomes) 
    
    #############  filter out the ones that we want ############
    answerDict,noMatch = findReadMatchesToGenomes(metaGenomics,listValid)
    validGenomes = [key for key in answerDict.keys() if answerDict[key][0]>=threshold_]

    ######### PART 2 RERUN THE SCRIPT WITH JUST THE PRESENT GENOMES  ##################

    ##### PART 2.1 MAKE DICT OF GENOMES HASHES AND SEQUENCES  ###########
    listValid = []
    # dict -- genomeName : [hashes : [indexes]]
    genomeHash = {}

    # dict -- genomeName : genomeSequence
    genomesH = {}
    numberOfGenomes = 50

    # Process each file based on its name
    for valid in validGenomes:
        # this is based on the naming convention staying the same
        numGenome = valid.split('_')[-1]
        file_name = f'project1d_genome_{numGenome}.fasta'
        file_path = os.path.join(folder_path, file_name)
        
        # Process the file based on its name
        if os.path.isfile(file_path) and 'project1d_genome_' in file_name and file_name.endswith('.fasta'):

                referenceGenome = read_fasta(file_path)
                for sequence_id, sequence_genome in referenceGenome.items():
                    if sequence_id.split()[0] in listValid:
                        print('Skipped duplicate read: ',sequence_id.split()[0])
                        continue
                    else: 
                        genomesH[sequence_id.split()[0]] = sequence_genome
                        genomeMini = hashGenome(sequence_genome) 
                        genomeHash[sequence_id.split()[0]] = genomeMini
                        listValid.append(sequence_id.split()[0])
                        numberOfGenomes -=1
        
        if numberOfGenomes%10 == 0:
            print('Genomes Left:',numberOfGenomes)
    
    
    #################  PART 2.2  MATCH THE READS USING HASHES   #############

    numberReads = len(reads_genomes.items())


    metaGenomicsHash = {}
    count = 0
 
    for sequence_id, sequence_read in reads_genomes.items():
        count +=1
        metaGenomicsHash[sequence_id] = []

        # split based on size of window
        # TODO allow this to be true for any size read
        nestedDict = [sequence_read[:window_],sequence_read[window_:window_*2],sequence_read[window_*2:window_*3]]

        if len(nestedDict)> 0:
            # loop through all of the genomes
            for genomeDict in genomeHash.keys():
                #loop through all of the hashes for that read
                for index,minimizer in enumerate(nestedDict):
                    # if the minimizer matches to a minimizer in the genome evaluate it
                    if minimizer in genomeHash[genomeDict]:
                        valid = NEWmatchReadMinimizers(genomesH[genomeDict],minimizer,index*window_,sequence_read, genomeHash[genomeDict][minimizer])
                        # if valid mark it as such in the metaGenomics hash and continue to the next genome
                        if valid:
                            metaGenomicsHash[sequence_id] += [genomeDict]
                            break
                            
        if count%10000 == 0:
            print('Genomes Left for hash:',numberReads-count)
            
    print('Finished hashing')


    ####### PART FINAL TALLY UP THE COUNTS AND OUTPUT THE FINAL LIST  ############
    answerDictFinal,listOfNoMatches = findReadMatchesToGenomes(metaGenomicsHash,validGenomes)

    # assign all reads not matched at all to the genome that appears most frequently
    answerDictSorted = dict(sorted(answerDictFinal.items(), key=lambda item: item[1]))

    # assigning the reads
    for genome in reversed(answerDictSorted.keys()):
        answerDictFinal[genome][0] += len(listOfNoMatches)
        answerDictFinal[genome][1] += listOfNoMatches
        break

    # name them for submitting
    answerList = []
    for genome in answerDictFinal.keys():
        for read in answerDictFinal[genome][1]:
            answerList.append([f'{read} {genome}'])

    return answerList

def main():
    parser = argparse.ArgumentParser(description='Find variants in genomic data.')
    parser.add_argument('--reference-genome', required=True, help='Path to the reference genomes file folder, with naming convention "project1d_genome_{TAGGING-ID}.fasta"')
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
        print("List of reads and genomes:", mutations)

if __name__ == "__main__":
    main()