import argparse
#Global Variables and Tunable Parameters
#Current Parameters seem to be the best functioning for 1b
#Finding window Parameters
window_ = 16
checkLength_ = 5

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
def hashGenome(genome,window = window_):
    hashDict = {}
    for i in range(len(genome)-window+1):
        kmer = genome[i:i+window]
        if kmer in hashDict:
            hashDict[kmer] += [i]
        else: 
            hashDict[kmer] = [i]
    return hashDict

# finds the mismatches stringwise
def misMatches(genomeMatchIndex,sequence_read,genome):
    misMatchList = []
    gaps = 0
    mismatches = 0
    i = 0
    j = 0
    
    while i < len(sequence_read) and j < len(genome):
        positionRefGenome = genomeMatchIndex + j
        baseRead = sequence_read[i]
        baseGenome = genome[j]
        if baseRead != baseGenome:
            # first check if it was a gap 
            if gaps < gapTolerance_:
                insertion = False
                deletion = False
                if i+1 < len(sequence_read):
                    if len(sequence_read[i+1:])>=checkLength_ and len(genome[j:])>=checkLength_:
                        insertion = checkIndel(sequence_read[i+1:],genome[j:])
                if j+1 <len(sequence_read):   
                    if len(genome[j+1:])>= checkLength_ and len(sequence_read[i:])>=checkLength_: 
                        deletion = checkIndel(sequence_read[i:],genome[j+1:])
                    
                if deletion and insertion: 
                    print('BOTH DELETION AND INSERTION!!')
                    return False,misMatchList
                elif deletion: 
                    misMatchList.append((positionRefGenome,'D',f' {baseGenome}')) 
                    gaps+=1
                    j+=1
                    continue 
                elif insertion: 
                    misMatchList.append((positionRefGenome,'I',f' {baseRead}'))
                    gaps+=1
                    i+=1
                    continue 
            
            # then check if it was a mismatch 
            if mismatches < SNPTolerance_:
                mistake = (positionRefGenome,'S',f' {baseGenome} {baseRead}')
                misMatchList.append(mistake)
                mismatches += 1
                i+=1
                j+=1
            else: 
                return False,misMatchList
        else: 
            i+=1
            j+=1

    #print(misMatchList)
    return True,misMatchList

#check if the base pairs match enough
def checkIndel(read,genome):
    for i in range(checkLength_):
        if read[i]!=genome[i]:
            return False
    return True

def findReadMatchesToGenomes(metaGenomics,listValid):
    answerDict = {}
    majorityDict = {}
    listOfNoMatches = []

    #go through all of the reads in the my answer dict
    for reads in metaGenomics.keys():
        # grab all of the genomes that are valid
        genomeList = [key for key, value in metaGenomics[reads].items() if value]
        
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

    #majority list in order
    genomeMajorityList = []
    for genome in answerDictSorted.keys():
        genomeMajorityList.append(genome)

    #if there are multiple genomes per dict, choose the one that appears most frequently
    for read in majorityDict.keys():
        for k in genomeMajorityList:
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
    #reads in the files that you want
    file_path_reads = reads
    reads_genomes = read_fasta(file_path_reads)

    import os
    genomeMinimizers = {}
    genomes = {}
    # Specify the folder path
    folder_path = reference_genome

    # List all files in the folder
    files = os.listdir(folder_path)
    listValid = []
    # Process each file based on its name
    for file_name in files:
        file_path = os.path.join(folder_path, file_name)

        # Check if the item is a file (not a directory)
        if os.path.isfile(file_path):
            # Process the file based on its name
            # You can add your processing logic here
            if 'project1c_genome_' in file_name and file_name.endswith('.fasta'): 
                reference_genome = read_fasta(file_path)
                print('Calculating minimizers for file:',file_name)
                for sequence_id, sequence_genome in reference_genome.items():
                    #reference genome that I am comparing to 
                    genomes[sequence_id.split()[0]] = sequence_genome 
                    #all of the minimizers in the sequence
                    #genomeMini = findMinimizers(sequence_genome)
                    genomeMini = hashGenome(sequence_genome) #TODO trying to figure out if kmers are better than minimizers
                    genomeMinimizers[sequence_id.split()[0]] = genomeMini
                    listValid.append(sequence_id.split()[0])




    # process the files
    count = 0
    metaGenomics = {}

    for sequence_id, sequence_read in reads_genomes.items():
        readDict = {}
        for i in genomeMinimizers.keys():
            readDict[i] = False
        metaGenomics[sequence_id] = readDict
        
        #find the minimizer for a given read 

        nestedDict = [sequence_read[:16],sequence_read[16:32],sequence_read[32:48]]
        
        count+=1
        
        if len(nestedDict)> 0:
            
            for genomeDict in genomeMinimizers.keys():
                
                
                #loop through all of the minimizers for that read
                for index,minimizer in enumerate(nestedDict):
                    #if the minimizer appears more than once deal with later
                    
                    #if the minimizer matches to a minimizer in the genome evaluate it
                    if minimizer in genomeMinimizers[genomeDict]:
                        #index for the minimizer in th read
                        
                        #match the minimizer to the location in the genome
                        #score,errors,alignmentIndex = NEWmatchReadMinimizers(genome,minimizer,minimizerIndex,sequence_read, genomeMinimizersDict[minimizer])
                        valid = NEWmatchReadMinimizers(genomes[genomeDict],minimizer,index*16,sequence_read, genomeMinimizers[genomeDict][minimizer])
                        #print(valid)
                        if valid: 
                            #print(genomeDict)

                            metaGenomics[sequence_id][genomeDict] = True
                            break
        
    print('Finished with processing, now filtering') 
    thresholdReadCounts_ = 10
    #check if they make the threshold
    listOfNoMatches = []
    answerDict,noMatch = findReadMatchesToGenomes(metaGenomics,listValid)

    validGenomes = [key for key in answerDict.keys() if answerDict[key][0]>=thresholdReadCounts_]

    # might not want to do this twice 
    answerDictFinal,listOfNoMatches = findReadMatchesToGenomes(metaGenomics,validGenomes)

    # assign all reads not matched at all to the genome that appears most frequently
    answerDictSorted = dict(sorted(answerDictFinal.items(), key=lambda item: item[1]))

    for genome in reversed(answerDictSorted.keys()):
        answerDictFinal[genome][0] += len(listOfNoMatches)
        answerDictFinal[genome][1] += listOfNoMatches
        break

    # name them 
    answerList = []
    for genome in answerDictFinal.keys():
        for read in answerDictFinal[genome][1]:
            answerList.append([f'{read} {genome}'])

    return answerList

def main():
    parser = argparse.ArgumentParser(description='Find variants in genomic data.')
    parser.add_argument('--reference-genome', required=True, help='Path to the reference genomes file folder, with naming convention "project1c_genome_"')
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