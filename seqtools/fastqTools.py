import random, gzip
import sys
import os
import numpy as np
from collections import defaultdict, deque


def bernoulli(probability):
    return np.random.binomial(1, float(probability))



def fastqToRaw(fqFile):
    ''' fastq structure is 4 lines per unit
        lines 1-4 == @name, sequence, +name, quality scores
        This extracts just the sequences (line 2s)
        it prints each so as to not store in memory
        Thus, this needs to be used from command line and redirected'''
    ctr = 0
    for line in open(fqFile, 'r'):
        ctr += 1
        if ctr%4 == 2:
            print(line[:-1])


def numReads(fqFile):
    with open(fqFile, 'r') as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return int((i+1)/4)


def extractRead(fqFile, readNum):
    '''the readNumTh is the read that occurs at line readNum*4'''
    ctr = 0
    if readNum == 0:
        totalNumRead = numReads(fqFile)
        readNum = random.randint(1,totalNumRead)
    f = open(fqFile, 'r')
    try:
        while f:
            ctr += 1
            if ctr == readNum:
                read = next(f) + next(f) + next(f) + next(f)
                sys.stdout.write(read)
                break
            else:
                for i in range(4): next(f)
    except StopIteration:
        pass


def shuffleFqFile(fqFile, outputFile=None, givenTotalNumReadsInFile=None):
    ''' 1 M reads will take up >= 3*8 MB = 24 MB memory
        10 M will take uup 240 MB
        100 M will take up 2.4 GB
        200 M will take up >= 4.8 GB
        Thus, this is probably only useful on files with < 100 M reads.
        A wrapper can be made that divides larger files up, shuffles smaller files, and randomly concats them back.'''
    if not givenTotalNumReadsInFile:
        totalNumReads = numReads(fqFile)
    else:
        totalNumReads = givenTotalNumReadsInFile
    readIndexes = deque(list(range(1, totalNumReads+1)))
        ## need this to be deque instead of list b/c lists are fast only for appending/popping items off end
        ## deques are fast at popping them off both sides: beginning and end (could also just choose to pop off end with list)
    random.shuffle(readIndexes)
    ## index file: open index file and read in forming dictionary with k:v == readNum:startInMem
    intIndex = internalIndex(fqFile)
    ## go through shuffled list (pop stuff off to shrink it), go to read location, extract read, and print read
    f = open(fqFile, 'r')
    if not outputFile:
        while readIndexes:
            f.seek(intIndex[readIndexes[0]])
            fastqEntry = next(f) + next(f) + next(f) + next(f)
            sys.stdout.write(fastqEntry)
            readIndexes.popleft()
    else: ## outputfile specified
        out = open(outputFile, 'w')
        while readIndexes:        
            f.seek(intIndex[readIndexes[0]])
            fastqEntry = next(f) + next(f) + next(f) + next(f)
            out.write(fastqEntry)
            readIndexes.popleft()      



def internalIndex(fqFile):
    '''This is the same as the "index" function only it is made for functions to use internally.
        Whereas index has O(constant) space req, this has O(n) space req.
        This is b/c "index" just goes through a file and writes to another file as it goes.
        And "internalIndex" stores a dictionary with n k,v pairs where n is number of reads.
        The k,v pairs are readIndex and byte it starts at (where a read's index is k if it is the kth read in file)
        "internalIndex" returns the dict
        Another difference is that index returns positions of actual sequence whereas this returns positions of read name
        '''
    f = open(fqFile, 'rb')
    f.seek(0)
    ## initial
    ctr = 0
    fqEntry = next(f) + next(f) + next(f) + next(f)
    size = sys.getsizeof(fqEntry) - 37
    sumSize = 1 
    intIndex = defaultdict(int)
    while fqEntry:
        ctr += 1
        intIndex[ctr] = sumSize ## the sum size when up to a 'line 1' is exact stating pos of a 'line2' b/c start is 0 index and bytes are 1-index
        try:
            fqEntry = next(f) + next(f) + next(f) + next(f)
        except StopIteration:
            break
        size = sys.getsizeof(fqEntry) - 37
        sumSize += size
    return intIndex
        

def open_fastq(fqFile):
    if fqFile.endswith('.gz'):
        return gzip.open(fqFile,'rb')
    else:
        return open(fqFile, 'r')

def downSampleReadsWithoutReplacement(fqFile, probability, outputFile=None):
    f = open_fastq(fqFile)
    if not outputFile:
        try:
            while f:
                read = next(f) + next(f) + next(f) + next(f)
                if bernoulli(probability):
                    sys.stdout.write(read)
        except StopIteration:
            pass
    else:
        out = open(outputFile, 'w')
        try:
            while f:
                read = next(f) + next(f) + next(f) + next(f)
                if bernoulli(probability):
                    out.write(read)
        except StopIteration:
            pass        


def downSamplePairsWithoutReplacement(fqFile1, fqFile2, probability, outputFile=None):
    if outputFile == None: ##not supporting stdout with pairs now, if I do in future, it would need to be interleaved.
        return
    f1 = open_fastq(fqFile1)
    f2 = open_fastq(fqFile2)
    outpre = (".").join(outputFile.split(".")[:-1])
    out1 = open(outpre+".1.fastq", 'w')
    out2 = open(outpre+".2.fastq", 'w')
    try:
        while f1 and f2:
            read1 = next(f1) + next(f1) + next(f1) + next(f1)
            read2 = next(f2) + next(f2) + next(f2) + next(f2)
            if bernoulli(probability):
                out1.write(read1)
                out2.write(read2)
    except StopIteration:
        pass        
    f1.close()
    f2.close()
    out1.close()
    out2.close()

       



def downSampleReadsWithReplacement(fqFile, numReadsToSample, outputFile=None):
    numReadsInFile = numReads(fqFile)
    getThisRead = defaultdict(int)
    for i in range(numReadsToSample):
        getThisRead[random.randint(1,numReadsInFile)] += 1
    readNum = 0
    f = open_fastq(fqFile)
    if not outputFile: ## stdout
        try:
            while f:
                readNum += 1
                read = next(f) + next(f) + next(f) + next(f)
                while getThisRead[readNum] > 0:
                    sys.stdout.write(read)
                    getThisRead[readNum] -= 1
        except StopIteration:
            pass
    else: ## write to file
        out = open(outputFile, 'w')
        try:
            while f:
                readNum += 1
                read = next(f) + next(f) + next(f) + next(f)
                while getThisRead[readNum] > 0:
                    out.write(read)
                    getThisRead[readNum] -= 1
        except StopIteration:
            pass
        
def sampleIntsWithReplacementDict(numReadsToSample):
    getThisRead = defaultdict(int)
    for i in range(numReadsToSample):
        getThisRead[random.randint(1,numReadsInFile)] += 1
    return getThisRead

def downSamplePairsWithReplacement(fqFile1, fqFile2, numReadsToSample, outputFile=None):
    if outputFile == None: ##not supporting stdout with pairs now, if I do in future, it would need to be interleaved.
        return
    numReadsInFile = numReads(fqFile)
    getThisRead = sampleWithReplacementDict(numReadsToSample)
    readNum = 0
    f1 = open_fastq(fqFile1)
    f2 = open_fastq(fqFile2)
    outpre = (".").join(outputFile.split(".")[:-1])
    out1 = open(outpre+"1.fastq", 'w')
    out2 = open(outpre+"3.fastq", 'w')
    try:
        while f1 and f2:
            readNum += 1
            read1 = next(f1) + next(f1) + next(f1) + next(f1)
            read2 = next(f2) + next(f2) + next(f2) + next(f2)
            while getThisRead[readNum] > 0:
                out1.write(read1)
                out2.write(read2)
                getThisRead[readNum] -= 1
    except StopIteration:
        pass
    f1.close()
    f2.close()
    out1.close()
    out2.close()
        
    
       
