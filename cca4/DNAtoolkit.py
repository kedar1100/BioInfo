import collections 
from collections import Counter
import random
from structure import DNACodon

Nucleotide=["A","C","G","T"]
NucleotideRNA=["T","G","C","A"]
DNAReverseComplement={"A":'T','T':'A','G':'C','C':'G'}


def validateSeq(dna_seq):
    tmpseq=dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotide:
            return False
    return tmpseq

# Q2
def countNucFrequency(seq):
    tmpFreqDict={"A": 0,"C":0,"G":0,"T":0}
    for nuc in seq:
        tmpFreqDict[nuc]+=1
    return tmpFreqDict
    # return dict(collections.Counter(seq))  

def countNucFrequencyPercentage(seq):
    a=countNucFrequency(seq)
    totalCount=sum(a.values())
    percentageDict={key:(value/totalCount*100)for key, value in a.items()}
    return percentageDict

# Q3
def validateAndDisplayDna(dna_seq):
    for nuc in dna_seq:
        if nuc not in Nucleotide:
            return False
    return dna_seq

# def splitIntoCodon(seq):
#     seq[]
#     seq=''.join("|"for _ in range(3))
#     return seq


# Random dna generator to practice upon operations 
def RandDNA(size):
    seq={}
    seq="".join(random.choice(Nucleotide) for _ in range(size))
    return seq
# DNA -> RNA Transcription thymine with uracil

def TranscriptionDNARNA(seq):
    return seq.replace("T","U")


def reverse_complement(seq):
    # return ''.join([DNAReverseComplement[nuc] for nuc in seq])[::-1]

    # optimizied and more pythonic code for the above function 
    mapping=str.maketrans('ATGC','TACG')
    return seq.translate(mapping)[::-1]



# CCA 4 functions
def gcContent(seq):
    return (round(((seq.count("C")+seq.count("G"))/len(seq))*100))


def gcContentSubSeq(seq,k):
    res =[]
    for i in range(0,len(seq)-k+1,k):
        subseq=seq[i:i+k]
        res.append(gcContent(subseq))
    return res

def GCSkew(seq):
    return ((seq.count("G")-seq.count("C"))/(seq.count("G")+seq.count("C")))

def ATSkew(seq):
    return((seq.count("A")-(seq.count("T")))/(seq.count("A"))+(seq.count("T")))

def compareGcContent(seq):
    g=seq.count("G")
    c=seq.count("C")
    return (((g+c)/len(seq))*100)


# this function can be re written by you my boy , you are fucking smart and can do it 😉

def dinucleotideFrequency(seq):
    # The total number of possible dinucleotides in a sequence of length L is L−1 (because each pair overlaps with the next, and you lose one potential pair at the end).
    seqLen=len(seq)-1
    count=[]
    i=0
    freqList={}
    for i in range(seqLen):
        count.append(seq[i]+seq[i+1])
        i+=1
    return count
    


def dinucleotide_frequency(sequence):
    seq = ''.join([base for base in sequence.upper() if base in "ACGT"])
    dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
    counts = Counter(dinucs)
    total = sum(counts.values())
    freqs = {k: v / total for k, v in counts.items()}
    
    return counts, freqs

def trinucleotideFrequency(sequence):
    seq = ''.join([base for base in sequence.upper() if base in "ACGT"])
    dinucs = [seq[i:i+2] for i in range(len(seq) - 2)]
    counts = Counter(dinucs)
    total = sum(counts.values())
    freqs = {k: v / total for k, v in counts.items()}
    return counts,freqs

def checkCpGIsland(seq):
    c=seq.count("C")
    g=seq.count("G")
    CpG=seq.count("CG")
    L=len(seq)
    expectedCpG=(c*g)/L
    CpGRatio=CpG/expectedCpG
    if CpGRatio<0.6:
        return f"Not a CpG Island : low cpg rati {CpGRatio}"
    
    return f"Is a CpG Island CpG ratio: {CpGRatio}"

def translateSeq(seq,init_pos=0):
    return[DNACodon[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]

def codonUsage(seq,aminoacid):
    a=aminoacid.upper()
    tmpList=[]
    for i in range (0,len(seq)-2,3):
        if DNACodon[seq[i:i+3]]==a:
            tmpList.append(seq[i:i+3])
    freqDict=dict(Counter(tmpList))
    totalWeight=sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq]=round(freqDict[seq]/totalWeight,2)
    return freqDict

def genReadingFrames(seq):
    frames=[]
    frames.append(translateSeq(seq,0))
    frames.append(translateSeq(seq,1))
    frames.append(translateSeq(seq,2))
    frames.append(translateSeq(reverse_complement(seq),0))
    frames.append(translateSeq(reverse_complement(seq),1))
    frames.append(translateSeq(reverse_complement(seq),2))
    return frames

def proteinsFromRF(aa_seq):
    currentProt=[]
    proteins=[]
    for aa in aa_seq:
        if aa=="_":
            if currentProt:
                for p in currentProt:
                    proteins.append(p)
                currentProt=[]
        else:
            if aa=="M":
                currentProt.append("")
            for i in range(len(currentProt)):
                currentProt[i]+=aa
    return proteins

# def allProteinFromORFS(seq,startReadPos=0,endReadPos=0):
#     if endReadPos>startReadPos:
#         rfs=genReadingFrames(seq[startReadPos:endReadPos])
#     else:
#         rfs=genReadingFrames(seq)
    
#     for rf in rfs:
        