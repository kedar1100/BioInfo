import collections
import random

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