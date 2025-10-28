import DNAtoolkit
from utilitiess import colored


seq=DNAtoolkit.RandDNA(200)

seq=DNAtoolkit.validateAndDisplayDna(seq)
colored(seq)
print(seq)
print(f"gc content {DNAtoolkit.gcContent(seq)} %")
print(f"gc content in subseq : {DNAtoolkit.gcContentSubSeq(seq,k=5)}")

print(f"gc skew in seq :{DNAtoolkit.gcSkew(seq)}")
print(f"comparing gc content {DNAtoolkit.compareGcContent(seq)}")
print(f"dinucleotide frequency {DNAtoolkit.dinucleotideFrequency(seq)}")
print(f"trinucleotide frequency {DNAtoolkit.trinucleotideFrequency(seq)}")
print(f"{DNAtoolkit.checkCpGIsland(seq)}")