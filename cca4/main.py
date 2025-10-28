import DNAtoolkit

seq=DNAtoolkit.RandDNA(5)
print(DNAtoolkit.validateAndDisplayDna(seq))
print(f"gc content {DNAtoolkit.gcContent(seq)} %")
print(f"gc content in subseq k=5: {DNAtoolkit.gcContentSubSeq(seq,k=5)}")

print(f"gc skew in seq :{DNAtoolkit.gcSkew(seq)}")
print(f"{DNAtoolkit.compareGcContent(seq)}")
print(f"{DNAtoolkit.dinucleotideFrequency(seq)}")
print(f"{DNAtoolkit.trinucleotideFrequency(seq)}")