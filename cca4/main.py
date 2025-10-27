import DNAtoolkit

seq=DNAtoolkit.RandDNA(50)
print(DNAtoolkit.validateAndDisplayDna(seq))
print(f"{DNAtoolkit.gcContent(seq)} %")
print(f"gc content in subseq k=5: {DNAtoolkit.gcContentSubSeq(seq,k=5)}")