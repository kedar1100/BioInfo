import DNAtoolkit

seq=DNAtoolkit.RandDNA(11)

print(DNAtoolkit.countNucFrequency(seq))
print(f"5' {DNAtoolkit.validateAndDisplayDna(seq)} 3'")
print(f"   {''.join(["|" for c in range (len(seq))])}")
print(f"3' {DNAtoolkit.reverse_complement(seq)} 5'")