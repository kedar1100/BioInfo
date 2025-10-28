"""
Assignment 3: Genomic Databases and Advanced Applications (CCA5)
Author: [Your Name]
Date: 28 Oct 2025

Part A - Q1: FASTA File Processing
---------------------------------
Implements reading, writing, and analyzing multiple FASTA sequences.
"""
"""
    Reads a FASTA file and returns a dictionary of sequences.
    
    Args:
        file_path (str): Path to the FASTA file.
    
    Returns:
        dict: {header: sequence}
    """

from typing import Dict, List

def read_fasta(file_path: str) -> Dict[str, str]:

    sequences = {}
    header = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                sequences[header] = ""
            else:
                sequences[header] += line
    return sequences


def write_fasta(sequences: Dict[str, str], output_path: str):
    """
    Writes sequences to a FASTA file with proper formatting.
    
    Args:
        sequences (dict): Dictionary of {header: sequence}.
        output_path (str): Output file path.
    """
    with open(output_path, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):  # Format: 80 chars per line
                f.write(seq[i:i+80] + "\n")

    """
    Extracts metadata (sequence length, GC content, etc.) for each sequence.
    
    Args:
        sequences (dict): {header: sequence}
    
    Returns:
        list: Metadata dictionaries for each sequence.
    """

def get_fasta_metadata(sequences: Dict[str, str]) -> List[Dict[str, str]]:

    metadata = []
    for header, seq in sequences.items():
        gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100 if len(seq) > 0 else 0
        metadata.append({
            "Header": header,
            "Length": len(seq),
            "GC_Content(%)": round(gc_content, 2)
        })
    return metadata


# Example usage and testing
if __name__ == "__main__":
    fasta_dict = read_fasta("example.fasta")
    print("Parsed Sequences:", fasta_dict)
    
    write_fasta(fasta_dict, "output.fasta")
    
    meta = get_fasta_metadata(fasta_dict)
    print("\nMetadata for sequences:")
    for m in meta:
        print(m)
