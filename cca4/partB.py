import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any
"""
    Calculates the Hamming distance between two sequences.
    
    The Hamming distance is defined as the number of positions at which the 
    corresponding symbols are different. If sequences are of unequal length, 
    the distance is calculated only over the length of the shorter sequence.

    Args:
        seq1 (str): The first biological sequence.
        seq2 (str): The second biological sequence.

    Returns:
        Tuple[int, str]: The calculated distance and a status message.
"""

# --- QUESTION 3: HAMMING DISTANCE IMPLEMENTATION ---

def hamming_distance(seq1: str, seq2: str) -> Tuple[int, str]:
    if not seq1 or not seq2:
        return 0, "Warning: One or both sequences are empty."

    n = min(len(seq1), len(seq2))
    distance = sum(c1 != c2 for c1, c2 in zip(seq1[:n], seq2[:n]))
    
    status = ""
    if len(seq1) != len(seq2):
        status = f"Note: Sequences have different lengths ({len(seq1)} vs {len(seq2)}). Distance calculated over the first {n} characters."
    else:
        status = f"Sequences are of equal length ({n})."
        
    return distance, status
    """
    Calculates the matrix of pairwise Hamming distances for multiple sequences.

    Args:
        sequences (Dict[str, str]): A dictionary mapping sequence names to sequence strings.

    Returns:
        pd.DataFrame: A symmetrical DataFrame representing the distance matrix.

    """
def calculate_distance_matrix(sequences: Dict[str, str]) -> pd.DataFrame:
    names = list(sequences.keys())
    seq_list = list(sequences.values())
    n = len(names)
    matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            dist, _ = hamming_distance(seq_list[i], seq_list[j])
            matrix[i, j] = dist
            matrix[j, i] = dist 
    distance_df = pd.DataFrame(matrix, index=names, columns=names)
    return distance_df

# --- QUESTION 4: MOTIF FINDING ALGORITHMS ---

def exact_pattern_matching(text: str, pattern: str) -> List[int]:
    """
    Finds all start indices of exact occurrences of a pattern in a text.
    Implements a simple sliding window/substring search algorithm.

    Args:
        text (str): The genomic sequence (haystack).
        pattern (str): The motif to search for (needle).

    Returns:
        List[int]: A list of starting indices where the pattern occurs exactly.
    """
    if not pattern:
        return []
        
    pattern_len = len(pattern)
    text_len = len(text)
    occurrences = []

    # Slide the window across the text
    for i in range(text_len - pattern_len + 1):
        if text[i:i + pattern_len] == pattern:
            occurrences.append(i)
            
    return occurrences

def fuzzy_pattern_matching(text: str, pattern: str, max_mismatches: int) -> List[Tuple[int, int]]:
    if not pattern:
        return []

    pattern_len = len(pattern)
    text_len = len(text)
    fuzzy_matches = []

    # Slide the window across the text
    for i in range(text_len - pattern_len + 1):
        substring = text[i:i + pattern_len]
        
        # Calculate mismatch count (Hamming distance)
        mismatch_count = sum(c1 != c2 for c1, c2 in zip(substring, pattern))
        
        if mismatch_count <= max_mismatches:
            fuzzy_matches.append((i, mismatch_count))

    return fuzzy_matches

def analyze_motif_conservation(sequences: Dict[str, str], pattern: str, max_mismatches: int) -> Dict[str, Any]:
    """
    Analyzes motif conservation by finding the best fuzzy match in each sequence.

    Args:
        sequences (Dict[str, str]): Sequences to search within.
        pattern (str): The motif pattern.
        max_mismatches (int): The maximum allowed mismatches.

    Returns:
        Dict[str, Any]: Analysis results including best match details and overall conservation metrics.
    """
    conservation_analysis = {}
    best_mismatch_counts = []

    for name, seq in sequences.items():
        matches = fuzzy_pattern_matching(seq, pattern, max_mismatches)
        
        if matches:
            # Find the match with the minimum number of mismatches
            best_match = min(matches, key=lambda x: x[1])
            index, mismatches = best_match
            
            best_mismatch_counts.append(mismatches)
            
            conservation_analysis[name] = {
                "Best Index": index,
                "Mismatches": mismatches,
                "Sequence": seq[index:index + len(pattern)],
                "Status": "Found best match"
            }
        else:
            conservation_analysis[name] = {
                "Status": f"No match found with <={max_mismatches} mismatches."
            }

    # Calculate overall conservation metrics
    total_mismatches = sum(best_mismatch_counts)
    num_found = len(best_mismatch_counts)
    
    if num_found > 0:
        avg_mismatches = total_mismatches / num_found
    else:
        avg_mismatches = np.nan
        
    return {
        "Conservation_Metrics": {
            "Total Sequences Analyzed": len(sequences),
            "Sequences with Match": num_found,
            "Average Mismatches (Best Match)": avg_mismatches,
            "Total Mismatches (Best Match)": total_mismatches
        },
        "Sequence_Details": conservation_analysis
    }

# --- EXAMPLE USAGE ---

if __name__ == '__main__':
    
    # Example Dataset (Transcription Factor Binding Sites)
    DNA_SEQUENCES = {
        "Seq_A": "GATTACAATGTCA",
        "Seq_B": "GATTACTCAGTGA",
        "Seq_C": "GATTACAGG",
        "Seq_D": "ATTCAGTTCA"
    }
    
    print("=====================================================")
    print("                QUESTION 3: HAMMING DISTANCE         ")
    print("=====================================================")

    # 1. Hamming Distance Implementation (Handling different lengths)
    seq1 = DNA_SEQUENCES["Seq_A"]
    seq2 = DNA_SEQUENCES["Seq_C"]
    distance, status = hamming_distance(seq1, seq2)
    print(f"\nComparing '{seq1}' and '{seq2}'")
    print(f"Hamming Distance: {distance}")
    print(f"Details: {status}")

    seq3 = DNA_SEQUENCES["Seq_A"]
    seq4 = DNA_SEQUENCES["Seq_B"]
    distance, status = hamming_distance(seq3, seq4)
    print(f"\nComparing '{seq3}' and '{seq4}'")
    print(f"Hamming Distance: {distance}")
    print(f"Details: {status}")
    
    # 2. Distance Matrix Calculation (Visualization Data Structure)
    print("\n--- Pairwise Hamming Distance Matrix (Data Structure) ---")
    matrix_df = calculate_distance_matrix(DNA_SEQUENCES)
    print(matrix_df.to_markdown(numalign="left", stralign="left"))
    print("\n*Note: This DataFrame provides the appropriate data structure for visualization (e.g., using a Seaborn heatmap or clustering dendrogram in a real environment).")


    print("\n=====================================================")
    print("             QUESTION 4: MOTIF FINDING ALGORITHMS    ")
    print("=====================================================")
    
    # Test sequence and motif
    genomic_text = "GATTACAATGTCAATGTCAATAGAT"
    motif_pattern = "ATGTCA"
    
    # 1. Exact Pattern Matching
    print(f"\n--- 4a: Exact Pattern Matching for Motif '{motif_pattern}' ---")
    exact_occurrences = exact_pattern_matching(genomic_text, motif_pattern)
    print(f"Text: {genomic_text}")
    print(f"Exact occurrences found at indices: {exact_occurrences}")
    
    # 2. Fuzzy Pattern Matching (with allowed mismatches)
    max_mismatches = 1
    print(f"\n--- 4b: Fuzzy Pattern Matching (Max Mismatches k={max_mismatches}) ---")
    fuzzy_occurrences = fuzzy_pattern_matching(genomic_text, "AATGTCA", max_mismatches)
    print(f"Text: {genomic_text}")
    print(f"Fuzzy matches (Index, Mismatch) found: {fuzzy_occurrences}")

    # 3. Motif Conservation Analysis across multiple sequences
    motif_to_find = "ATTACA"
    max_mismatches_for_conservation = 1
    print(f"\n--- 4c: Motif Conservation Analysis for Motif '{motif_to_find}' (Max Mismatches={max_mismatches_for_conservation}) ---")
    conservation_results = analyze_motif_conservation(DNA_SEQUENCES, motif_to_find, max_mismatches_for_conservation)
    
    # Print Conservation Metrics
    print("\n[Conservation Metrics]")
    metrics = conservation_results["Conservation_Metrics"]
    for k, v in metrics.items():
        print(f"- {k}: {v}")
        
    # Print Sequence Details
    print("\n[Sequence Details (Best Match)]")
    details_df = pd.DataFrame(conservation_results["Sequence_Details"]).T
    print(details_df.to_markdown(numalign="left", stralign="left"))
