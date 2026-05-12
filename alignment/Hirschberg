def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence += line
    return sequence

def nw_score(seq1, seq2):
    """
    Computes only the last row of the Needleman-Wunsch matrix.
    Space complexity: O(len(seq2))
    """
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    m = len(seq2)
    # prev_row stores the scores for the previous character in seq1
    prev_row = [j * gap_penalty for j in range(m + 1)]
    curr_row = [0] * (m + 1)

    for i in range(1, len(seq1) + 1):
        curr_row[0] = i * gap_penalty
        for j in range(1, m + 1):
            score_sub = prev_row[j - 1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            score_del = prev_row[j] + gap_penalty
            score_ins = curr_row[j - 1] + gap_penalty
            curr_row[j] = max(score_sub, score_del, score_ins)
        # Move current to previous for the next iteration
        prev_row[:] = curr_row[:]
        
    return prev_row

def hirschberg(seq1, seq2):
    """
    Recursive Divide and Conquer alignment.
    """
    aligned1 = ""
    aligned2 = ""

    # Base Cases
    if len(seq1) == 0:
        for char in seq2:
            aligned1 += '-'
            aligned2 += char
    elif len(seq2) == 0:
        for char in seq1:
            aligned1 += char
            aligned2 += '-'
    elif len(seq1) == 1 or len(seq2) == 1:
        # For very small strings, use the standard NW to finish it off
        aligned1, aligned2, _ = needleman_wunsch_basic(seq1, seq2)
    else:
        # 1. Divide seq1 in half
        mid1 = len(seq1) // 2

        # 2. Find the score of the first half (Forward)
        score_left = nw_score(seq1[:mid1], seq2)
        
        # 3. Find the score of the second half (Backward)
        # We reverse both strings to calculate scores from the end to the middle
        score_right = nw_score(seq1[mid1:][::-1], seq2[::-1])
        
        # 4. Find the split point in seq2 (mid2) that maximizes the total score
        # Total score at index j = score_left[j] + score_right[len(seq2) - j]
        max_v = -float('inf')
        mid2 = 0
        for j in range(len(seq2) + 1):
            if score_left[j] + score_right[len(seq2) - j] > max_v:
                max_v = score_left[j] + score_right[len(seq2) - j]
                mid2 = j

        # 5. Conquer: Recursively solve the two sub-problems
        top_left_1, top_left_2 = hirschberg(seq1[:mid1], seq2[:mid2])
        bottom_right_1, bottom_right_2 = hirschberg(seq1[mid1:], seq2[mid2:])
        
        aligned1 = top_left_1 + bottom_right_1
        aligned2 = top_left_2 + bottom_right_2

    return aligned1, aligned2

def needleman_wunsch_basic(seq1, seq2):
    """
    Small-scale NW used as a base case for the recursion.
    """
    match_score, mismatch_penalty, gap_penalty = 2, -1, -2
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score = [[0]*cols for _ in range(rows)]
    for i in range(1, rows): score[i][0] = score[i-1][0] + gap_penalty
    for j in range(1, cols): score[0][j] = score[0][j-1] + gap_penalty
    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            score[i][j] = max(diag, score[i-1][j] + gap_penalty, score[i][j-1] + gap_penalty)
    
    # Simple traceback for small sequences
    res1, res2 = "", ""
    i, j = len(seq1), len(seq2)
    final_s = score[i][j]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            res1, res2 = seq1[i-1] + res1, seq2[j-1] + res2
            i, j = i-1, j-1
        elif i > 0 and score[i][j] == score[i-1][j] + gap_penalty:
            res1, res2 = seq1[i-1] + res1, "-" + res2
            i -= 1
        else:
            res1, res2 = "-" + res1, seq2[j-1] + res2
            j -= 1
    return res1, res2, final_s

# --- Execution ---
human_sequence = read_fasta("DNA Sequences/human_hba.fasta")
chimp_sequence = read_fasta("DNA Sequences/chimp_hba.fasta")

aligned_human, aligned_chimp = hirschberg(human_sequence, chimp_sequence)

print(f"Aligned Human:\n{aligned_human}")
print(f"\nAligned Chimp:\n{aligned_chimp}")
