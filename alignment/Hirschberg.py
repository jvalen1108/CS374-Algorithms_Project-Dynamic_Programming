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

    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    m = len(seq2)
    
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
        # For very small strings, use the standard NW 
        aligned1, aligned2, _ = needleman_wunsch_basic(seq1, seq2)
    else:
        # Divide seq1 in half
        mid1 = len(seq1) // 2

        # Find  score of the first half (Forward)
        score_left = nw_score(seq1[:mid1], seq2)
        
        # Find the score of the second half (Backward)
        # Reverse both strings to calculate scores from the end to the middle
        score_right = nw_score(seq1[mid1:][::-1], seq2[::-1])
        
        # Find the split point in seq2 (mid2) that max the total score
        # Total score at index j = score_left[j] + score_right[len(seq2) - j]
        max_v = -float('inf')
        mid2 = 0
        for j in range(len(seq2) + 1):
            if score_left[j] + score_right[len(seq2) - j] > max_v:
                max_v = score_left[j] + score_right[len(seq2) - j]
                mid2 = j

        # Recursively solve the two sub-problems
        top_left_1, top_left_2 = hirschberg(seq1[:mid1], seq2[:mid2])
        bottom_right_1, bottom_right_2 = hirschberg(seq1[mid1:], seq2[mid2:])
        
        aligned1 = top_left_1 + bottom_right_1
        aligned2 = top_left_2 + bottom_right_2

    return aligned1, aligned2

def needleman_wunsch_basic(seq1, seq2):

    match_score, mismatch_penalty, gap_penalty = 2, -1, -2
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score = [[0]*cols for _ in range(rows)]
    for i in range(1, rows): score[i][0] = score[i-1][0] + gap_penalty
    for j in range(1, cols): score[0][j] = score[0][j-1] + gap_penalty
    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            score[i][j] = max(diag, score[i-1][j] + gap_penalty, score[i][j-1] + gap_penalty)
    
    
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


human_sequence = read_fasta("DNA Sequences/human_hba.fasta")
import time

def score_from_alignment(a1, a2):
    """Compute alignment score from the aligned strings."""
    match_score, mismatch_penalty, gap_penalty = 2, -1, -2
    score = 0
    for c1, c2 in zip(a1, a2):
        if c1 == '-' or c2 == '-':
            score += gap_penalty
        elif c1 == c2:
            score += match_score
        else:
            score += mismatch_penalty
    return score

def identity(a1, a2):
    if len(a1) == 0: return 0.0
    matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != '-')
    return 100.0 * matches / len(a1)

def run(label, seq1, seq2):
    print(f"\n=== Hirschberg: {label} ===")
    print(f"Lengths: {len(seq1)} aa vs {len(seq2)} aa")
    start = time.perf_counter()
    a1, a2 = hirschberg(seq1, seq2)
    elapsed = time.perf_counter() - start
    sc = score_from_alignment(a1, a2)
    ident = identity(a1, a2)
    # Compare against full NW to verify same score
    _, _, nw_sc = needleman_wunsch_basic(seq1, seq2)
    match = "MATCH" if sc == nw_sc else f"MISMATCH (NW={nw_sc})"
    print(f"Score:           {sc}  (vs NW: {match})")
    print(f"Identity:        {ident:.2f}%")
    print(f"Alignment len:   {len(a1)} aa")
    print(f"Runtime:         {elapsed:.4f} seconds")

# Hemoglobin comparisons (verifies correctness against NW)
human = read_fasta("DNA Sequences/human_hba.fasta")
chimp = read_fasta("DNA Sequences/chimp_hba.fasta")
mouse = read_fasta("DNA Sequences/mouse_hba.fasta")

run("Human vs Chimp HBA",  human, chimp)
run("Human vs Mouse HBA",  human, mouse)
