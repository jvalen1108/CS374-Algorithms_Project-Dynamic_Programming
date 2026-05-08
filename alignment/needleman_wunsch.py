def read_fasta(filename):

    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence += line
    return sequence

def needleman_wunsch(seq1, seq2):
    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -1

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Create scoring table
    score = []
    for i in range(rows):
        score.append([0] * cols)

    # Fill first column with gap penalties
    for i in range(1, rows):
        score[i][0] = score[i - 1][0] + gap_penalty

    # Fill first row with gap penalties
    for j in range(1, cols):
        score[0][j] = score[0][j - 1] + gap_penalty

    # Fill the rest of the table
    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal = score[i - 1][j - 1] + match_score
            else:
                diagonal = score[i - 1][j - 1] + mismatch_penalty

            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty

            score[i][j] = max(diagonal, up, left)

    # Traceback to build the alignment
    aligned_seq1 = ""
    aligned_seq2 = ""

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:

        if i > 0 and j > 0:
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = match_score
            else:
                diagonal_score = mismatch_penalty

        if score[i][j] == score[i - 1][j - 1] + diagonal_score:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        continue

    if i > 0 and score[i][j] == score[i - 1][j] + gap_penalty:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1

    else:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2, score[len(seq1)][len(seq2)]

# Main Program
human_sequence = read_fasta("human_hba.fasta")
chimp_sequence = read_fasta("chimp_hba.fasta")
print("Human Sequence:", human_sequence)
print("Chimp Sequence:", chimp_sequence)

aligned_human, aligned_chimp, final_score = needleman_wunsch(
    human_sequence,
    chimp_sequence
)

print("Human sequence:")
print(human_sequence)

print("\nChimp sequence:")
print(chimp_sequence)

print("\nAligned Human:")
print(aligned_human)

print("\nAligned Chimp:")
print(aligned_chimp)

print("\nFinal alignment score:")
print(final_score)