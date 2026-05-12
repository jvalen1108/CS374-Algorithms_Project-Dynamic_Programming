import time


def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence += line
    return sequence


def smith_waterman(seq1, seq2):
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Create scoring table. First row/col stay at 0 (key SW difference)
    score = [[0] * cols for _ in range(rows)]

    # Track the highest-scoring cell as we fill the matrix
    max_score = 0
    max_i = 0
    max_j = 0

    # Fill the table
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal = score[i - 1][j - 1] + match_score
            else:
                diagonal = score[i - 1][j - 1] + mismatch_penalty

            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty

            # Diffference from Needleman-Wunsch: floor at 0 instead of allowing negatives
            score[i][j] = max(0, diagonal, up, left)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_i = i
                max_j = j

    # Traceback from the highest-scoring cell, stopping when we hit a 0
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = max_i
    j = max_j
    end_i = max_i
    end_j = max_j

    while i > 0 and j > 0 and score[i][j] > 0:
        if seq1[i - 1] == seq2[j - 1]:
            diagonal_score = match_score
        else:
            diagonal_score = mismatch_penalty

        if score[i][j] == score[i - 1][j - 1] + diagonal_score:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    start_i = i + 1  # convert to 1-indexed for reporting
    start_j = j + 1

    return aligned_seq1, aligned_seq2, max_score, (start_i, end_i), (start_j, end_j)


def alignment_identity(aligned1, aligned2):
    if len(aligned1) == 0:
        return 0.0
    matches = 0
    for a, b in zip(aligned1, aligned2):
        if a == b and a != "-":
            matches += 1
    return 100.0 * matches / len(aligned1)


def format_alignment(a1, a2, offset_w, offset_o, width=60):
    lines = []
    for start in range(0, len(a1), width):
        chunk1 = a1[start:start + width]
        chunk2 = a2[start:start + width]
        track = ""
        for c1, c2 in zip(chunk1, chunk2):
            if c1 == c2 and c1 != "-":
                track += "|"
            else:
                track += " "
        lines.append(f"Wuhan   {offset_w + start:>5}  {chunk1}")
        lines.append(f"                {track}")
        lines.append(f"Omicron {offset_o + start:>5}  {chunk2}")
        lines.append("")
    return "\n".join(lines)


# Main Program
wuhan = read_fasta("DNA Sequences/sarscov2_wuhan_spike.fasta")
omicron = read_fasta("DNA Sequences/sarscov2_omicron_spike.fasta")

print(f"Wuhan spike length:   {len(wuhan)} aa")
print(f"Omicron spike length: {len(omicron)} aa")

start = time.perf_counter()
aligned_w, aligned_o, score, w_range, o_range = smith_waterman(wuhan, omicron)
elapsed = time.perf_counter() - start

identity = alignment_identity(aligned_w, aligned_o)

print()
print("=== Smith-Waterman: Wuhan vs Omicron spike protein ===")
print(f"Runtime:                 {elapsed:.4f} seconds")
print(f"Best local score:        {score}")
print(f"Alignment length:        {len(aligned_w)} aa")
print(f"Identity in region:      {identity:.2f}%")
print(f"Conserved region (Wuhan):   positions {w_range[0]}-{w_range[1]}")
print(f"Conserved region (Omicron): positions {o_range[0]}-{o_range[1]}")

print()
print("=== First 180 aa of the conserved local alignment ===")
print(format_alignment(aligned_w[:180], aligned_o[:180], w_range[0], o_range[0]))