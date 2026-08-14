import random
import time
from collections import defaultdict

def generate_sequence_with_repeat(seq_len=100000, repeat_len=200, mutations=6, min_dist=10000):
    """Generates a random DNA sequence and implants a distant approximate repeat."""
    random.seed(42)
    alphabet = ['A', 'C', 'G', 'T']
    seq = [random.choice(alphabet) for _ in range(seq_len)]
    
    src_pos = 10000
    dst_pos = src_pos + min_dist
    
    # Copy segment and introduce mutations
    repeat_region = seq[src_pos : src_pos + repeat_len]
    mutated = list(repeat_region)
    
    mut_positions = random.sample(range(repeat_len), mutations)
    for pos in mut_positions:
        curr = mutated[pos]
        mutated[pos] = random.choice([b for b in alphabet if b != curr])
        
    seq[dst_pos : dst_pos + repeat_len] = mutated
    return "".join(seq), (src_pos, dst_pos)

# -------------------------------------------------------------------
# 1. Naive Indexing: Stores EVERY k-mer in the sequence
# -------------------------------------------------------------------
def build_naive_index(seq: str, k: int):
    index = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        index[kmer].append(i)
    return index

# -------------------------------------------------------------------
# 2. Minimizer Indexing: Stores only the minimum-hash k-mer per window
# -------------------------------------------------------------------
def build_minimizer_index(seq: str, k: int, w: int):
    index = defaultdict(list)
    seen = set()  # Avoid adding identical (pos, kmer) tuples across overlapping windows
    
    for i in range(len(seq) - (w + k - 1) + 1):
        min_hash = float('inf')
        min_pos = -1
        min_kmer = ""
        
        # Find lexicographical / hash minimum k-mer in sliding window of size w
        for j in range(w):
            pos = i + j
            kmer = seq[pos : pos + k]
            h = hash(kmer)
            if h < min_hash:
                min_hash = h
                min_pos = pos
                min_kmer = kmer
                
        if (min_pos, min_kmer) not in seen:
            seen.add((min_pos, min_kmer))
            index[min_kmer].append(min_pos)
            
    return index

# -------------------------------------------------------------------
# Seed Querying & Extension Phase
# -------------------------------------------------------------------
def find_distant_repeats(index, seq, k, min_dist=5000, max_mismatches=10, ext_len=200):
    candidates_evaluated = 0
    detected_pairs = set()
    
    for kmer, positions in index.items():
        if len(positions) < 2:
            continue
            
        # Compare positional pairs sharing the same seed
        for x in range(len(positions)):
            for y in range(x + 1, len(positions)):
                pos1, pos2 = positions[x], positions[y]
                
                # Filter out immediate / close repeats
                if pos2 - pos1 < min_dist:
                    continue
                    
                candidates_evaluated += 1
                
                # Extend seed and check Hamming distance over the repeat window
                if pos1 + ext_len <= len(seq) and pos2 + ext_len <= len(seq):
                    s1 = seq[pos1 : pos1 + ext_len]
                    s2 = seq[pos2 : pos2 + ext_len]
                    mismatches = sum(c1 != c2 for c1, c2 in zip(s1, s2))
                    
                    if mismatches <= max_mismatches:
                        detected_pairs.add((pos1, pos2))
                        
    return detected_pairs, candidates_evaluated

# -------------------------------------------------------------------
# Main Execution & Benchmark
# -------------------------------------------------------------------
if __name__ == "__main__":
    K = 15          # Seed length
    W = 10          # Window size for minimizers
    MUTATIONS = 6   # Substituted bases in 200bp repeat
    SEQ_LEN = 100000

    seq, (src, dst) = generate_sequence_with_repeat(seq_len=SEQ_LEN, repeat_len=200, mutations=MUTATIONS)
    print(f"Sequence length: {SEQ_LEN} bp | Planted repeat at positions ({src}, {dst})\n")

    # --- Benchmark Naive Seed-and-Extend ---
    t0 = time.time()
    naive_idx = build_naive_index(seq, K)
    t1 = time.time()
    naive_repeats, naive_candidates = find_distant_repeats(naive_idx, seq, K, min_dist=5000, max_mismatches=10)
    t2 = time.time()

    # --- Benchmark Minimizer Seed-and-Extend ---
    t3 = time.time()
    mini_idx = build_minimizer_index(seq, K, W)
    t4 = time.time()
    mini_repeats, mini_candidates = find_distant_repeats(mini_idx, seq, K, min_dist=5000, max_mismatches=10)
    t5 = time.time()

    # Aggregate entry counts
    naive_total_entries = sum(len(v) for v in naive_idx.values())
    mini_total_entries = sum(len(v) for v in mini_idx.values())

    # Results Comparison Table
    print("=" * 75)
    print(f"{'Metric':<30} | {'Naive (All k-mers)':<18} | {'Minimizer (w=10)':<18}")
    print("=" * 75)
    print(f"{'Indexed Position Entries':<30} | {naive_total_entries:<18} | {mini_total_entries:<18}")
    print(f"{'Index Build Time (s)':<30} | {t1 - t0:<18.4f} | {t4 - t3:<18.4f}")
    print(f"{'Candidate Pairs Evaluated':<30} | {naive_candidates:<18} | {mini_candidates:<18}")
    print(f"{'Repeat Search Time (s)':<30} | {t2 - t1:<18.4f} | {t5 - t4:<18.4f}")
    print(f"{'Planted Repeat Detected?':<30} | {str(any(abs(p[0]-src) < 50 and abs(p[1]-dst) < 50 for p in naive_repeats)):<18} | {str(any(abs(p[0]-src) < 50 and abs(p[1]-dst) < 50 for p in mini_repeats)):<18}")
    print("=" * 75)