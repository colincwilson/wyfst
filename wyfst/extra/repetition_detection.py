# Find (quasi-)repeated substrings within string using sliding window comparison.
import re, sys
import polars as pl
import edlib

sys.path.append("~/Code/Python/phonopy/")
from phonopy import str_util

# # # # # # # # # #


def find_internal_repeats(text,
                          chunk_size=4,
                          max_distance=0,
                          min_gap=None,
                          sep=" "):
    """Find (quasi-)repeated substrings within string."""
    if min_gap is None:
        min_gap = chunk_size

    char_vec = text.split(sep)
    n_chars = len(char_vec)
    n_chunks = n_chars - chunk_size + 1

    if n_chunks < 2:
        return None

    chunks = [sep.join(char_vec[i:(i + chunk_size)]) for i in range(n_chunks)]

    results = []

    for i in range(n_chunks - 1):
        start_j = i + min_gap
        if start_j > n_chunks:
            continue

        for j in range(start_j, n_chunks):
            d = edlib.align(chunks[i], chunks[j])['editDistance']

            if d <= max_distance:
                results.append({
                    'seq1': chunks[i],
                    'pos1': i + 1,  # R uses 1-based indexing
                    'seq2': chunks[j],
                    'pos2': j + 1,  # R uses 1-based indexing
                    'gap': j - (i + chunk_size),
                    'lv_dist': d
                })

    if not results:
        return None

    df = pl.DataFrame(results)
    return df.unique()


def find_all_repeats(text,
                     min_chunk_size=2,
                     max_chunk_size=10,
                     max_distance=0,
                     min_gap=None,
                     sep=" "):
    """Find all repeated chunks: returns all repeats for all chunk sizes."""
    # Find repeats for all chunk sizes
    all_repeats = []
    for size in range(min_chunk_size, max_chunk_size + 1):
        repeats = find_internal_repeats(text,
                                        chunk_size=size,
                                        max_distance=max_distance,
                                        min_gap=min_gap,
                                        sep=sep)

        if repeats is not None and len(repeats) > 0:
            repeats = repeats.with_columns(pl.lit(size).alias('chunk_size'))
            all_repeats.append(repeats)

    if len(all_repeats) == 0:
        return None

    combined = pl.concat(all_repeats)
    combined = combined.sort(['chunk_size', 'pos1', 'pos2'],
                             descending=[True, False, False])
    return combined


def find_maximal_repeats(text,
                         min_chunk_size=2,
                         max_chunk_size=10,
                         max_distance=0,
                         min_gap=None,
                         sep=" "):
    """Find maximal (longest) repeated chunks: filters out repeated
    chunks that are completely contained within longer chunks"""
    # Find repeats for all chunk sizes
    combined = find_all_repeats(text,
                                min_chunk_size=min_chunk_size,
                                max_chunk_size=max_chunk_size,
                                max_distance=max_distance,
                                min_gap=min_gap,
                                sep=sep)
    if combined is None or len(combined) == 0:
        return None

    # Filter: keep only chunks that are not subsumed by longer chunks
    # A chunk (i1, j1) with size s1 is subsumed by (i2, j2) with size s2 if:
    # s2 > s1 AND i2 <= i1 AND j2 <= j1 AND (i2 + s2) >= (i1 + s1) AND (j2 + s2) >= (j1 + s1)
    maximal = [combined.row(0, named=True)]  # Start with longest

    for i in range(1, len(combined)):
        current = combined.row(i, named=True)
        is_subsumed = False

        for larger in maximal:
            # Check if current is contained in larger
            if (larger['chunk_size'] > current['chunk_size']
                    and larger['pos1'] <= current['pos1']
                    and larger['pos2'] <= current['pos2']
                    and (larger['pos1'] + larger['chunk_size'])
                    >= (current['pos1'] + current['chunk_size'])
                    and (larger['pos2'] + larger['chunk_size'])
                    >= (current['pos2'] + current['chunk_size'])):
                is_subsumed = True
                break

        if not is_subsumed:
            maximal.append(current)

    result = pl.DataFrame(maximal)
    return result.sort(['chunk_size', 'pos1'], descending=[True, False])


def find_best_repeats(
    text,
    min_chunk_size=2,
    max_chunk_size=10,
    max_distance=0,
    min_gap=None,
    sep=" ",
    length_weight=2,  # Weight for chunk length
    distance_penalty=1,  # Penalty for Levenshtein distance
    gap_penalty=0.01,  # Penalty for gap size
    top_n=None  # Return only top N chunks (None = all)
):
    """Score-based selection: assigns a score to each repeat
    based on weighted length, distance, and gap."""
    all_repeats = []

    for size in range(min_chunk_size, max_chunk_size + 1):
        repeats = find_internal_repeats(text,
                                        chunk_size=size,
                                        max_distance=max_distance,
                                        min_gap=min_gap,
                                        sep=sep)

        if repeats is not None and len(repeats) > 0:
            repeats = repeats.with_columns(pl.lit(size).alias('chunk_size'))
            all_repeats.append(repeats)

    if len(all_repeats) == 0:
        return None

    combined = pl.concat(all_repeats)

    # Calculate score: reward length, penalize distance and gap
    combined = combined.with_columns(
        ((pl.col('chunk_size') * length_weight) -
         (pl.col('lv_dist') * distance_penalty) -
         (pl.col('gap') * gap_penalty)).alias('score'))

    combined = combined.sort('score', descending=True)

    if top_n is not None:
        combined = combined.head(min(top_n, len(combined)))

    return combined


def find_nonoverlapping_repeats(text,
                                min_chunk_size=2,
                                max_chunk_size=10,
                                max_distance=0,
                                min_gap=None,
                                sep=" ",
                                length_weight=2,
                                distance_penalty=1,
                                gap_penalty=0.01):
    """Greedy non-overlapping selection: selects the best
    non-overlapping set of repeats."""
    # First get all scored repeats
    all_scored = find_best_repeats(text,
                                   max_chunk_size,
                                   min_chunk_size,
                                   max_distance,
                                   min_gap,
                                   sep,
                                   length_weight,
                                   distance_penalty,
                                   gap_penalty,
                                   top_n=None)

    if all_scored is None or len(all_scored) == 0:
        return None

    # Greedy selection: pick highest-scoring non-overlapping chunks
    selected = [all_scored.row(0, named=True)]

    for i in range(1, len(all_scored)):
        current = all_scored.row(i, named=True)
        overlaps = False

        for prev in selected:
            # Check if current overlaps with any previously selected chunk
            current_end1 = current['pos1'] + current['chunk_size'] - 1
            current_end2 = current['pos2'] + current['chunk_size'] - 1
            prev_end1 = prev['pos1'] + prev['chunk_size'] - 1
            prev_end2 = prev['pos2'] + prev['chunk_size'] - 1

            # Overlap occurs if ranges intersect
            overlap1 = not (current['pos1'] > prev_end1
                            or prev['pos1'] > current_end1)
            overlap2 = not (current['pos2'] > prev_end2
                            or prev['pos2'] > current_end2)

            if overlap1 or overlap2:
                overlaps = True
                break

        if not overlaps:
            selected.append(current)

    result = pl.DataFrame(selected)
    return result.sort('score', descending=True)


def expand_repeat(text, pos1, pos2, chunk_size, max_distance, sep):
    """Expand from seeds: start with small chunks and try to extend them."""
    char_vec = text.split(sep)
    n_chars = len(char_vec)

    # Try to extend forward
    extended_size = chunk_size

    while True:
        # Convert from R's 1-based to Python's 0-based indexing
        idx1 = pos1 - 1
        idx2 = pos2 - 1

        if idx1 + extended_size > n_chars or idx2 + extended_size > n_chars:
            break

        chunk1 = sep.join(char_vec[idx1:(idx1 + extended_size)])
        chunk2 = sep.join(char_vec[idx2:(idx2 + extended_size)])

        d = edlib.align(chunk1, chunk2)['editDistance']

        if d <= max_distance:
            extended_size = extended_size + 1
        else:
            break

    return extended_size - 1  # Return the last valid size


def find_extended_repeats(text,
                          seed_chunk_size=3,
                          max_distance=0,
                          min_gap=None,
                          sep=" "):
    """Find extended repeats starting from seed chunks."""
    # Find seed repeats
    seeds = find_internal_repeats(text,
                                  chunk_size=seed_chunk_size,
                                  max_distance=max_distance,
                                  min_gap=min_gap,
                                  sep=sep)

    if seeds is None or len(seeds) == 0:
        return None

    char_vec = text.split(sep)

    # Try to extend each seed
    extended_data = []

    for row in seeds.iter_rows(named=True):
        extended_size = expand_repeat(text,
                                      row['pos1'],
                                      row['pos2'],
                                      chunk_size=seed_chunk_size,
                                      max_distance=max_distance,
                                      sep=sep)

        # Convert from R's 1-based to Python's 0-based indexing
        idx1 = row['pos1'] - 1
        idx2 = row['pos2'] - 1

        seq1 = sep.join(char_vec[idx1:(idx1 + extended_size)])
        seq2 = sep.join(char_vec[idx2:(idx2 + extended_size)])

        extended_data.append({
            'seq1': seq1,
            'pos1': row['pos1'],
            'seq2': seq2,
            'pos2': row['pos2'],
            'gap': row['pos2'] - (row['pos1'] + extended_size),
            'lv_dist': edlib.align(seq1, seq2)['editDistance'],
            'chunk_size': extended_size
        })

    extended = pl.DataFrame(extended_data)
    extended = extended.unique()
    extended = extended.sort(['chunk_size', 'pos1'], descending=[True, False])

    return extended


# # # # # # # # # #

if __name__ == "__main__":
    # Example wordforms.
    examples = [
        "satu-satuɲa", "kəkaseh-kəkaseh", "sə-səpet", "asal-usol",
        "llama-llama", "doggy-oggy", "piggy-wiggy", "snalnal", "snalfak"
    ]

    x = examples[1]
    x = re.sub("[-]", " ", x)
    x = str_util.str_split(x)
    print(x)

    print("\n=== All repeats ===")
    reps_all = find_all_repeats(x, max_chunk_size=8, sep=" ")
    print(reps_all)

    print("\n=== Maximal (longest) repeats ===")
    reps_maximal = find_maximal_repeats(x, max_chunk_size=8, sep=" ")
    print(reps_maximal)

    print("\n=== Score-based (top 5) ===")
    reps_scored = find_best_repeats(x, max_chunk_size=8, sep=" ", top_n=5)
    print(reps_scored)

    print("\n=== Non-overlapping greedy selection ===")
    reps_nonoverlap = find_nonoverlapping_repeats(x, max_chunk_size=8, sep=" ")
    print(reps_nonoverlap)

    print("\n=== Extended from seeds ===")
    reps_extended = find_extended_repeats(x, seed_chunk_size=2, sep=" ")
    print(reps_extended)
