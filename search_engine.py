from database import load_fasta_database
from alignment import (
    local_alignment,
    calculate_alignment_identity,
    calculate_query_coverage
)


def generate_kmers(sequence, k=4):
    """
    Generate k-mers from a sequence.
    """

    sequence = sequence.upper()
    kmers = {}

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]

        if kmer not in kmers:
            kmers[kmer] = []

        kmers[kmer].append(i)

    return kmers


def count_shared_kmers(query, target, k=4):
    """
    Count unique k-mers shared between query and target.
    """

    query_kmers = generate_kmers(query, k)
    target_kmers = generate_kmers(target, k)

    shared = set(query_kmers.keys()) & set(target_kmers.keys())

    return len(shared)


def similarity_search(query, database, k=4):
    """
    BLAST-inspired sequence similarity pipeline.

    Steps:
    1. K-mer based seed screening
    2. Local alignment
    3. Alignment score calculation
    4. Identity calculation
    5. Query coverage calculation
    6. Final ranking
    """

    query = query.upper()

    results = []

    for entry in database:

        target = entry["sequence"]

        # -----------------------------
        # STEP 1: K-mer seed screening
        # -----------------------------

        shared_kmers = count_shared_kmers(
            query,
            target,
            k
        )

        # -----------------------------
        # STEP 2: Local alignment
        # -----------------------------

        alignment_score, aligned_query, aligned_target = local_alignment(
            query,
            target
        )

        # -----------------------------
        # STEP 3: Alignment statistics
        # -----------------------------

        identity = calculate_alignment_identity(
            aligned_query,
            aligned_target
        )

        coverage = calculate_query_coverage(
            query,
            aligned_query
        )

        # -----------------------------
        # STEP 4: Final ranking score
        # -----------------------------

        final_score = (
            alignment_score
            + (shared_kmers * 0.5)
            + (identity * 0.1)
            + (coverage * 0.05)
        )

        results.append({
            "id": entry["id"],
            "description": entry["description"],
            "sequence": target,
            "shared_kmers": shared_kmers,
            "alignment_score": alignment_score,
            "identity": identity,
            "coverage": coverage,
            "final_score": round(final_score, 2),
            "aligned_query": aligned_query,
            "aligned_target": aligned_target
        })

    # Rank best hits first

    results.sort(
        key=lambda x: x["final_score"],
        reverse=True
    )

    return results


# --------------------------------------------------
# TESTING
# --------------------------------------------------

if __name__ == "__main__":

    database = load_fasta_database(
        "data/sequences.fasta"
    )

    query = (
        "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    )

    results = similarity_search(
        query,
        database,
        k=4
    )

    print("\n")
    print("=" * 70)
    print("             BioBLAST-X SIMILARITY RESULTS")
    print("=" * 70)

    for rank, result in enumerate(results, start=1):

        print(f"\nRank #{rank}")
        print("-" * 50)

        print("Sequence ID       :", result["id"])
        print("Shared k-mers     :", result["shared_kmers"])
        print("Alignment Score   :", result["alignment_score"])
        print("Identity          :", result["identity"], "%")
        print("Query Coverage    :", result["coverage"], "%")
        print("Final Score       :", result["final_score"])

        print("\nAlignment:")
        print("Query  :", result["aligned_query"])
        print("Target :", result["aligned_target"])