def local_alignment(
    query,
    target,
    match_score=2,
    mismatch_score=-1,
    gap_penalty=-2
):
    """
    Smith-Waterman-style local sequence alignment.

    Returns:
        best_score
        aligned_query
        aligned_target
    """

    query = query.upper()
    target = target.upper()

    rows = len(query) + 1
    cols = len(target) + 1

    # Dynamic programming matrix
    matrix = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    best_score = 0
    best_position = (0, 0)

    # Fill matrix
    for i in range(1, rows):
        for j in range(1, cols):

            if query[i - 1] == target[j - 1]:
                diagonal = (
                    matrix[i - 1][j - 1]
                    + match_score
                )
            else:
                diagonal = (
                    matrix[i - 1][j - 1]
                    + mismatch_score
                )

            up = matrix[i - 1][j] + gap_penalty
            left = matrix[i][j - 1] + gap_penalty

            matrix[i][j] = max(
                0,
                diagonal,
                up,
                left
            )

            if matrix[i][j] > best_score:
                best_score = matrix[i][j]
                best_position = (i, j)

    # Traceback
    i, j = best_position

    aligned_query = []
    aligned_target = []

    while i > 0 and j > 0 and matrix[i][j] > 0:

        current = matrix[i][j]

        if query[i - 1] == target[j - 1]:
            diagonal_score = (
                matrix[i - 1][j - 1]
                + match_score
            )
        else:
            diagonal_score = (
                matrix[i - 1][j - 1]
                + mismatch_score
            )

        if current == diagonal_score:

            aligned_query.append(query[i - 1])
            aligned_target.append(target[j - 1])

            i -= 1
            j -= 1

        elif current == matrix[i - 1][j] + gap_penalty:

            aligned_query.append(query[i - 1])
            aligned_target.append("-")

            i -= 1

        else:

            aligned_query.append("-")
            aligned_target.append(target[j - 1])

            j -= 1

    aligned_query.reverse()
    aligned_target.reverse()

    return (
        best_score,
        "".join(aligned_query),
        "".join(aligned_target)
    )


def calculate_alignment_identity(
    aligned_query,
    aligned_target
):
    """
    Calculate identity percentage
    from an alignment.
    """

    if not aligned_query:
        return 0.0

    matches = 0
    aligned_positions = 0

    for q, t in zip(
        aligned_query,
        aligned_target
    ):

        if q != "-" and t != "-":

            aligned_positions += 1

            if q == t:
                matches += 1

    if aligned_positions == 0:
        return 0.0

    return round(
        (matches / aligned_positions) * 100,
        2
    )


def calculate_query_coverage(
    query,
    aligned_query
):
    """
    Calculate percentage of query
    participating in the local alignment.
    """

    if len(query) == 0:
        return 0.0

    aligned_query_bases = sum(
        1 for base in aligned_query
        if base != "-"
    )

    return round(
        (aligned_query_bases / len(query)) * 100,
        2
    )
if __name__ == "__main__":

    query = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    target = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAC"

    score, aligned_query, aligned_target = local_alignment(
        query,
        target
    )

    identity = calculate_alignment_identity(
        aligned_query,
        aligned_target
    )

    coverage = calculate_query_coverage(
        query,
        aligned_query
    )

    print("\nBioBLAST-X Local Alignment")
    print("-" * 40)

    print("Alignment Score:", score)
    print("Identity:", identity, "%")
    print("Query Coverage:", coverage, "%")

    print("\nQuery :")
    print(aligned_query)

    print("\nTarget:")
    print(aligned_target)