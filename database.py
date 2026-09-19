from Bio import SeqIO


def load_fasta_database(file_path):
    """
    Load sequences from a FASTA database.
    """

    database = []

    for record in SeqIO.parse(file_path, "fasta"):
        database.append({
            "id": record.id,
            "description": record.description,
            "sequence": str(record.seq).upper()
        })

    return database


if __name__ == "__main__":
    database = load_fasta_database("data/sequences.fasta")

    print("Number of sequences:", len(database))

    for entry in database:
        print(entry["id"], "->", entry["sequence"])