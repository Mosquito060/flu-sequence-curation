from Bio import SeqIO
from collections import Counter
import csv
import argparse

def normalize_sequence(seq):
    return str(seq).upper().replace("\n", "").replace(" ", "")

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    seq_dict = {}

    for record in SeqIO.parse(args.input, "fasta"):

        seq = normalize_sequence(record.seq)

        if seq not in seq_dict:
            seq_dict[seq] = []

        seq_dict[seq].append(record.id)

    total = sum(len(v) for v in seq_dict.values())

    fasta_out = args.output
    csv_out = args.output.replace(".fasta", "_counts.csv")

    with open(fasta_out, "w") as fout, open(csv_out, "w", newline="") as csvfile:

        writer = csv.writer(csvfile)
        writer.writerow(["Variant_ID", "Count", "Frequency"])

        for i, (seq, ids) in enumerate(seq_dict.items(), start=1):

            count = len(ids)
            freq = count / total
            variant_id = f"HAP_{i}"

            header = f"{variant_id}|{ids[0]}|count={count}|freq={freq:.4f}"

            fout.write(f">{header}\n{seq}\n")

            writer.writerow([variant_id, count, round(freq, 6)])

    print(f"Total sequences: {total}")
    print(f"Unique variants: {len(seq_dict)}")

if __name__ == "__main__":
    main()