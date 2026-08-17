import argparse
import re
from Bio import SeqIO

def calculate_ambiguity(seq):
    return seq.count("X") + seq.count("Z") + seq.count("*")

def ambiguity_percentage(seq):
    return (calculate_ambiguity(seq) / len(seq)) * 100 if len(seq) > 0 else 100

def check_start_motif(seq, motif):
    return seq.startswith(motif)

def check_regex_pattern(seq, pattern):
    if not pattern:
        return True
    return re.search(pattern, seq) is not None

def check_ha1_coverage(seq, min_length):
    return len(seq) >= min_length

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="qc_filtered.fasta")
    parser.add_argument("--qc_report", default="qc_report.csv")
    parser.add_argument("--min_length", type=int, default=550)
    parser.add_argument("--max_ambiguity", type=float, default=5.0)
    parser.add_argument("--start_motif", default="MK")
    parser.add_argument("--regex", default=None)

    args = parser.parse_args()

    with open(args.output, "w") as fout, open(args.qc_report, "w") as qc:

        qc.write("ID,Status,Reason,Length,Ambiguity\n")

        for record in SeqIO.parse(args.input, "fasta"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            amb = ambiguity_percentage(seq)

            reason = None

            if not check_ha1_coverage(seq, args.min_length):
                reason = "Coverage"
            elif amb > args.max_ambiguity:
                reason = "Ambiguity"
            elif not check_start_motif(seq, args.start_motif):
                reason = "Motif"
            elif not check_regex_pattern(seq, args.regex):
                reason = "Regex"

            if reason:
                qc.write(f"{record.id},EXCLUDED,{reason},{seq_len},{amb:.2f}\n")
                continue

            SeqIO.write(record, fout, "fasta")
            qc.write(f"{record.id},KEPT,OK,{seq_len},{amb:.2f}\n")

if __name__ == "__main__":
    main()
