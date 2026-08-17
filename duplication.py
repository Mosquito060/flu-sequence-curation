import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Remoção de sequências duplicadas (exatas)")
    
    parser.add_argument("-i", "--input", required=True, help="Arquivo FASTA de entrada")
    parser.add_argument("-o", "--output", default="deduplicated.fasta", help="FASTA sem duplicatas")
    parser.add_argument("--report", default="duplicates_report.csv", help="Relatório de duplicatas")

    args = parser.parse_args()

    seen = {}
    
    with open(args.output, "w") as fout, open(args.report, "w") as rep:
        rep.write("ID,Status,Representative\n")

        for record in SeqIO.parse(args.input, "fasta"):
            seq = str(record.seq).upper()

            if seq in seen:
                rep.write(f"{record.id},DUPLICATE,{seen[seq]}\n")
            else:
                seen[seq] = record.id
                SeqIO.write(record, fout, "fasta")
                rep.write(f"{record.id},KEPT,{record.id}\n")

    print(f"Arquivo gerado: {args.output}")
    print(f"Relatório: {args.report}")

if __name__ == "__main__":
    main()