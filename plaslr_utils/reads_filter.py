import argparse 
import sys

from Bio import SeqIO

def filter_reads(reads_path, output, *,plasflow_path=None, plasclass_path=None, ground_truth_path=None):
    extension = reads_path.split(".")[-1]
    if extension in ["fq", "fastq"]:
        fmt = "fastq"
    else:
        fmt = "fasta"

    proper_indices = []
    if ground_truth_path:
        with open(f"{output}/filtered.fasta", "w+") as output_fasta, open(f"{output}/filtered.truth.txt", "w+") as output_truth, open(f"{output}/filtered.ids.txt", "w+") as output_ids:
            for n, (record, truth) in enumerate(zip(SeqIO.parse(reads_path, fmt), open(ground_truth_path))):
                if len(record.seq) >= 1000:
                    output_truth.write(truth)
                    output_fasta.write(f">{record.id}\n{str(record.seq)}\n")
                    output_ids.write(f"{record.id.strip()}\n")
                    proper_indices.append(n)
    else:
        with open(f"{output}/filtered.fasta", "w+") as output_fasta, open(f"{output}/filtered.ids.txt", "w+") as output_ids:
            for n, record in enumerate(SeqIO.parse(reads_path, fmt)):
                if len(record.seq) >= 1000:
                    output_fasta.write(f">{record.id}\n{str(record.seq)}\n")
                    output_ids.write(f"{record.id.strip()}\n")
                    proper_indices.append(n)
    
    # hack for faster lookup with O(1) time
    proper_indices = set(proper_indices)
    
    if plasclass_path:
        with open(f"{output}/filtered.pc", "w+") as filtered_pc:
            for n, line in enumerate(open(plasclass_path)):
                if n in proper_indices:
                    filtered_pc.write(line)

    if plasflow_path:
        with open(f"{output}/filtered.pf", "w+") as filtered_pf:
            for n, line in enumerate(open(plasflow_path)):
                if n-1 in proper_indices or n==0:
                    filtered_pf.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Filter reads and ground truth with threshold 1000bp""")

    parser.add_argument('--input', '-i',
                        help="Input reads",
                        type=str,
                        required=True)
    parser.add_argument('--plasclass', '-pc',
                        help="Plasclass result",
                        type=str,
                        default=None,
                        required=False)
    parser.add_argument('--plasflow', '-pf',
                        help="Plasclass result",
                        type=str,
                        default=None,
                        required=False)                        
    parser.add_argument('--ground-truth',
                        help="Ground truth path",
                        type=str,
                        default=None,
                        required=False)                 
    parser.add_argument('--output', '-o',
                        help="Output directory",
                        type=str,
                        required=True)    


    args = parser.parse_args()

    if args.plasclass==None and args.plasflow==None:
        print("Atleast one tool's result must be provided")
        sys.exit(1)

    filter_reads(args.input, args.output, plasflow_path=args.plasflow, plasclass_path=args.plasclass, ground_truth_path=args.ground_truth)
