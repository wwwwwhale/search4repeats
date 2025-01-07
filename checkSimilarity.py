import argparse
import csv
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from collections import defaultdict


def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Analyze sequence similarities between groups')

    # Required arguments
    parser.add_argument('--genome', '-g', required=True,
                        help='Input FASTA file path containing genome sequence')
    parser.add_argument('--intervals', '-i', required=True,
                        help='Input CSV file containing merged intervals')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file path for similarity analysis')

    # Optional alignment parameters
    parser.add_argument('--match-score', type=float, default=2,
                        help='Score for matching bases (default: 2)')
    parser.add_argument('--mismatch-score', type=float, default=-1,
                        help='Score for mismatching bases (default: -1)')
    parser.add_argument('--open-gap-score', type=float, default=-0.5,
                        help='Penalty for opening a gap (default: -0.5)')
    parser.add_argument('--extend-gap-score', type=float, default=-0.1,
                        help='Penalty for extending a gap (default: -0.1)')

    return parser.parse_args()


def load_genome_sequence(fasta_file):
    """
    Load genome sequence from FASTA file.

    Args:
        fasta_file (str): Path to FASTA file

    Returns:
        str: Genome sequence
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq)
    return None


import csv
from collections import defaultdict


def load_intervals(csv_file, reference_genome):
    """
    Load and validate intervals from CSV file.

    Args:
        csv_file (str): Path to CSV file containing intervals
        reference_genome (str): Reference genome sequence

    Returns:
        dict: Dictionary of group sequences, with keys starting from 1
    """
    group_sequences = defaultdict(list)

    with open(csv_file, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)

        group_id = None
        for row in reader:
            # Skip empty rows
            if not row:
                continue

            # Detect group header (e.g., 'Group 1')
            if row[0].startswith("Group"):
                group_id = int(row[0].split()[1])  # Extract group number (e.g., 1 for 'Group 1')
                continue  # Skip the 'Start End Length' header row

            # Read interval data after group header
            if group_id is not None:
                try:
                    start = int(row[0])
                    end = int(row[1])
                    length = int(row[2])

                    # Validate interval
                    if start < 1 or end > len(reference_genome):
                        print(f"Warning: Interval {start}-{end} out of genome bounds. Skipping.")
                        continue

                    # Extract sequence from reference genome
                    seq = reference_genome[start - 1:end]
                    group_sequences[group_id].append({
                        'start': start,
                        'end': end,
                        'sequence': seq,
                        'length': length
                    })
                    print(f"Extracted sequence for group {group_id}: length {len(seq)}")

                except ValueError:
                    # Handle invalid row (could happen due to empty values or malformed data)
                    print(f"Skipping invalid row: {row}")

    return group_sequences


def setup_aligner(args):
    """
    Configure PairwiseAligner with provided parameters.

    Args:
        args (Namespace): Command line arguments

    Returns:
        PairwiseAligner: Configured aligner
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match = args.match_score
    aligner.mismatch = args.mismatch_score
    aligner.open_gap_score = args.open_gap_score
    aligner.extend_gap_score = args.extend_gap_score
    return aligner


def analyze_sequences(group_sequences, aligner, output_file):
    """
    Analyze sequence similarities and write results to file.

    Args:
        group_sequences (dict): Dictionary of group sequences
        aligner (PairwiseAligner): Configured aligner
        output_file (str): Path to output file
    """
    similarities = {}

    with open(output_file, 'w') as out:
        for group, sequences in group_sequences.items():
            if len(sequences) < 2:
                continue

            out.write(f"\n{'=' * 80}\n")
            out.write(f"Group {group} Analysis:\n")
            out.write(f"Number of sequences in group: {len(sequences)}\n")

            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    seq1 = sequences[i]['sequence']
                    seq2 = sequences[j]['sequence']
                    pair_id = f"Group{group}_Seq{i + 1}_vs_Seq{j + 1}"
                    print(f"Comparing {pair_id}")

                    out.write(f"\nComparing sequences in {pair_id}:\n")
                    out.write(f"Sequence 1 position: {sequences[i]['start']}-{sequences[i]['end']} "
                              f"(length: {sequences[i]['length']})\n")
                    out.write(f"Sequence 2 position: {sequences[j]['start']}-{sequences[j]['end']} "
                              f"(length: {sequences[j]['length']})\n")

                    try:
                        alignments = aligner.align(seq1, seq2)
                        score = alignments.score
                        max_score = max(len(seq1), len(seq2)) * aligner.match_score
                        similarity_percentage = (score / max_score) * 100 if max_score != 0 else 0.0

                        out.write(f"\nSimilarity: {similarity_percentage:.2f}%\n")
                        similarities[pair_id] = (score, similarity_percentage)
                    except Exception as e:
                        out.write(f"Error during alignment: {str(e)}\n")

            # Write group summary
            out.write(f"\nGroup {group} Similarity Scores Summary:\n")
            group_similarities = {k: v for k, v in similarities.items() if k.startswith(f"Group{group}")}
            for pair, (score, percentage) in group_similarities.items():
                out.write(f"{pair}: Score = {score:.4f}, Similarity = {percentage:.2f}%\n")

            out.write(f"\nGroup {group} analysis complete\n")


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Load genome sequence
    reference_genome = load_genome_sequence(args.genome)
    if not reference_genome:
        print(f"Error: Could not load genome sequence from {args.genome}")
        return

    # Load and validate intervals
    group_sequences = load_intervals(args.intervals, reference_genome)
    if not group_sequences:
        print("Error: No valid intervals found in CSV file")
        return

    # Setup aligner with provided parameters
    aligner = setup_aligner(args)

    # Perform analysis
    analyze_sequences(group_sequences, aligner, args.output)

    print("Analysis complete. Results written to", args.output)


if __name__ == '__main__':
    main()