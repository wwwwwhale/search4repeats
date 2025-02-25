import argparse
from collections import defaultdict
from Bio import SeqIO


def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Find and analyze repeated subsequences in genome')

    # Required arguments
    parser.add_argument('--input', '-i', required=True,
                        help='Input GBFF file path containing genome sequence')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file path for repeated subsequences')
    parser.add_argument('--length', '-l', type=int, required=True,
                        help='Length of subsequences to search for')

    return parser.parse_args()


def read_genome_sequence(gbff_file):
    """
    Read genome sequence from GBFF file

    Args:
        gbff_file (str): Path to GBFF file

    Returns:
        str: Genome sequence
    """
    try:
        # Parse GBFF file using BioPython's GenBank parser
        for record in SeqIO.parse(gbff_file, "genbank"):
            # Return the first sequence found (assuming single sequence per file)
            return str(record.seq).strip().upper()
    except Exception as e:
        print(f"Error reading GBFF file: {e}")
        return None
    return None


def find_repeated_subsequences(sequence, sub_len):
    """
    Find repeated subsequences of specified length in sequence

    Args:
        sequence (str): Input sequence to search
        sub_len (int): Length of subsequences to find

    Returns:
        list: List of tuples containing (subsequence, count, positions)
    """
    seq_len = len(sequence)
    seen = defaultdict(list)
    repeated_subsequences = []

    for i in range(seq_len - sub_len + 1):
        sub_seq = sequence[i:i + sub_len]
        seen[sub_seq].append(i)

    for subseq, positions in seen.items():
        if len(positions) > 1:
            repeated_subsequences.append((subseq, len(positions), positions))

    return repeated_subsequences


def get_merged_subsequence(sequence, merged_positions):
    """
    Extract subsequence based on merged positions

    Args:
        sequence (str): Input sequence
        merged_positions (list): List of merged position tuples

    Returns:
        str: Extracted subsequence
    """
    start, end = merged_positions[0]
    return sequence[start:end + 1]


def compare_and_group_subsequences(repeated_subsequences, length):
    """
    Compare and group subsequences based on repeat patterns

    Args:
        repeated_subsequences (list): List of subsequence tuples
        length (int): Length of subsequences

    Returns:
        list: List of grouped subsequences
    """
    repeated_subsequences.sort(key=lambda x: x[1], reverse=True)
    grouped_subsequences = []
    current_group = []

    for i in range(len(repeated_subsequences)):
        subseq, count, positions = repeated_subsequences[i]

        if not current_group:
            current_group.append((subseq, count, positions))
            continue

        last_subseq, last_count, last_positions = current_group[-1]

        if count == last_count:
            can_merge = True
            distance = None
            for pos1, pos2 in zip(positions, last_positions):
                if abs(pos1 - pos2) <= length:
                    if distance is None:
                        distance = pos2 - pos1
                    elif pos2 - pos1 != distance:
                        can_merge = False
                        break
                else:
                    can_merge = False
                    break

            if can_merge:
                current_group.append((subseq, count, positions))
            else:
                grouped_subsequences.append(current_group)
                current_group = [(subseq, count, positions)]
        else:
            grouped_subsequences.append(current_group)
            current_group = [(subseq, count, positions)]

    if current_group:
        grouped_subsequences.append(current_group)

    return grouped_subsequences


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Read genome sequence from GBFF file
    genome_sequence = read_genome_sequence(args.input)
    if not genome_sequence:
        print(f"Error: Could not read genome sequence from {args.input}")
        return

    # Find repeated subsequences
    repeated_subsequences = find_repeated_subsequences(genome_sequence, args.length)

    # Group subsequences
    grouped_subsequences = compare_and_group_subsequences(repeated_subsequences, args.length)

    # Process and write results
    with open(args.output, 'w', encoding='utf-8') as output_file:
        merged_repeated_subsequences = []

        for group in grouped_subsequences:
            first_subseq, _, first_positions = group[0]
            last_subseq, _, last_positions = group[-1]

            merged_positions = [
                (start, end + len(last_subseq) - 1)
                for start, end in zip(first_positions, last_positions)
            ]

            merged_subseq = get_merged_subsequence(genome_sequence, merged_positions)
            count = group[0][1]
            merged_repeated_subsequences.append((merged_subseq, count, merged_positions))

        for subseq, count, positions in merged_repeated_subsequences:
            positions_str = ", ".join([f"[{pos[0]},{pos[1]}]" for pos in positions])
            output_file.write(f"重复出现次数：{count}\n重复序列：{subseq}\nPositions: {positions_str}\n\n")


if __name__ == '__main__':
    main()