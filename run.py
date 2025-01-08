import os
import sys
import subprocess


def run_command(command, error_message):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        print(error_message, file=sys.stderr)
        sys.exit(1)


def main():
    # 检查是否提供了基因组文件作为命令行参数
    if len(sys.argv) < 2:
        print("Usage: python run.py <genome_file>", file=sys.stderr)
        sys.exit(1)

    genome_file = sys.argv[1]

    # 检查基因组文件是否存在
    if not os.path.isfile(genome_file):
        print(f"Error: Genome file '{genome_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # 动态生成输出文件名
    base_name = os.path.splitext(genome_file)[0]
    repeated_seq_file = f"{base_name}_repeated_seq.txt"
    repeated_len = 500
    extended_file = f"{base_name}_similar_seq.csv"
    group_file = f"{base_name}_grouping_results.csv"
    similarity_file = f"{base_name}_check_similarity.txt"

    # Step 1: 运行 findRepeatedSeq.py
    print("Running findRepeatedSeq.py...")
    find_repeated_seq_cmd = (
        f"python findRepeatedSeq.py -i {genome_file} -o {repeated_seq_file} -l {repeated_len}"
    )
    run_command(find_repeated_seq_cmd, "Error: findRepeatedSeq.py failed.")
    print(f"Output saved to {repeated_seq_file}")

    # Step 2: 运行 maxSimilarBoundary.py
    print("Running maxSimilarBoundary.py...")
    max_similar_boundary_cmd = (
        f"python maxSimilarBoundary.py -i {repeated_seq_file} -g {genome_file} -o {extended_file}"
    )
    run_command(max_similar_boundary_cmd, "Error: maxSimilarBoundary.py failed.")
    print(f"Output saved to {extended_file}")

    # Step 3: 运行 mergeSeq.py
    print("Running mergeSeq.py...")
    merge_seq_cmd = (
        f"python mergeSeq.py --interval-file {extended_file} "
        f"--sequence-file {repeated_seq_file} --output-file {group_file}"
    )
    run_command(merge_seq_cmd, "Error: mergeSeq.py failed.")
    print(f"Output saved to {group_file}")

    # Step 4: 运行 checkSimilarity.py
    print("Running checkSimilarity.py...")
    check_similarity_cmd = (
        f"python checkSimilarity.py -g {genome_file} -i {group_file} -o {similarity_file}"
    )
    run_command(check_similarity_cmd, "Error: checkSimilarity.py failed.")
    print(f"Output saved to {similarity_file}")

    print("All steps completed successfully!")


if __name__ == "__main__":
    main()
