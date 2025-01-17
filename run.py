import os
import sys
import subprocess


def run_command(command, error_message, output_file=None):
    try:
        result = subprocess.run(command, check=True, shell=True, capture_output=True, text=True)
        
        # 检查输出文件是否存在且非空
        if output_file:
            if not os.path.exists(output_file):
                raise subprocess.CalledProcessError(1, command, 
                    f"Output file '{output_file}' was not created")
            if os.path.getsize(output_file) == 0:
                raise subprocess.CalledProcessError(1, command, 
                    f"Output file '{output_file}' is empty")
                
        return True
    except subprocess.CalledProcessError as e:
        print(f"{error_message}", file=sys.stderr)
        print(f"Command output:", file=sys.stderr)
        if e.stdout:
            print(e.stdout, file=sys.stderr)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        sys.exit(1)



def main():
    if len(sys.argv) < 2:
        print("Usage: python run.py <genome_file(.gbff)> [output_dir]", file=sys.stderr)
        sys.exit(1)

    genome_file = sys.argv[1]

    # 检查文件是否存在
    for file_path in [genome_file]:
        if not os.path.isfile(file_path):
            print(f"Error: File '{file_path}' does not exist.", file=sys.stderr)
            sys.exit(1)

    # 获取输出目录（可选参数），默认为当前目录
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    # 检查输出目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        print(f"Output directory '{output_dir}' does not exist. Creating it...")
        os.makedirs(output_dir)

    # 动态生成输出文件名
    base_name = os.path.splitext(os.path.basename(genome_file))[0]
    repeated_seq_file = os.path.join(output_dir, f"{base_name}_repeated_seq.txt")
    extended_file = os.path.join(output_dir, f"{base_name}_similar_seq.csv")
    group_file = os.path.join(output_dir, f"{base_name}_grouped_results.csv")
    summary_file = os.path.join(output_dir, f"{base_name}_summary.csv")
    repeated_len = 100

    # Step 1: 运行 findRepeatedSeq.py
    print("Running findRepeatedSeq.py...")
    find_repeated_seq_cmd = (
        f"python findRepeatedSeq.py -i {genome_file} -o {repeated_seq_file} -l {repeated_len}"
    )
    run_command(find_repeated_seq_cmd, "Error: findRepeatedSeq.py failed.", repeated_seq_file)

    # 检查文件是否存在且非空
    if not os.path.exists(repeated_seq_file) or os.path.getsize(repeated_seq_file) == 0:
        print(f"Error: {repeated_seq_file} was not created or is empty", file=sys.stderr)
        sys.exit(1)

    # Step 2: 运行 maxSimilarBoundary.py
    print("Running maxSimilarBoundary.py...")
    max_similar_boundary_cmd = (
        f"python maxSimilarBoundary.py -i {repeated_seq_file} -g {genome_file} -o {extended_file}"
    )
    run_command(max_similar_boundary_cmd, "Error: maxSimilarBoundary.py failed.", extended_file)

    # 检查文件是否存在且非空
    if not os.path.exists(extended_file) or os.path.getsize(extended_file) == 0:
        print(f"Error: {extended_file} was not created or is empty", file=sys.stderr)
        sys.exit(1)

    # Step 3: 运行 groupSeq.py
    print("Running groupSeq.py...")
    group_seq_cmd = (
        f"python groupSeq.py --csv-file {extended_file} "
        f"--genome {genome_file} --output-file {group_file}"
    )
    run_command(group_seq_cmd, "Error: groupSeq.py failed.", group_file)

    # Step 4: 运行序列注释（如果有GFF文件）
    if genome_file and os.path.isfile(genome_file):
        print("Running sequence annotation...")
        annotation_cmd = (
            f"python seqAnnotation.py -i {group_file} -g {genome_file} "
        )
        run_command(annotation_cmd, "Error: annotateSequences.py failed.")
        print(f"Annotation results saved to {group_file}")
    else:
        print("Skipping annotation step: GFF file not found")

    # Step 5: 运行序列注释（如果有GFF文件）
    print("Running summary.py...")
    group_seq_cmd = (
        f"python summary.py --input {group_file} "
        f"--output {summary_file}"
    )
    run_command(group_seq_cmd, "Error: summary.py failed.")
    print(f"Output saved to {summary_file}")

    print("All steps completed successfully!")


if __name__ == "__main__":
    main()
