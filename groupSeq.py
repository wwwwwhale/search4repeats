import csv
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import multiprocessing

def read_intervals_from_csv(csv_file_path):
    """
    读取CSV文件，返回区间对和length。

    Returns:
        list: [((start1, end1), (start2, end2), length), ...]
    """
    seqpair = []
    with open(csv_file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # 跳过标题行

        for row in reader:
            if not row:
                continue
            # 假设CSV的格式为：ID, start1, end1, start2, end2, length
            try:
                current_start1, current_end1, current_start2, current_end2, length = map(int, row[1:])
                seqpair.append(((current_start1, current_end1), (current_start2, current_end2), length))
            except ValueError:
                print(f"跳过无效行: {row}")
                continue

        # 按照length值进行排序
        seqpair.sort(key=lambda x: x[2])

    return seqpair


def calculate_overlap(interval_a, interval_b):
    """
    计算两个区间的重叠比例。
    :param interval_a: 第一个区间 (start, end)
    :param interval_b: 第二个区间 (start, end)
    :return: 如果重叠比例 >= 80%，返回 True，否则返回 False
    """
    start_a, end_a = interval_a
    start_b, end_b = interval_b
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    overlap_length = max(0, overlap_end - overlap_start)  # 重叠长度
    total_length = max(end_a - start_a, end_b - start_b)  # 两区间的最大长度
    return (overlap_length / total_length) >= 0.80



def group_intervals_by_overlap(seqpairs):
    """
    将区间对按重叠关系分组。

    Returns:
        list: 分组后的序列对列表
    """
    grouped_results = []

    for seqpair in seqpairs:
        interval1, interval2, length = seqpair
        added_to_group = False

        for group in grouped_results:
            # 检查当前序列对是否与组内任何序列对有重叠
            if any(
                calculate_overlap(interval1, g[0]) >= 0.85 or
                calculate_overlap(interval1, g[1]) >= 0.85 or
                calculate_overlap(interval2, g[0]) >= 0.85 or
                calculate_overlap(interval2, g[1]) >= 0.85
                for g in group
            ):
                group.append(seqpair)
                added_to_group = True
                break

        if not added_to_group:
            grouped_results.append([seqpair])

    return grouped_results

def merge_intervals_in_group(intervals):
    """
    合并组内重叠的序列区间。

    Args:
        intervals: 需要合并的区间列表 [(start, end), ...]

    Returns:
        list: 合并后的区间列表 [(start, end), ...]
    """
    if not intervals:
        return []

    # 按起始位置排序
    sorted_intervals = sorted(intervals)
    merged = []
    used = set()

    for i, interval1 in enumerate(sorted_intervals):
        if i in used:
            continue

        current_start = interval1[0]
        current_end = interval1[1]

        # 不断检查和合并重叠窗口
        change_made = True
        while change_made:
            change_made = False
            for j, interval2 in enumerate(sorted_intervals):
                if j in used:
                    continue

                if calculate_overlap((current_start, current_end), interval2) >= 0.85:
                    current_start = min(current_start, interval2[0])
                    current_end = max(current_end, interval2[1])
                    used.add(j)
                    change_made = True

        merged.append((current_start, current_end))

    return sorted(merged)

def merge_groups(grouped_results):
    """
    处理每个组内的所有序列，不区分第一序列和第二序列。

    Args:
        grouped_results: 分组后的序列对列表

    Returns:
        list: 处理后的分组列表，每组包含合并后的序列
    """
    merged_groups = []

    for group in grouped_results:
        # 收集组内所有序列
        all_intervals = []
        for seq in group:
            all_intervals.append(seq[0])  # 添加第一个序列
            all_intervals.append(seq[1])  # 添加第二个序列

        # 合并该组的所有序列
        merged = merge_intervals_in_group(all_intervals)
        merged_groups.append(merged)

    return merged_groups

def load_genome_sequence(fasta_file):
    """
    加载基因组序列。

    参数:
    - fasta_file (str): FASTA文件路径

    返回:
    - str: 基因组序列
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_seq = str(record.seq)
        return genome_seq

def initialize_aligner():
    """
    初始化PairwiseAligner对象，每个进程独立拥有自己的对齐器实例。

    返回:
    - PairwiseAligner: 配置好的对齐器
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # 使用全局比对模式
    aligner.match_score = 2  # 匹配的得分
    aligner.mismatch_score = -1  # 错配的惩罚
    aligner.open_gap_score = -0.5  # 打开间隙的惩罚
    aligner.extend_gap_score = -0.1  # 延长间隙的惩罚
    return aligner

def compare_sequence_pair(start1, end1, start2, end2, genome):
    """
    比对两个序列并计算相似度。

    参数:
    - start1, end1: 第一个序列的起始和结束位置
    - start2, end2: 第二个序列的起始和结束位置
    - genome (str): 基因组序列

    返回:
    - tuple: (Seq1 Start, Seq1 End, Seq2 Start, Seq2 End, Similarity)
    """
    aligner = initialize_aligner()

    seq1 = genome[start1:end1]
    seq2 = genome[start2:end2]

    score = aligner.score(seq1, seq2)
    max_possible_score = len(seq1) * aligner.match_score
    similarity = round(score / max_possible_score if max_possible_score > 0 else 0, 3)

    return (start1, end1, start2, end2, similarity)

def compare_sequences_in_groups(merged_groups, genome_seq):
    """
    对每个组内的序列进行比对，并返回比对结果。

    Args:
        merged_groups: 合并后的分组列表
        genome_seq: 基因组序列

    Returns:
        dict: {组ID: 比对结果列表}
    """
    comparison_results = {}

    for group_id, merged_intervals in enumerate(merged_groups, start=1):
        # 收集组内所有序列
        sequences = sorted(merged_intervals, key=lambda x: x[0])

        # 生成所有可能的序列对进行比对
        comparison_args = []
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                start1, end1 = sequences[i]
                start2, end2 = sequences[j]
                comparison_args.append((start1, end1, start2, end2, genome_seq))

        # 使用多进程加速比对
        if comparison_args:
            with multiprocessing.Pool() as pool:
                comparisons = pool.starmap(compare_sequence_pair, comparison_args)
            comparison_results[group_id] = comparisons
        else:
            comparison_results[group_id] = []

    return comparison_results


def write_results_to_csv(merged_groups, comparison_results, output_file: str):
    """
    将分组结果和比对结果写入CSV文件。

    Args:
        merged_groups: 合并后的分组列表
        comparison_results: 比对结果字典 {组ID: 比对结果列表}
        output_file: 输出文件路径
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 按组ID排序并写入，从1开始编号
        for group_id, merged_intervals in enumerate(merged_groups, start=1):
            writer.writerow([f'Group {group_id}'])
            writer.writerow(['Start', 'End', 'Length'])

            # 写入该组的所有区间
            for start, end in sorted(merged_intervals, key=lambda x: x[0]):
                length = end - start + 1
                writer.writerow([start, end, length])

            writer.writerow([])

            # 写入组内序列比对结果
            writer.writerow(['Sequence Comparisons'])
            writer.writerow(['Seq1 Start', 'Seq1 End', 'Seq2 Start', 'Seq2 End', 'Similarity'])

            comparisons = comparison_results.get(group_id, [])
            for comp in comparisons:
                s1, e1, s2, e2, similarity = comp
                writer.writerow([s1, e1, s2, e2, similarity])

            writer.writerow([])

def main():
    parser = argparse.ArgumentParser(description='处理基因组序列与区间比对')

    parser.add_argument('--csv-file', type=str, help='输入CSV文件路径')
    parser.add_argument('--genome', type=str, help='基因组FASTA文件路径')
    parser.add_argument('--output-file', type=str, help='输出结果CSV文件路径')

    args = parser.parse_args()

    # 1. 读取序列对
    seqpairs = read_intervals_from_csv(args.csv_file)
    print(f"读取到 {len(seqpairs)} 个序列对")

    # 2. 分组
    grouped_results = group_intervals_by_overlap(seqpairs)
    print(f"分成 {len(grouped_results)} 组")

    # 3. 合并每组内的序列
    merged_groups = merge_groups(grouped_results)

    # 4. 加载基因组序列
    genome_seq = load_genome_sequence(args.genome)
    if genome_seq is None:
        print("未能加载基因组序列。请检查FASTA文件。")
        return
    print(f"基因组序列长度: {len(genome_seq)}")

    # 5. 比对组内序列
    comparison_results = compare_sequences_in_groups(merged_groups, genome_seq)
    print("完成组内序列比对")

    # 6. 输出结果到CSV文件，并包含比对结果
    write_results_to_csv(merged_groups, comparison_results, args.output_file)
    print(f"\n结果已写入文件: {args.output_file}")

if __name__ == "__main__":
    main()
