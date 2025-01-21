import csv
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import multiprocessing


def calculate_group_avg_length(group):
    """
    计算分组中序列的平均长度。

    Args:
        group: 序列组 [(interval1, interval2, length), ...]

    Returns:
        float: 平均长度
    """
    if not group:
        return 0
    lengths = [seq[2] for seq in group]
    return sum(lengths) / len(lengths)


def is_length_compatible(seq_length, group_avg_length, threshold=0.9):
    """
    判断序列长度是否与组平均长度相近。

    Args:
        seq_length: 当前序列长度
        group_avg_length: 组平均长度
        threshold: 允许的长度差异比例

    Returns:
        bool: 长度是否相近
    """
    if group_avg_length == 0:
        return True

    length_diff_ratio = abs(seq_length - group_avg_length) / group_avg_length
    return length_diff_ratio <= 1 - threshold


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
            try:
                current_start1, current_end1, current_start2, current_end2, length = map(int, row[1:])
                seqpair.append(((current_start1, current_end1), (current_start2, current_end2), length))
            except ValueError:
                print(f"跳过无效行: {row}")
                continue

        seqpair.sort(key=lambda x: x[2])

    return seqpair


def calculate_overlap(interval_a, interval_b):
    """
    计算两个区间的重叠比例。

    Args:
        interval_a: 第一个区间 (start, end)
        interval_b: 第二个区间 (start, end)

    Returns:
        float: 两个区间的重叠比例（0 到 1 之间的浮点数）
    """
    start_a, end_a = interval_a
    start_b, end_b = interval_b
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    overlap_length = max(0, overlap_end - overlap_start)
    total_length = min(end_a - start_a, end_b - start_b)
    return overlap_length / total_length if total_length > 0 else 0


def group_intervals_by_overlap(seqpairs):
    """
    将区间对按长度和重叠关系分组，并记录分组过程日志。

    Args:
        seqpairs: 序列对列表
        log_file: 日志文件路径

    Returns:
        list: 分组后的序列对列表
    """
    grouped_results = []



    for seqpair in seqpairs:
        interval1, interval2, length = seqpair
        added_to_group = False

        for i, group in enumerate(grouped_results):
            group_avg_length = calculate_group_avg_length(group)

            if not is_length_compatible(length, group_avg_length):
                continue


            for g in group:
                g_interval1, g_interval2, _ = g
                overlap1 = calculate_overlap(interval1, g_interval1)
                overlap2 = calculate_overlap(interval1, g_interval2)
                overlap3 = calculate_overlap(interval2, g_interval1)
                overlap4 = calculate_overlap(interval2, g_interval2)


                if (
                        overlap1 >= 0.9 or
                        overlap2 >= 0.9 or
                        overlap3 >= 0.9 or
                        overlap4 >= 0.9
                ):
                    group.append(seqpair)
                    added_to_group = True
                    break

            if added_to_group:
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

    sorted_intervals = sorted(intervals)
    merged = []
    used = set()

    for i, interval1 in enumerate(sorted_intervals):
        if i in used:
            continue

        current_start = interval1[0]
        current_end = interval1[1]

        change_made = True
        while change_made:
            change_made = False
            for j, interval2 in enumerate(sorted_intervals):
                if j in used:
                    continue

                if calculate_overlap((current_start, current_end), interval2) >= 0.9:
                    current_start = min(current_start, interval2[0])
                    current_end = max(current_end, interval2[1])
                    used.add(j)
                    change_made = True

        merged.append((current_start, current_end))

    return sorted(merged)


def merge_groups(grouped_results):
    """
    处理每个组内的所有序列。

    Args:
        grouped_results: 分组后的序列对列表

    Returns:
        list: 处理后的分组列表
    """
    merged_groups = []

    for group in grouped_results:
        all_intervals = []
        for seq in group:
            all_intervals.append(seq[0])
            all_intervals.append(seq[1])

        merged = merge_intervals_in_group(all_intervals)
        merged_groups.append(merged)

    return merged_groups


def load_genome_sequence(gbff_file):
    """
    从GBFF文件加载基因组序列。

    Args:
        gbff_file (str): GBFF文件路径

    Returns:
        str: 基因组序列，如果加载失败则返回None
    """
    try:
        for record in SeqIO.parse(gbff_file, "genbank"):
            print(f"Loading genome sequence: {record.id}")
            print(f"Sequence length: {len(record.seq)} bp")

            # 记录额外的序列信息
            if hasattr(record, 'annotations'):
                source = record.annotations.get('source', 'Unknown source')
                print(f"Source: {source}")

            return str(record.seq)

        print("No sequence found in GBFF file")
        return None

    except Exception as e:
        print(f"Error loading genome sequence: {str(e)}")
        return None


def initialize_aligner():
    """
    初始化序列比对器。

    Returns:
        PairwiseAligner: 配置好的比对器实例
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    return aligner


def compare_sequence_pair(start1, end1, start2, end2, genome):
    """
    比对两个序列并计算相似度。
    优化：将aligner的初始化移到函数内部，避免多进程共享问题
    """
    try:
        aligner = initialize_aligner()
        
        seq1 = genome[start1:end1]
        seq2 = genome[start2:end2]

        score = aligner.score(seq1, seq2)
        max_possible_score = max(len(seq1), len(seq2)) * aligner.match_score
        similarity = round(score / max_possible_score if max_possible_score > 0 else 0, 3)

        return (start1, end1, start2, end2, similarity)
    except Exception as e:
        print(f"序列比对出错: {str(e)}")
        return (start1, end1, start2, end2, 0.0)


def compare_sequences_in_groups(merged_groups, genome_seq):
    """
    对每个组内的序列进行比对。

    Args:
        merged_groups: 合并后的分组列表
        genome_seq: 基因组序列

    Returns:
        dict: {组ID: 比对结果列表}
    """
    comparison_results = {}
    
    # 获取CPU核心数并设置进程数
    cpu_count = multiprocessing.cpu_count()
    print(f"系统CPU核心数: {cpu_count}")

    # 创建一个进程池在整个处理过程中重用
    with multiprocessing.Pool() as pool:
        # 收集所有组的比对任务
        all_comparison_args = []
        group_indices = []  # 记录每个任务属于哪个组
        
        for group_id, merged_intervals in enumerate(merged_groups, start=1):
            sequences = sorted(merged_intervals, key=lambda x: x[0])
            
            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    start1, end1 = sequences[i]
                    start2, end2 = sequences[j]
                    all_comparison_args.append((start1, end1, start2, end2, genome_seq))
                    group_indices.append(group_id)

        # 批量处理所有比对任务
        if all_comparison_args:
            print(f"开始并行处理 {len(all_comparison_args)} 个序列比对任务...")
            all_results = pool.starmap(compare_sequence_pair, all_comparison_args)
            
            # 将结果按组整理
            for result, group_id in zip(all_results, group_indices):
                if group_id not in comparison_results:
                    comparison_results[group_id] = []
                comparison_results[group_id].append(result)
        
    return comparison_results


def write_results_to_csv(merged_groups, comparison_results, output_file):
    """
    将分组结果和比对结果写入CSV文件。

    Args:
        merged_groups: 合并后的分组列表
        comparison_results: 比对结果字典
        output_file: 输出文件路径
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        for group_id, merged_intervals in enumerate(merged_groups, start=1):
            writer.writerow([f'Group {group_id}'])
            writer.writerow(['Start', 'End', 'Length'])

            for start, end in sorted(merged_intervals, key=lambda x: x[0]):
                length = end - start + 1
                writer.writerow([start, end, length])

            writer.writerow([])
            writer.writerow(['Sequence Comparisons'])
            writer.writerow(['Seq1 Start', 'Seq1 End', 'Seq2 Start', 'Seq2 End', 'Similarity'])

            comparisons = comparison_results.get(group_id, [])
            for comp in comparisons:
                writer.writerow(comp)

            writer.writerow([])


def main():
    parser = argparse.ArgumentParser(description='处理基因组序列与区间比对')
    parser.add_argument('--csv-file', type=str, help='输入CSV文件路径')
    parser.add_argument('--genome', type=str, help='基因组FASTA文件路径')
    parser.add_argument('--output-file', type=str, help='输出结果CSV文件路径')

    args = parser.parse_args()

    print("开始处理序列比对...")

    seqpairs = read_intervals_from_csv(args.csv_file)
    print(f"读取到 {len(seqpairs)} 个序列对")

    grouped_results = group_intervals_by_overlap(seqpairs)
    print(f"分成 {len(grouped_results)} 组")

    merged_groups = merge_groups(grouped_results)

    genome_seq = load_genome_sequence(args.genome)
    if genome_seq is None:
        print("未能加载基因组序列。请检查GBFF文件。")
        return
    print(f"基因组序列长度: {len(genome_seq)}")

    comparison_results = compare_sequences_in_groups(merged_groups, genome_seq)
    print("完成组内序列比对")

    write_results_to_csv(merged_groups, comparison_results, args.output_file)
    print(f"\n结果已写入文件: {args.output_file}")


if __name__ == "__main__":
    main()