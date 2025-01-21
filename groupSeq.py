import csv
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import multiprocessing
import logging  # 导入日志模块

# 设置日志配置，将日志输出到文件
logging.basicConfig(level=logging.INFO,  # 设置日志级别为INFO
                    format='%(asctime)s - %(levelname)s - %(message)s',  # 设置日志格式
                    handlers=[logging.FileHandler('process.log', mode='w', encoding='utf-8')])  # 输出到文件


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
    avg_length = sum(lengths) / len(lengths)
    logging.debug(f"Calculated average length: {avg_length}")
    return avg_length


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
    compatible = length_diff_ratio <= 1 - threshold
    logging.debug(f"Checking length compatibility: {seq_length} vs {group_avg_length}, compatible: {compatible}")
    return compatible


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
                logging.warning(f"跳过无效行: {row}")
                continue

        seqpair.sort(key=lambda x: x[2])

    logging.info(f"Loaded {len(seqpair)} sequence pairs from CSV.")
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
    overlap_ratio = overlap_length / total_length if total_length > 0 else 0
    logging.debug(f"Calculated overlap: {overlap_ratio} between intervals {interval_a} and {interval_b}")
    return overlap_ratio


def group_intervals_by_overlap(seqpairs):
    """
    将区间对按长度和重叠关系分组，并记录分组过程日志。

    Args:
        seqpairs: 序列对列表

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
                    logging.info(f"Added sequence pair {seqpair} to group {i + 1}")
                    break

            if added_to_group:
                break

        if not added_to_group:
            grouped_results.append([seqpair])
            logging.info(f"Created new group for sequence pair {seqpair}")

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
            logging.info(f"Loading genome sequence: {record.id}")
            logging.info(f"Sequence length: {len(record.seq)} bp")

            # 记录额外的序列信息
            if hasattr(record, 'annotations'):
                source = record.annotations.get('source', 'Unknown source')
                logging.info(f"Source: {source}")

            return str(record.seq)

        logging.error("No sequence found in GBFF file")
        return None

    except Exception as e:
        logging.error(f"Error loading genome sequence: {str(e)}")
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
    logging.info("Aligner initialized")
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
        logging.error(f"序列比对出错: {str(e)}")
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
    pool = multiprocessing.Pool(cpu_count)

    all_comparison_args = []
    group_indices = []

    for group_id, group in enumerate(merged_groups, start=1):
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                interval1 = group[i]
                interval2 = group[j]
                all_comparison_args.append(
                    (interval1[0], interval1[1], interval2[0], interval2[1], genome_seq)
                )
                group_indices.append(group_id)

        # 使用pool.map进行并行比对
        results = pool.starmap(compare_sequence_pair, all_comparison_args)

        # 将比对结果按照组ID归类
        for result, group_id in zip(results, group_indices):
            if group_id not in comparison_results:
                comparison_results[group_id] = []
            comparison_results[group_id].append(result)

    return comparison_results


def main():
    parser = argparse.ArgumentParser(description="处理序列对")
    parser.add_argument("-i", "--input", required=True, help="输入的CSV文件路径")
    parser.add_argument("-g", "--genome", required=True, help="基因组GBFF文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出文件路径")
    args = parser.parse_args()

    # 读取CSV文件
    seqpairs = read_intervals_from_csv(args.input)

    # 加载基因组序列
    genome_sequence = load_genome_sequence(args.genome)
    if genome_sequence is None:
        return

    # 按照重叠关系分组
    grouped_results = group_intervals_by_overlap(seqpairs)

    # 合并组内的重叠序列区间
    merged_groups = merge_groups(grouped_results)

    # 对每个组进行序列比对
    comparison_results = compare_sequences_in_groups(merged_groups, genome_sequence)

    # 将比对结果输出到文件
    with open(args.output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Group", "Start1", "End1", "Start2", "End2", "Similarity"])
        for group_id, results in comparison_results.items():
            for result in results:
                writer.writerow([group_id, *result])

    logging.info(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()