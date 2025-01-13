import csv
import re
from collections import defaultdict, namedtuple
import argparse
import numpy as np
from typing import List, Dict, NamedTuple, Optional, Set
from collections import defaultdict
import maxSimilarBoundary
from Bio import SeqIO

class IntervalOperations:
    """区间操作类，处理区间的并集和交集运算"""

    def __init__(self):
        self.Interval = namedtuple('Interval', ['start', 'end'])

    def is_overlapping(self, interval1: namedtuple, interval2: namedtuple) -> bool:
        """检查两个区间是否重叠"""
        return not (interval1.end <= interval2.start or interval2.end <= interval1.start)

    def read_intervals_from_csv(self, file_path: str, start_col1: int, end_col1: int,
                                start_col2: Optional[int] = None, end_col2: Optional[int] = None,
                                skip_header: bool = True) -> List[namedtuple]:
        """从CSV文件读取区间数据"""
        intervals = []
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file)
            if skip_header:
                next(reader)

            for row in reader:
                intervals.append(self.Interval(int(row[start_col1]), int(row[end_col1])))
                if start_col2 is not None and end_col2 is not None:
                    intervals.append(self.Interval(int(row[start_col2]), int(row[end_col2])))

        return sorted(intervals, key=lambda x: x.start)

    def merge_intervals(self, intervals: List[namedtuple]) -> List[namedtuple]:
        """合并重叠度超过80%的区间"""
        if not intervals:
            return []

        # 按起始位置排序
        intervals.sort(key=lambda x: x.start)
        merged = [intervals[0]]

        for current in intervals[1:]:
            last = merged[-1]

            # 计算重叠部分
            if current.start <= last.end:
                overlap_length = min(last.end, current.end) - current.start
                last_length = last.end - last.start
                current_length = current.end - current.start

                # 计算重叠比例(取两个区间中较小的那个作为基准)
                overlap_ratio = overlap_length / min(last_length, current_length)
                # # 计算重叠比例(取两个区间中较小的那个作为基准)
                # overlap_ratio = overlap_length / max(last_length, current_length)

                # 重叠度超过80%才合并
                if overlap_ratio > 0.8:
                    merged[-1] = self.Interval(
                        last.start,
                        max(last.end, current.end)
                    )
                else:
                    merged.append(current)
            else:
                merged.append(current)

        return merged

    def find_intersections(self, intervals: List[namedtuple]) -> List[namedtuple]:
        """找出区间的交集"""
        if not intervals:
            return []

        points = []
        for interval in intervals:
            points.append((interval.start, 'start'))
            points.append((interval.end, 'end'))
        points.sort()

        intersections = []
        count = 0
        intersection_start = None

        for pos, point_type in points:
            if point_type == 'start':
                count += 1
                if count >= 2 and intersection_start is None:
                    intersection_start = pos
            else:
                count -= 1
                if count < 2 and intersection_start is not None:
                    intersections.append(self.Interval(intersection_start, pos))
                    intersection_start = None

        return intersections

    def write_intervals_to_csv(self, intervals: List[namedtuple], output_file: str,
                               operation_type: str = 'merged'):
        """将区间结果写入CSV文件"""
        with open(output_file, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            header = [f'{operation_type}_start', f'{operation_type}_end']
            writer.writerow(header)

            for interval in intervals:
                writer.writerow([interval.start, interval.end])

def load_genome_sequence(fasta_file: str) -> str:
    """加载基因组序列"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_seq = str(record.seq)
        return genome_seq

def extract_sequence(genome_seq: str, start: int, end: int) -> str:
    """从基因组中提取指定区间的序列"""
    if not genome_seq:
        return ""
    return genome_seq[start:end]

def parse_sequences(file_path):
    """解析TXT文件，提取重复序列信息"""
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
            print("文件内容读取成功")

        count_pattern = r"重复出现次数：(\d+)"
        sequence_pattern = r"重复序列：([ATGC]+)"
        positions_pattern = r"Positions:\s*((?:\[\d+,\d+\],?\s*)+)"

        counts = re.findall(count_pattern, content)
        sequences = re.findall(sequence_pattern, content)
        positions_blocks = re.findall(positions_pattern, content)

        all_positions = []
        for block in positions_blocks:
            pos_pairs = re.findall(r"\[(\d+),(\d+)\]", block)
            pos_pairs = [(int(start), int(end)) for start, end in pos_pairs]
            all_positions.append(pos_pairs)

        data = []
        for count, seq, pos in zip(counts, sequences, all_positions):
            entry = {
                'count': int(count),
                'sequence': seq,
                'positions': pos
            }
            data.append(entry)

        return data

    except Exception as e:
        print(f"解析文件时发生错误: {str(e)}")
        return []

def group_by_length(merged_intervals: List[NamedTuple], length_threshold: float = 0.2) -> Dict[int, List[NamedTuple]]:
    """
    将区间按照长度分组

    Args:
        merged_intervals: 合并后的区间列表
        length_threshold: 长度相似度阈值，默认0.2表示长度差异在20%以内认为相似

    Returns:
        按长度分组后的区间字典，key为组ID，value为区间列表
    """
    if not merged_intervals:
        return {}

    # 计算所有区间的长度
    interval_lengths = [(interval, interval.end - interval.start) for interval in merged_intervals]

    # 按长度排序
    interval_lengths.sort(key=lambda x: x[1])

    length_groups = defaultdict(list)
    current_group_id = 0
    current_group_base_length = interval_lengths[0][1]

    for interval, length in interval_lengths:
        # 计算与当前组基准长度的相对差异
        length_diff_ratio = abs(length - current_group_base_length) / current_group_base_length

        # 如果长度差异超过阈值，创建新组
        if length_diff_ratio > length_threshold:
            current_group_id += 1
            current_group_base_length = length

        length_groups[current_group_id].append(interval)

    return length_groups

def write_grouping_results(
        final_groups: Dict,
        sequence_data: List[Dict],
        output_file: str
):
    """写入分组结果到CSV文件"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 生成扁平化的组列表
        flat_groups = []
        group_counter = 0

        # 遍历所有组并扁平化
        for length_group_id in sorted(final_groups.keys()):
            for sim_group_id in sorted(final_groups[length_group_id].keys()):
                group = final_groups[length_group_id][sim_group_id]
                flat_groups.append((group_counter, group))
                group_counter += 1

        # 按组ID排序并写入
        for group_id, group in sorted(flat_groups):
            writer.writerow([f'Group {group_id}'])
            writer.writerow(['Start', 'End', 'Length'])

            # 写入该组的所有区间
            for merged in sorted(group['merged'], key=lambda x: x.start):
                writer.writerow([
                    merged.start,
                    merged.end,
                    merged.end - merged.start
                ])

            # 每个组后添加空行
            writer.writerow([])
            writer.writerow([])

def compare_sequence_pair(seq1,seq2):
    aligner = maxSimilarBoundary.initialize_aligner()

    score = aligner.score(seq1, seq2)
    max_possible_score = max(len(seq1),len(seq2)) * aligner.match_score
    similarity = round(score / max_possible_score if max_possible_score > 0 else 0, 3)

    return similarity

def group_sequences_by_similarity(sequences: List[Dict], similarity_threshold: float = 0.7) -> Dict[int, Set[int]]:
    """
    将序列按相似度分组，使用传递关系

    Args:
        sequences: 序列数据列表
        similarity_threshold: 相似度阈值

    Returns:
        Dict[int, Set[int]]: 分组结果，key为组ID，value为序列索引集合
    """
    n = len(sequences)
    groups = {}
    current_group_id = 0
    processed = set()

    for i in range(n):
        if i in processed:
            continue

        current_group = {i}
        seq1 = sequences[i]['sequence']

        # 考虑传递关系
        to_check = {i}
        while to_check:
            current_idx = to_check.pop()
            seq1 = sequences[current_idx]['sequence']

            for j in range(n):
                if j in processed or j in current_group:
                    continue

                seq2 = sequences[j]['sequence']
                similarity = compare_sequence_pair(seq1, seq2)

                if similarity >= similarity_threshold:
                    current_group.add(j)
                    to_check.add(j)

        processed.update(current_group)
        groups[current_group_id] = current_group
        current_group_id += 1

    return groups

def process_two_level_grouping(
        sequence_data: List[Dict],
        merged_intervals: List[NamedTuple],
        length_threshold: float = 0.2,
        similarity_threshold: float = 0.7,
        genome: Optional[str] = None
) -> Dict[int, Dict[int, Dict]]:
    """
    两级分组：先按长度分组，再按序列相似度分组

    Args:
        sequence_data: 序列数据列表
        merged_intervals: 合并后的区间列表
        length_threshold: 长度相似度阈值
        similarity_threshold: 序列相似度阈值

    Returns:
        Dict: 两级分组结果
    """
    genome_seq = load_genome_sequence(genome)
    # 创建日志文件
    with open('grouping_process.txt', 'w', encoding='utf-8') as f:
        f.write("=============== 序列分组过程 ===============\n")

    def log_debug(message):
        with open('grouping_process.txt', 'a', encoding='utf-8') as f:
            f.write(message + '\n')
        print(message)

    # 第一级：按长度分组
    log_debug("开始按长度分组...")
    length_groups = group_by_length(merged_intervals, length_threshold)
    log_debug(f"形成了 {len(length_groups)} 个长度组")

    final_groups = {}

    # 在每个长度组内进行序列相似度分组
    for length_group_id, intervals in length_groups.items():
        log_debug(f"\n处理长度组 {length_group_id}")
        final_groups[length_group_id] = {}

        # 找出与当前长度组相关的序列
        relevant_sequences = []
        relevant_sequence_indices = []

        # 找出所有与当前长度组区间重叠的序列
        for seq_idx, seq in enumerate(sequence_data):
            for interval in intervals:
                for pos in seq['positions']:
                    if not (interval.end <= pos[0] or pos[1] <= interval.start):
                        if seq_idx not in relevant_sequence_indices:
                            # 从基因组中提取实际序列
                            if genome_seq:
                                actual_sequence = extract_sequence(genome_seq, pos[0], pos[1])
                                seq_copy = seq.copy()
                                seq_copy['sequence'] = actual_sequence
                                relevant_sequences.append(seq_copy)
                            else:
                                relevant_sequences.append(seq)
                            relevant_sequence_indices.append(seq_idx)
                        break
                if seq_idx in relevant_sequence_indices:
                    break

        log_debug(f"- 找到 {len(relevant_sequences)} 个相关序列")

        if not relevant_sequences:
            continue

        # 对相关序列进行相似度分组
        similarity_groups = group_sequences_by_similarity(
            relevant_sequences,
            similarity_threshold
        )

        log_debug(f"- 形成 {len(similarity_groups)} 个相似度组")

        # 保存分组结果
        for sim_group_id, seq_indices in similarity_groups.items():
            # 转换回原始序列索引
            original_indices = {relevant_sequence_indices[i] for i in seq_indices}

            final_groups[length_group_id][sim_group_id] = {
                'sequences': original_indices,
                'merged': intervals
            }

            log_debug(f"  - 相似度组 {sim_group_id}: {len(original_indices)} 个序列")
            # 记录实际序列内容
            if genome_seq:
                log_debug("  实际序列内容:")
                for i in seq_indices:
                    seq = relevant_sequences[i]
                    log_debug(f"    序列 {relevant_sequence_indices[i]}: {seq['sequence']}")
    return final_groups


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='区间操作和序列分析工具')

    parser.add_argument('--interval-file',
                        type=str,
                        required=True,
                        help='包含区间数据的CSV文件路径')

    parser.add_argument('--sequence-file',
                        type=str,
                        required=True,
                        help='包含重复序列数据的TXT文件路径')

    parser.add_argument('--output-file',
                        type=str,
                        default='sequential_grouping_results.csv',
                        help='输出结果的CSV文件路径 (默认: sequential_grouping_results.csv)')

    parser.add_argument('--genome', '-g', required=True,type=str,
                        help='Input FASTA file path containing genome sequence')
    # 添加长度阈值参数
    parser.add_argument(
        '--length-threshold',
        type=float,
        default=0.2,
        help='长度相似度阈值(0-1之间)，默认0.2表示长度差异在20%以内认为相似'
    )
    parser.add_argument('--start-col1',
                        type=int,
                        default=1,
                        help='第一个起始列索引 (默认: 1)')

    parser.add_argument('--end-col1',
                        type=int,
                        default=2,
                        help='第一个结束列索引 (默认: 2)')

    parser.add_argument('--start-col2',
                        type=int,
                        default=3,
                        help='第二个起始列索引 (默认: 3)')

    parser.add_argument('--end-col2',
                        type=int,
                        default=4,
                        help='第二个结束列索引 (默认: 4)')

    return parser.parse_args()


def main():
    args = parse_arguments()

    # 初始化区间操作类
    ops = IntervalOperations()
    print("读取merged区间...")

    # 读取区间
    intervals = ops.read_intervals_from_csv(
        args.interval_file,
        args.start_col1,
        args.end_col1,
        args.start_col2,
        args.end_col2
    )

    # 计算并集
    merged = ops.merge_intervals(intervals)
    print(f"读取到 {len(merged)} 个merged区间")

    # 读取序列数据
    print("\n读取核心序列...")
    sequence_data = parse_sequences(args.sequence_file)
    print(f"读取到 {len(sequence_data)} 个核心序列")

    # 使用新的两级分组处理
    final_groups = process_two_level_grouping(
        sequence_data,
        merged,
        length_threshold=args.length_threshold,
        similarity_threshold=0.7,
        genome=args.genome
    )
    # 保存结果
    write_grouping_results(final_groups, sequence_data, args.output_file)
    print(f"\n结果已保存到: {args.output_file}")


if __name__ == "__main__":
    main()