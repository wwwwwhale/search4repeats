import csv
from typing import List, Optional
import re
from collections import defaultdict, namedtuple
from typing import List, Dict, Tuple, Set
import numpy as np
from sklearn.cluster import KMeans
import argparse


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

                # 重叠度超过40%才合并
                if overlap_ratio > 0.4:
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


def process_merged_intervals_sequentially(sequence_data: List[Dict],
                                          merged_intervals: List[namedtuple],
                                          interval_ops) -> Dict[int, List[Dict]]:
    """逐个处理merged区间，构建分组"""
    groups = defaultdict(lambda: {'sequences': set(), 'merged': []})
    current_group_id = 0

    # 创建或清空日志文件
    with open('sequence_grouping_process.txt', 'w', encoding='utf-8') as f:
        f.write("=============== 序列分组详细过程 ===============\n")

    def log_debug(message):
        """写入调试信息到文件"""
        with open('sequence_grouping_process.txt', 'a', encoding='utf-8') as f:
            f.write(message + '\n')
        print(message)  # 保持控制台输出

    def find_overlapping_sequences(merged_interval, sequences):
        overlaps = []
        log_debug(f"\n检查重叠情况:")
        log_debug(f"  当前merged区间: [{merged_interval.start}, {merged_interval.end}]")

        for seq_idx, seq in enumerate(sequences):
            log_debug(f"\n  检查序列 {seq_idx}:")
            for pos_idx, pos in enumerate(seq['positions']):
                core_interval = interval_ops.Interval(pos[0], pos[1])

                # 计算重叠度
                overlap_start = max(merged_interval.start, core_interval.start)
                overlap_end = min(merged_interval.end, core_interval.end)

                if overlap_start <= overlap_end:  # 有重叠
                    # 计算重叠长度
                    overlap_length = overlap_end - overlap_start
                    # 计算两个区间长度
                    merged_length = merged_interval.end - merged_interval.start
                    core_length = core_interval.end - core_interval.start
                    # 使用较短区间作为基准计算重叠度
                    overlap_ratio = overlap_length / min(merged_length, core_length)

                    log_debug(f"    位置 {pos_idx}: [{core_interval.start}, {core_interval.end}]")
                    log_debug(f"      重叠区间: [{overlap_start}, {overlap_end}]")
                    log_debug(f"      重叠长度: {overlap_length}")
                    log_debug(f"      重叠度: {overlap_ratio:.3f}")

                    if overlap_ratio > 0.5:
                        overlaps.append((seq_idx, core_interval))
                        log_debug(f"      → 重叠度超过0.5，添加到组!")
                    else:
                        log_debug(f"      → 重叠度不足0.5，忽略")
                else:
                    log_debug(f"    位置 {pos_idx}: [{core_interval.start}, {core_interval.end}] → 不重叠")

        return overlaps

    log_debug("\n============= 开始处理merged区间 =============")
    log_debug(f"初始数据:\n- 序列总数: {len(sequence_data)}\n- 待处理merged区间数: {len(merged_intervals)}")

    for merged_idx, merged in enumerate(merged_intervals):
        log_debug(f"\n【处理第 {merged_idx + 1}/{len(merged_intervals)} 个merged区间】")
        log_debug(f"区间位置: [{merged.start}, {merged.end}]")

        overlaps = find_overlapping_sequences(merged, sequence_data)

        if not overlaps:
            log_debug(f"\n❌ 没有找到重叠的核心序列，跳过此区间")
            continue

        log_debug(f"\n✓ 找到 {len(overlaps)} 个重叠的核心序列:")
        for seq_idx, core_interval in overlaps:
            log_debug(f"  - 序列 {seq_idx}: [{core_interval.start}, {core_interval.end}]")

        found_groups = set()
        log_debug("\n检查序列所属组:")
        for seq_idx, _ in overlaps:
            for group_id, group in groups.items():
                if seq_idx in group['sequences']:
                    found_groups.add(group_id)
                    log_debug(f"  序列 {seq_idx} 属于组 {group_id}")

        if not found_groups:
            group_id = current_group_id
            current_group_id += 1
            log_debug(f"\n➡️ 创建新组 {group_id}")
        else:
            group_id = min(found_groups)
            log_debug(f"\n➡️ 使用现有组 {group_id}")

            if len(found_groups) > 1:
                log_debug(f"\n⚠️ 需要合并组: {sorted(found_groups)}")
                main_group_id = min(found_groups)

                for other_group_id in sorted(found_groups - {main_group_id}):
                    log_debug(f"\n合并组 {other_group_id} 到组 {main_group_id}:")
                    log_debug(f"- 组 {other_group_id} 原有序列: {sorted(groups[other_group_id]['sequences'])}")
                    log_debug(f"- 组 {other_group_id} merged区间数: {len(groups[other_group_id]['merged'])}")

                    groups[main_group_id]['sequences'].update(groups[other_group_id]['sequences'])
                    groups[main_group_id]['merged'].extend(groups[other_group_id]['merged'])
                    del groups[other_group_id]

                    log_debug(f"合并后组 {main_group_id} 状态:")
                    log_debug(f"- 序列: {sorted(groups[main_group_id]['sequences'])}")
                    log_debug(f"- merged区间数: {len(groups[main_group_id]['merged'])}")

        groups[group_id]['sequences'].update(seq_idx for seq_idx, _ in overlaps)
        groups[group_id]['merged'].append(merged)

        log_debug(f"\n组 {group_id} 最终状态:")
        log_debug(f"- 序列: {sorted(groups[group_id]['sequences'])}")
        log_debug(f"- merged区间数: {len(groups[group_id]['merged'])}")

    log_debug("\n============= 分组完成 =============")
    log_debug(f"最终形成 {len(groups)} 个组:")

    for group_id in sorted(groups.keys()):
        group = groups[group_id]
        log_debug(f"\n组 {group_id}:")
        log_debug(f"- 序列 ({len(group['sequences'])}个): {sorted(group['sequences'])}")
        log_debug(f"- merged区间 ({len(group['merged'])}个):")
        for idx, interval in enumerate(group['merged']):
            log_debug(f"  {idx + 1}. [{interval.start}, {interval.end}]")

    return groups


def write_sequential_results(groups: Dict[int, Dict],
                             sequence_data: List[Dict],
                             output_file: str):
    """将顺序处理的结果写入CSV文件"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for group_id, group in sorted(groups.items()):
            writer.writerow([f'Group {group_id + 1}'])
            writer.writerow(['Start', 'End', 'Length'])
            for merged in group['merged']:
                writer.writerow([
                    merged.start,
                    merged.end,
                    merged.end - merged.start
                ])
            writer.writerow([])
            writer.writerow([])


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
    # 解析命令行参数
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

    # 顺序处理merged区间
    groups = process_merged_intervals_sequentially(
        sequence_data,
        merged,
        ops
    )

    # 保存结果
    write_sequential_results(groups, sequence_data, args.output_file)
    print(f"\n结果已保存到: {args.output_file}")


if __name__ == "__main__":
    main()