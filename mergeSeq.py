import csv
from typing import List, Optional
import re
from collections import defaultdict, namedtuple
import argparse
from collections import defaultdict
from typing import List, Dict, NamedTuple, Tuple, Set
import numpy as np

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

                    if overlap_ratio == 1:
                        overlaps.append((seq_idx, core_interval))
                        log_debug(f"      → 重叠度超过1，添加到组!")
                    else:
                        log_debug(f"      → 重叠度不足1，忽略")
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


def process_merged_intervals_with_length_groups(
        sequence_data: List[Dict],
        merged_intervals: List[NamedTuple],
        interval_ops,
        length_threshold: float = 0.2
) -> Dict[int, Dict]:
    """
    两级分组处理：先按长度分组，再在组内按重叠关系分组

    Args:
        sequence_data: 序列数据
        merged_intervals: 合并后的区间列表
        interval_ops: 区间操作工具类
        length_threshold: 长度相似度阈值

    Returns:
        两级分组结果，格式为：
        {
            length_group_id: {
                overlap_group_id: {
                    'sequences': set(),
                    'merged': []
                }
            }
        }
    """
    # 创建日志文件
    with open('sequence_grouping_process.txt', 'w', encoding='utf-8') as f:
        f.write("=============== 序列分组详细过程 ===============\n")

    def log_debug(message):
        with open('sequence_grouping_process.txt', 'a', encoding='utf-8') as f:
            f.write(message + '\n')
        print(message)

    # 第一级：按长度分组
    length_groups = group_by_length(merged_intervals, length_threshold)
    log_debug(f"\n共形成 {len(length_groups)} 个长度组")

    # 记录每个长度组的信息
    for group_id, intervals in length_groups.items():
        lengths = [interval.end - interval.start for interval in intervals]
        log_debug(f"\n长度组 {group_id}:")
        log_debug(f"- 区间数量: {len(intervals)}")
        log_debug(f"- 长度范围: {min(lengths):.2f} - {max(lengths):.2f}")
        log_debug(f"- 平均长度: {np.mean(lengths):.2f}")

    # 第二级：在每个长度组内进行重叠分组
    final_groups = {}

    for length_group_id, length_group_intervals in length_groups.items():
        log_debug(f"\n处理长度组 {length_group_id}")

        # 使用原有的分组逻辑处理当前长度组的区间
        overlap_groups = process_merged_intervals_sequentially(
            sequence_data,
            length_group_intervals,
            interval_ops
        )

        final_groups[length_group_id] = overlap_groups

        log_debug(f"\n长度组 {length_group_id} 处理完成，形成 {len(overlap_groups)} 个重叠组")

    # 输出最终分组统计
    log_debug("\n=============== 最终分组统计 ===============")
    total_overlap_groups = sum(len(overlap_groups) for overlap_groups in final_groups.values())
    log_debug(f"长度组数量: {len(final_groups)}")
    log_debug(f"总重叠组数量: {total_overlap_groups}")

    for length_group_id, overlap_groups in final_groups.items():
        log_debug(f"\n长度组 {length_group_id}:")
        log_debug(f"- 包含 {len(overlap_groups)} 个重叠组")
        for overlap_group_id, group in overlap_groups.items():
            log_debug(f"  - 重叠组 {overlap_group_id}: {len(group['sequences'])} 个序列, "
                      f"{len(group['merged'])} 个merged区间")

    return final_groups


def write_sequential_results(final_groups: Dict[int, Dict],
                             sequence_data: List[Dict],
                             output_file: str):
    """将分组结果写入CSV文件，保持原始格式"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 生成扁平化的组列表
        flat_groups = []
        group_counter = 0

        # 遍历所有组并扁平化
        for length_groups in final_groups.values():
            for overlap_group in length_groups.values():
                flat_groups.append((group_counter, overlap_group))
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


def write_summary_file(final_groups: Dict[int, Dict],
                       output_file: str):
    """写入分组摘要信息到单独的文件

    Args:
        final_groups: 两级分组结果
        output_file: 输出文件路径
    """
    summary_file = output_file.replace('.csv', '_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 写入总体统计
        writer.writerow(['Groups Summary'])
        writer.writerow(['Number of Length Groups', len(final_groups)])
        total_overlap_groups = sum(len(groups) for groups in final_groups.values())
        writer.writerow(['Total Overlap Groups', total_overlap_groups])
        writer.writerow([])

        # 写入每个长度组的统计
        writer.writerow(['Length Group', 'Overlap Groups', 'Total Intervals', 'Avg Length'])
        for length_group_id in sorted(final_groups.keys()):
            overlap_groups = final_groups[length_group_id]

            # 统计信息
            total_intervals = sum(len(group['merged'])
                                  for group in overlap_groups.values())
            all_lengths = []
            for group in overlap_groups.values():
                all_lengths.extend(interval.end - interval.start
                                   for interval in group['merged'])
            avg_length = sum(all_lengths) / len(all_lengths) if all_lengths else 0

            writer.writerow([
                length_group_id,
                len(overlap_groups),
                total_intervals,
                f'{avg_length:.2f}'
            ])


# 使用示例：
def save_results(final_groups: Dict[int, Dict],
                 sequence_data: List[Dict],
                 output_file: str):
    """保存所有结果文件

    Args:
        final_groups: 两级分组结果
        sequence_data: 序列数据
        output_file: 基础输出文件路径
    """
    # 写入详细结果
    write_sequential_results(final_groups, sequence_data, output_file)

    # 写入摘要信息
    write_summary_file(final_groups, output_file)

    print(f"\n详细结果已保存到: {output_file}")
    print(f"摘要信息已保存到: {output_file.replace('.csv', '_summary.csv')}")


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
    # 使用命令行参数中的阈值
    final_groups = process_merged_intervals_with_length_groups(
        sequence_data,
        merged,
        ops,
        length_threshold=args.length_threshold
    )
    # 保存结果
    write_sequential_results(final_groups, sequence_data, args.output_file)
    print(f"\n结果已保存到: {args.output_file}")


if __name__ == "__main__":
    main()