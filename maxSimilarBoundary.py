import re
import csv
import multiprocessing as mp
import time
from collections import defaultdict
from itertools import combinations
import argparse

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio import SeqIO

# Global variables at module top
STOPPED_EXTENSIONS = []

def parse_sequences(file_path):
    """
    解析TXT文件，提取重复序列、出现次数和位置。

    参数:
    - file_path (str): TXT文件路径

    返回:
    - list: 包含序列信息的字典列表
    """
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

def convert_to_2d_array(sequence_data):
    """
    将序列数据转换为窗口对的形式

    返回:
    - positions_2d: 二维数组，每个元素是该序列的窗口对列表
    - sequence_info: 对应的序列信息字典
    """
    positions_2d = []
    sequence_info = {}

    for idx, data in enumerate(sequence_data):
        positions = data['positions']
        window_pairs = list(combinations(positions, 2))
        positions_2d.append(window_pairs)
        sequence_info[idx] = (data['sequence'], data['count'])

    return positions_2d, sequence_info

def load_genome_sequence(gbff_file):
    """
    加载基因组序列。

    参数:
    - gbff_file (str): GBFF文件路径

    返回:
    - str: 基因组序列
    """
    try:
        for record in SeqIO.parse(gbff_file, "genbank"):
            genome_seq = str(record.seq)
            print(f"基因组序列 '{record.id}' 加载完成，长度为 {len(genome_seq)} 个碱基。")
            return genome_seq

        print("没有找到基因组序列。")
        return ""
    except Exception as e:
        print(f"加载基因组序列时发生错误: {e}")
        return ""

def initialize_aligner():
    """
    初始化PairwiseAligner对象。

    返回:
    - PairwiseAligner: 配置好的对齐器
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    return aligner

def compare_sequence_pair(args):
    """
    比对一个窗口对中的两个序列并计算相似度
    """
    window_pair, genome_seq, seq_idx, pair_idx = args
    aligner = initialize_aligner()

    try:
        start1, end1 = window_pair[0]
        start2, end2 = window_pair[1]

        seq1 = genome_seq[start1:end1]
        seq2 = genome_seq[start2:end2]

        score = aligner.score(seq1, seq2)
        max_possible_score = len(seq1) * aligner.match_score
        similarity = round(score / max_possible_score if max_possible_score > 0 else 0, 3)

        return (seq_idx, similarity, (start1, end1), (start2, end2), pair_idx)
    except Exception as e:
        print(f"比对序列时出错 (序列{seq_idx}): {e}")
        return (seq_idx, 0, None, None, pair_idx)

def analyze_sequence_similarities(positions_2d, genome_seq, pool, threshold=0.8):
    """
    使用外部传入的进程池进行并行计算
    """
    comparison_tasks = []
    sequences_need_extension = set()
    input_pair_indices = defaultdict(set)

    for seq_idx, window_pairs in enumerate(positions_2d):
        for pair_idx, window_pair in enumerate(window_pairs):
            comparison_tasks.append((window_pair, genome_seq, seq_idx, pair_idx))
            input_pair_indices[seq_idx].add(pair_idx)

    # 使用传入的进程池
    results = pool.map(compare_sequence_pair, comparison_tasks)

    high_similarity_pairs = defaultdict(list)
    low_similarity_pairs = defaultdict(list)
    processed_pair_indices = defaultdict(set)

    for result in results:
        seq_idx, similarity, pos1, pos2, pair_idx = result
        if pos1 and pos2:
            processed_pair_indices[seq_idx].add(pair_idx)

            if similarity > threshold:
                high_similarity_pairs[seq_idx].append((pos1, pos2, similarity, pair_idx))
                sequences_need_extension.add(seq_idx)
            else:
                low_similarity_pairs[seq_idx].append({
                    'positions': (pos1, pos2),
                    'similarity': similarity,
                    'stop_reason': 'similarity_below_threshold',
                    'pair_idx': pair_idx
                })

    return high_similarity_pairs, low_similarity_pairs, sequences_need_extension

def compare_direction_similarity(genome_seq, pos1, pos2, extend_length, direction='left'):
    """
    比较特定方向扩展的序列相似度
    """
    aligner = initialize_aligner()

    if direction == 'left':
        start1 = max(0, pos1[0] - extend_length)
        start2 = max(0, pos2[0] - extend_length)
        seq1 = genome_seq[start1:pos1[0]]
        seq2 = genome_seq[start2:pos2[0]]
    else:
        end1 = min(len(genome_seq), pos1[1] + extend_length)
        end2 = min(len(genome_seq), pos2[1] + extend_length)
        seq1 = genome_seq[pos1[1]:end1]
        seq2 = genome_seq[pos2[1]:end2]

    if not seq1 or not seq2:
        return 0.0

    score = aligner.score(seq1, seq2)
    max_possible_score = len(seq1) * aligner.match_score
    return score / max_possible_score if max_possible_score > 0 else 0

def extend_sequence_positions(window_pairs, extend_length, genome_length, genome_seq):
    """
    扩展窗口对
    """
    extended = []
    for pos1, pos2, _, pair_idx in window_pairs:
        left_similarity = compare_direction_similarity(genome_seq, pos1, pos2, extend_length, 'left')
        right_similarity = compare_direction_similarity(genome_seq, pos1, pos2, extend_length, 'right')

        new_start1, new_end1 = pos1[0], pos1[1]
        new_start2, new_end2 = pos2[0], pos2[1]

        if left_similarity > right_similarity:
            new_start1 = max(0, pos1[0] - extend_length)
            new_start2 = max(0, pos2[0] - extend_length)
        else:
            new_end1 = min(genome_length, pos1[1] + extend_length)
            new_end2 = min(genome_length, pos2[1] + extend_length)

        extended.append(((new_start1, new_end1), (new_start2, new_end2), pair_idx))

    return extended

def iterative_sequence_extension(positions_2d, genome_seq, extend_length, threshold=0.8, max_iterations=10, origin_array=None):
    """
    在外层创建进程池
    """
    genome_length = len(genome_seq)
    iteration = 0
    all_stopped_pairs = defaultdict(list)

    # 获取CPU核心数并创建进程池
    cpu_count = mp.cpu_count()
    print(f"系统CPU核心数: {cpu_count}")

    # 创建一个进程池在整个迭代过程中重用
    with mp.Pool() as pool:
        while iteration < max_iterations:
            iteration_start_time = time.time()
            
            high_similarity_pairs, low_similarity_pairs, sequences_need_extension = analyze_sequence_similarities(
                positions_2d,
                genome_seq,
                pool,  # 传入进程池
                threshold
            )

            # 处理低相似度对
            for seq_idx, pairs_info in low_similarity_pairs.items():
                for pair_info in pairs_info:
                    pos1, pos2 = pair_info['positions']
                    pair_idx = pair_info['pair_idx']
                    original_pair = origin_array[seq_idx][pair_idx]

                    all_stopped_pairs[seq_idx].append({
                        'current_start1': pos1[0],
                        'current_end1': pos1[1],
                        'current_start2': pos2[0],
                        'current_end2': pos2[1],
                        'length': pos1[1] - pos1[0] + 1
                    })

            # 检查是否还有需要扩展的序列
            if not high_similarity_pairs:
                print(f"没有满足相似度阈值的窗口对，迭代结束。共进行了 {iteration + 1} 次扩展")
                break

            # 扩展序列
            extension_count = 0
            for seq_idx, window_pairs in high_similarity_pairs.items():
                if seq_idx < len(positions_2d):
                    positions_2d[seq_idx] = extend_sequence_positions(
                        window_pairs,
                        extend_length,
                        genome_length,
                        genome_seq
                    )
                    extension_count += 1

            iteration_time = time.time() - iteration_start_time
            print(f"第 {iteration + 1} 次迭代完成，扩展了 {extension_count} 个序列")
            print(f"耗时: {iteration_time:.2f} 秒")

            iteration += 1

        if iteration >= max_iterations:
            print(f"达到最大迭代次数 {max_iterations}，强制停止扩展")

    return all_stopped_pairs

def save_results_to_csv(stopped_pairs, output_file="result.csv"):
    """
    保存结果到CSV
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'seq_idx',
            'current_start1', 'current_end1', 'current_start2', 'current_end2',
            'length'
        ])

        for seq_idx, pairs in stopped_pairs.items():
            for pair in pairs:
                writer.writerow([
                    seq_idx,
                    pair['current_start1'], pair['current_end1'],
                    pair['current_start2'], pair['current_end2'],
                    pair['length']
                ])

def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Sequence Analysis Tool')

    parser.add_argument('--input', '-i', required=True,
                        help='Input TXT file path containing repeated subsequences')
    parser.add_argument('--genome', '-g', required=True,
                        help='Input GBFF file path containing genome sequence')
    parser.add_argument('--output', '-o', required=True,
                        help='Output CSV file path for results')
    parser.add_argument('--extend-length', type=int, default=100,
                        help='Length to extend sequences (default: 100)')
    parser.add_argument('--threshold', type=float, default=0.9,
                        help='Similarity threshold (default: 0.9)')
    parser.add_argument('--max-iterations', type=int, default=100000,
                        help='Maximum number of iterations (default: 100000)')

    return parser.parse_args()

def main():
    args = parse_arguments()

    sequence_data = parse_sequences(args.input)
    if not sequence_data:
        print("No sequence data parsed, terminating program.")
        return

    positions_array, seq_info = convert_to_2d_array(sequence_data)
    origin_array = positions_array.copy()

    genome_seq = load_genome_sequence(args.genome)
    if not genome_seq:
        print("Genome sequence loading failed, terminating program.")
        return

    stopped_pairs = iterative_sequence_extension(
        positions_array,
        genome_seq,
        extend_length=args.extend_length,
        threshold=args.threshold,
        max_iterations=args.max_iterations,
        origin_array=origin_array
    )

    save_results_to_csv(stopped_pairs, args.output)
    print("Program execution completed.")

if __name__ == '__main__':
    main()