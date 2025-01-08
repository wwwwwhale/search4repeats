import re
import logging
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

def initialize_logger(log_file="sequence_analysis.log", console_level=logging.INFO, file_level=logging.DEBUG):
    """Initialize logger with specified parameters."""
    # Convert string level to logging level if necessary
    if isinstance(console_level, str):
        console_level = getattr(logging, console_level.upper())
    if isinstance(file_level, str):
        file_level = getattr(logging, file_level.upper())

    # Clear log file
    with open(log_file, 'w', encoding='utf-8') as f:
        f.write('')

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    if not logger.handlers:
        # File handler
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(file_level)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(console_level)
        console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

    return logger

def parse_sequences(file_path, logger):
    """
    解析TXT文件，提取重复序列、出现次数和位置。

    参数:
    - file_path (str): TXT文件路径
    - logger (Logger): 日志记录器

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

            # 转换为整数对，注意不要进行字符串拼接
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
        # 生成所有可能的窗口对组合
        window_pairs = list(combinations(positions, 2))
        positions_2d.append(window_pairs)
        sequence_info[idx] = (data['sequence'], data['count'])

    return positions_2d, sequence_info


def load_genome_sequence(fasta_file, logger):
    """
    加载基因组序列。

    参数:
    - fasta_file (str): FASTA文件路径
    - logger (Logger): 日志记录器

    返回:
    - str: 基因组序列
    """
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            genome_seq = str(record.seq)
            logger.info(f"基因组序列 '{record.id}' 加载完成，长度为 {len(genome_seq)} 个碱基。")
            return genome_seq
        logger.error("没有找到基因组序列。")
    except Exception as e:
        logger.error(f"加载基因组序列时发生错误: {e}")
    return ""


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


def compare_sequence_pair(args):
    """
    比对一个窗口对中的两个序列并计算相似度，现在包含位置索引
    """
    window_pair, genome_seq, seq_idx, pair_idx = args  # 添加pair_idx参数
    aligner = initialize_aligner()

    try:
        start1, end1 = window_pair[0]
        start2, end2 = window_pair[1]

        seq1 = genome_seq[start1:end1]
        seq2 = genome_seq[start2:end2]

        score = aligner.score(seq1, seq2)
        max_possible_score = len(seq1) * aligner.match_score
        similarity = round(score / max_possible_score if max_possible_score > 0 else 0, 3)

        # 返回结果中包含pair_idx
        return (seq_idx, similarity, (start1, end1), (start2, end2), pair_idx)
    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.error(f"比对序列时出错 (序列{seq_idx}): {e}")
        return (seq_idx, 0, None, None, pair_idx)


def analyze_sequence_similarities(positions_2d, genome_seq, threshold=0.8):
    """
    分析序列相似度，并验证窗口对完整性
    """
    logger = logging.getLogger(__name__)
    comparison_tasks = []
    sequences_need_extension = set()

    # 记录输入的窗口对索引
    input_pair_indices = defaultdict(set)  # {seq_idx: set(pair_idx)}
    for seq_idx, window_pairs in enumerate(positions_2d):
        for pair_idx, window_pair in enumerate(window_pairs):
            comparison_tasks.append((
                window_pair,
                genome_seq,
                seq_idx,
                pair_idx
            ))
            input_pair_indices[seq_idx].add(pair_idx)

    # 并行执行比对
    with mp.Pool() as pool:
        results = pool.map(compare_sequence_pair, comparison_tasks)

    # 处理结果
    high_similarity_pairs = defaultdict(list)
    low_similarity_pairs = defaultdict(list)
    processed_pair_indices = defaultdict(set)  # {seq_idx: set(pair_idx)}

    for result in results:
        seq_idx, similarity, pos1, pos2, pair_idx = result
        if pos1 and pos2:  # 确保位置有效
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

    # 验证每个序列的窗口对完整性
    for seq_idx in input_pair_indices:
        input_pairs = input_pair_indices[seq_idx]
        processed_pairs = processed_pair_indices[seq_idx]

        # 检查是否有丢失的窗口对
        missing_pairs = input_pairs - processed_pairs
        if missing_pairs:
            logger.warning(f"序列{seq_idx}丢失了以下窗口对: {missing_pairs}")

        # 检查是否有多出的窗口对
        extra_pairs = processed_pairs - input_pairs
        if extra_pairs:
            logger.warning(f"序列{seq_idx}多出了以下窗口对: {extra_pairs}")

        # 统计该序列的窗口对分布
        high_sim_pairs = {p[3] for p in high_similarity_pairs[seq_idx]}  # 获取pair_idx
        low_sim_pairs = {p['pair_idx'] for p in low_similarity_pairs[seq_idx]}

        # 检查是否所有处理过的窗口对都被正确分类
        unclassified_pairs = processed_pairs - (high_sim_pairs | low_sim_pairs)
        if unclassified_pairs:
            logger.warning(f"序列{seq_idx}的以下窗口对未被分类: {unclassified_pairs}")

        # 输出详细统计
        logger.info(f"序列{seq_idx}统计:")
        logger.info(f"  - 输入窗口对数量: {len(input_pairs)}")
        logger.info(f"  - 处理窗口对数量: {len(processed_pairs)}")
        logger.info(f"  - 高相似度窗口对: {len(high_sim_pairs)}")
        logger.info(f"  - 低相似度窗口对: {len(low_sim_pairs)}")

    return high_similarity_pairs, low_similarity_pairs, sequences_need_extension


def compare_direction_similarity(genome_seq, pos1, pos2, extend_length, direction='left', logger=None):
    """
    比较特定方向扩展的序列相似度，并记录详细日志
    """
    aligner = initialize_aligner()

    if direction == 'left':
        start1 = max(0, pos1[0] - extend_length)
        start2 = max(0, pos2[0] - extend_length)
        seq1 = genome_seq[start1:pos1[0]]
        seq2 = genome_seq[start2:pos2[0]]
        if logger:
            logger.debug(f"左向扩展: 位置1 [{start1}-{pos1[0]}], 位置2 [{start2}-{pos2[0]}]")
    else:  # right
        end1 = min(len(genome_seq), pos1[1] + extend_length)
        end2 = min(len(genome_seq), pos2[1] + extend_length)
        seq1 = genome_seq[pos1[1]:end1]
        seq2 = genome_seq[pos2[1]:end2]
        if logger:
            logger.debug(f"右向扩展: 位置1 [{pos1[1]}-{end1}], 位置2 [{pos2[1]}-{end2}]")

    if not seq1 or not seq2:
        if logger:
            logger.warning(f"序列为空: seq1长度={len(seq1)}, seq2长度={len(seq2)}")
        return 0.0

    score = aligner.score(seq1, seq2)
    max_possible_score = len(seq1) * aligner.match_score
    similarity = score / max_possible_score if max_possible_score > 0 else 0

    if logger:
        logger.debug(f"{direction}向扩展相似度: {similarity:.3f}")

    return similarity


def extend_sequence_positions(window_pairs, extend_length, genome_length, genome_seq):
    """
    扩展窗口对，保持位置索引的跟踪
    """
    extended = []
    for pos1, pos2, _, pair_idx in window_pairs:  # 现在包含pair_idx
        # 计算这对窗口的左右扩展相似度
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

        # 在扩展结果中保持pair_idx
        extended.append(((new_start1, new_end1), (new_start2, new_end2), pair_idx))

    return extended


def iterative_sequence_extension(positions_2d, genome_seq, extend_length, threshold=0.8, max_iterations=10, origin_array=None):
    """
    迭代扩展序列，跟踪窗口对的原始位置
    """
    logger = logging.getLogger(__name__)
    genome_length = len(genome_seq)
    iteration = 0
    all_stopped_pairs = defaultdict(list)

    while iteration < max_iterations:
        iteration_start_time = time.time()
        logger.info(f"\n===== 开始第 {iteration + 1} 次迭代 =====")

        # 分析序列相似度
        high_similarity_pairs, low_similarity_pairs, sequences_need_extension = analyze_sequence_similarities(
            positions_2d,
            genome_seq,
            threshold
        )

        logger.info(f"高相似度窗口对数量: {sum(len(pairs) for pairs in high_similarity_pairs.values())}")
        logger.info(f"低相似度窗口对数量: {sum(len(pairs) for pairs in low_similarity_pairs.values())}")

        # 记录本次迭代中停止扩展的窗口对
        for seq_idx, pairs_info in low_similarity_pairs.items():
            for pair_info in pairs_info:
                pos1, pos2 = pair_info['positions']
                pair_idx = pair_info['pair_idx']  # 获取原始位置索引
                original_pair = origin_array[seq_idx][pair_idx]  # 直接使用pair_idx获取原始对

                all_stopped_pairs[seq_idx].append({
                    'current_start1': pos1[0],
                    'current_end1': pos1[1],
                    'current_start2': pos2[0],
                    'current_end2': pos2[1],
                    'length': pos1[1] - pos1[0] + 1
                })

        # 如果没有高相似度的窗口对，结束迭代
        if not high_similarity_pairs:
            logger.info(f"没有满足相似度阈值的窗口对，迭代结束。共进行了 {iteration + 1} 次扩展")
            break

        # 扩展高相似度的窗口对
        extension_count = 0
        for seq_idx, window_pairs in high_similarity_pairs.items():
            if seq_idx < len(positions_2d):
                logger.debug(f"开始扩展序列 {seq_idx} (包含 {len(window_pairs)} 个窗口对)")
                positions_2d[seq_idx] = extend_sequence_positions(
                    window_pairs,
                    extend_length,
                    genome_length,
                    genome_seq
                )
                extension_count += 1

        iteration_time = time.time() - iteration_start_time
        logger.info(f"第 {iteration + 1} 次迭代完成")
        logger.info(f"本次迭代扩展了 {extension_count} 个序列")
        logger.info(f"耗时: {iteration_time:.2f} 秒")

        iteration += 1

    if iteration >= max_iterations:
        logger.warning(f"达到最大迭代次数 {max_iterations}，强制停止扩展")
        # ... [处理达到最大迭代次数的逻辑保持不变]

    return all_stopped_pairs




def save_results_to_csv(stopped_pairs, output_file="result.csv"):
    """
    保存结果到CSV，现在包含窗口对的位置索引
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
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Sequence Analysis Tool')

    # Required arguments
    parser.add_argument('--input', '-i', required=True,
                        help='Input TXT file path containing repeated subsequences')
    parser.add_argument('--genome', '-g', required=True,
                        help='Input FASTA file path containing genome sequence')
    parser.add_argument('--output', '-o', required=True,
                        help='Output CSV file path for results')

    # Optional arguments
    parser.add_argument('--log-file', default="sequence_analysis.log",
                        help='Log file path (default: sequence_analysis.log)')
    parser.add_argument('--extend-length', type=int, default=100,
                        help='Length to extend sequences (default: 100)')
    parser.add_argument('--threshold', type=float, default=0.9,
                        help='Similarity threshold (default: 0.9)')
    parser.add_argument('--max-iterations', type=int, default=100000,
                        help='Maximum number of iterations (default: 100000)')
    parser.add_argument('--console-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Console logging level (default: INFO)')
    parser.add_argument('--file-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='DEBUG', help='File logging level (default: DEBUG)')

    return parser.parse_args()


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Initialize logger with command line parameters
    logger = initialize_logger(
        log_file=args.log_file,
        console_level=args.console_level,
        file_level=args.file_level
    )

    # Parse sequence file
    sequence_data = parse_sequences(args.input, logger)

    if not sequence_data:
        logger.error("No sequence data parsed, terminating program.")
        return

    # Convert data
    positions_array, seq_info = convert_to_2d_array(sequence_data)
    origin_array = positions_array.copy()

    # Load genome sequence
    genome_seq = load_genome_sequence(args.genome, logger)

    if not genome_seq:
        logger.error("Genome sequence loading failed, terminating program.")
        return

    # Execute sequence extension and get results
    stopped_pairs = iterative_sequence_extension(
        positions_array,
        genome_seq,
        extend_length=args.extend_length,
        threshold=args.threshold,
        max_iterations=args.max_iterations,
        origin_array=origin_array
    )

    # Save results to CSV file
    save_results_to_csv(stopped_pairs, args.output)

    logger.info("Program execution completed.")

if __name__ == '__main__':
    main()