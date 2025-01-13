import argparse
import csv
from collections import defaultdict
import groupSeq


def format_feature_info(attrs):
    """
    格式化特征的详细信息
    返回一个包含所有重要信息的字典
    """
    info = {
        'name': None,
        'product': None,
        'note': None,
        'pseudo': False,
        'partial': False,
        'protein_id': None,
        'ec_number': None
    }

    # 获取基因名称（按优先级）
    if 'gene' in attrs:
        info['name'] = attrs['gene']
    elif 'Name' in attrs:
        info['name'] = attrs['Name']
    elif 'locus_tag' in attrs:
        info['name'] = attrs['locus_tag']

    # 获取产物信息
    if 'product' in attrs:
        info['product'] = attrs['product']

    # 获取注释信息
    if 'Note' in attrs:
        info['note'] = attrs['Note'].replace('%2C', ',').replace('%3B', ';')

    # 检查是否是假基因
    if 'pseudo' in attrs or 'pseudogene' in attrs:
        info['pseudo'] = True

    # 检查是否是部分基因
    if 'partial' in attrs:
        info['partial'] = True

    # 获取蛋白质ID
    if 'protein_id' in attrs:
        info['protein_id'] = attrs['protein_id']

    # 获取EC号
    if 'EC_number' in attrs:
        info['ec_number'] = attrs['EC_number']

    return info
def read_gff(gff_file):
    """
    读取GFF文件并返回所有特征信息字典
    """
    gff_dict = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # 解析属性字段
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=')
                    attr_dict[key.strip()] = value.strip()

            # 使用位置作为键
            gff_dict[(start, end)] = {
                'strand': strand,
                'attributes': attr_dict
            }

    return gff_dict
def is_completely_contained(seq_start, seq_end, feat_start, feat_end):
    """
    判断特征是否完全包含在序列中
    """
    return seq_start <= feat_start and seq_end >= feat_end


def load_intervals(csv_file, reference_genome):
    """
    加载所有Group部分的序列信息
    """
    group_sequences = defaultdict(list)

    with open(csv_file, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)

        in_group_section = False
        group_id = None

        for row in reader:
            if not row:
                continue

            # 检测组标题
            if row[0].startswith("Group"):
                in_group_section = True
                group_id = int(row[0].split()[1])
                continue

            # 遇到"Sequence Comparisons"时重置标志，等待下一个Group
            if row[0].strip() == "Sequence Comparisons":
                in_group_section = False
                group_id = None
                continue

            # 只在Group部分读取数据
            if in_group_section and group_id is not None:
                try:
                    start = int(row[0])
                    end = int(row[1])
                    length = int(row[2])

                    if start < 1 or end > len(reference_genome):
                        print(f"Warning: Interval {start}-{end} out of genome bounds. Skipping.")
                        continue

                    group_sequences[group_id].append({
                        'start': start,
                        'end': end,
                        'length': length,
                        'original_row': row
                    })
                except ValueError:
                    continue

    return group_sequences
def get_feature_name_with_direction_and_info(feat_info):
    """
    从特征信息中提取并格式化完整的注释信息
    """
    attrs = feat_info['attributes']
    info = format_feature_info(attrs)

    if not info['name']:
        return None

    # 构建注释字符串
    parts = []

    # 基本信息
    name_part = info['name']
    if info['pseudo']:
        name_part += '[pseudo]'
    if info['partial']:
        name_part += '[partial]'
    direction = '+' if feat_info['strand'] == '+' else '-'
    parts.append(f"{name_part}({direction})")

    # 产物信息
    if info['product']:
        parts.append(f"Product: {info['product']}")

    # 蛋白质ID
    if info['protein_id']:
        parts.append(f"Protein ID: {info['protein_id']}")

    # EC号
    if info['ec_number']:
        parts.append(f"EC: {info['ec_number']}")

    # 注释信息
    if info['note']:
        parts.append(f"Note: {info['note']}")

    return ' | '.join(parts)


def find_max_features(group_sequences, gff_dict):
    """
    找出所有序列中最大的特征数量
    """
    max_features = 0
    for sequences in group_sequences.values():
        for seq_data in sequences:
            feature_count = 0
            for (feat_start, feat_end), feat_info in gff_dict.items():
                if is_completely_contained(seq_data['start'], seq_data['end'], feat_start, feat_end):
                    if get_feature_name_with_direction_and_info(feat_info):
                        feature_count += 1
            max_features = max(max_features, feature_count)
    return max_features


def write_formatted_csv(results, fieldnames, output_csv):
    """
    将结果写入格式化的CSV文件
    """
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in results:
            # 如果是组标题行，添加分隔行
            if 'Group' in str(row['Start']):
                # 添加空行
                empty_row = {field: '' for field in fieldnames}
                writer.writerow(empty_row)
                # 添加组标题行
                writer.writerow(row)
            else:
                writer.writerow(row)


def annotate_sequences(input_csv, gff_file, reference_genome):
    """
    逐个group处理序列注释
    """
    gff_dict = read_gff(gff_file)

    # 读取原始CSV内容
    with open(input_csv, 'r', newline='', encoding='utf-8') as f:
        original_content = list(csv.reader(f))

    # 创建新的输出内容
    output_content = []
    current_group_rows = []
    in_group = False

    for row in original_content:
        if not row:  # 空行
            # 如果当前在处理group，先处理完当前group
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gff_dict, reference_genome)
                output_content.extend(processed_rows)
                current_group_rows = []
                in_group = False
            output_content.append(row)
            continue

        if row[0].startswith("Group"):  # 新的Group开始
            # 如果之前有未处理的group，先处理完
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gff_dict, reference_genome)
                output_content.extend(processed_rows)

            # 开始新的group
            current_group_rows = [row]
            in_group = True
            continue

        if row[0].strip() == "Sequence Comparisons":  # 遇到比较部分
            # 处理最后一个group
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gff_dict, reference_genome)
                output_content.extend(processed_rows)
            # 添加比较部分及后续所有内容
            output_content.append(row)
            in_group = False
            continue

        if in_group:  # 在Group处理过程中
            current_group_rows.append(row)
        else:  # 不在Group处理过程中的其他行直接添加
            output_content.append(row)

    # 处理最后一个group（如果有的话）
    if in_group and current_group_rows:
        processed_rows = process_group_rows(current_group_rows, gff_dict, reference_genome)
        output_content.extend(processed_rows)

    # 写回原始文件
    with open(input_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerows(output_content)

    print(f"注释完成，结果已更新至 {input_csv}")


def process_group_rows(group_rows, gff_dict, reference_genome):
    """
    处理单个group中的所有行
    """
    processed_rows = [group_rows[0]]  # 添加group标题行

    for row in group_rows[1:]:  # 跳过标题行
        try:
            start = int(row[0])
            end = int(row[1])

            # 查找此序列区间内的所有特征
            features = []
            for (feat_start, feat_end), feat_info in gff_dict.items():
                if is_completely_contained(start, end, feat_start, feat_end):
                    feature_name = get_feature_name_with_direction_and_info(feat_info)
                    if feature_name:
                        features.append(feature_name)

            # 在原始行后添加注释
            processed_rows.append(row + features)
        except ValueError:  # 非数据行直接添加
            processed_rows.append(row)

    return processed_rows


def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description='对序列进行GFF注释分析')

    parser.add_argument('-i', '--input',
                        required=True,
                        help='输入的CSV文件路径 (包含序列分组信息)')

    parser.add_argument('-g', '--gff',
                        required=True,
                        help='GFF注释文件路径')

    parser.add_argument('-f', '--fasta',
                        required=True,
                        help='参考基因组FASTA文件路径')


    return parser.parse_args()


def main():
    """
    主函数：处理命令行参数并执行注释
    """
    args = parse_arguments()

    try:
        reference_genome = groupSeq.load_genome_sequence(args.fasta)
    except Exception as e:
        print(f"错误：无法加载参考基因组文件 {args.fasta}")
        print(f"错误信息：{str(e)}")
        return

    try:
        # 移除output参数，直接修改输入文件
        annotate_sequences(args.input, args.gff, reference_genome)
    except Exception as e:
        print(f"错误：注释过程中出现异常")
        print(f"错误信息：{str(e)}")


if __name__ == "__main__":
    main()