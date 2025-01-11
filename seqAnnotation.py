import argparse
import csv
from collections import defaultdict
import checkSimilarity


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


def annotate_sequences(input_csv, gff_file, output_csv, reference_genome):
    """
    主函数：读取CSV文件，注释序列，并输出结果
    """
    # 读取GFF文件
    gff_dict = read_gff(gff_file)

    # 读取分组数据
    group_sequences = checkSimilarity.load_intervals(input_csv, reference_genome)

    # 找出最大特征数量
    max_features = find_max_features(group_sequences, gff_dict)

    # 准备列名（动态生成特征列）
    fieldnames = ['Start', 'End', 'Length'] + [f'Feature_{i + 1}' for i in range(max_features)]

    # 存储结果
    results = []

    # 处理每个分组
    for group_id in sorted(group_sequences.keys()):
        # 添加组标题
        group_row = {field: '' for field in fieldnames}
        group_row['Start'] = f'Group {group_id}'
        results.append(group_row)

        # 处理组内的每个序列
        for seq_data in group_sequences[group_id]:
            start = seq_data['start']
            end = seq_data['end']

            # 收集该序列的所有特征
            features = []
            for (feat_start, feat_end), feat_info in gff_dict.items():
                if is_completely_contained(start, end, feat_start, feat_end):
                    feature_name = get_feature_name_with_direction_and_info(feat_info)
                    if feature_name:
                        features.append(feature_name)

            # 创建结果行
            result_row = {field: '' for field in fieldnames}
            result_row['Start'] = start
            result_row['End'] = end
            result_row['Length'] = end - start

            # 填充特征
            for i, feature in enumerate(sorted(features)):
                result_row[f'Feature_{i + 1}'] = feature

            results.append(result_row)

    # 将结果写入CSV文件
    write_formatted_csv(results, fieldnames, output_csv)
    print(f"注释完成，结果已保存至 {output_csv}")


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

    parser.add_argument('-o', '--output',
                        required=True,
                        help='输出的CSV文件路径')

    return parser.parse_args()


def main():
    """
    主函数：处理命令行参数并执行注释
    """
    args = parse_arguments()

    try:
        reference_genome = checkSimilarity.load_genome_sequence(args.fasta)
    except Exception as e:
        print(f"错误：无法加载参考基因组文件 {args.fasta}")
        print(f"错误信息：{str(e)}")
        return

    try:
        annotate_sequences(args.input, args.gff, args.output, reference_genome)
    except Exception as e:
        print(f"错误：注释过程中出现异常")
        print(f"错误信息：{str(e)}")


if __name__ == "__main__":
    main()