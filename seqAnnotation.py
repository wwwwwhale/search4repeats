import argparse
import csv
from collections import defaultdict
from Bio import SeqIO

def format_feature_info(feature):
    """
    从 GBFF 特征对象中提取并格式化特征信息

    Args:
        feature: SeqFeature对象

    Returns:
        dict: 包含所有重要信息的字典
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

    qualifiers = feature.qualifiers

    # 获取基因名称（按优先级）
    if 'gene' in qualifiers:
        info['name'] = qualifiers['gene'][0]
    elif 'locus_tag' in qualifiers:
        info['name'] = qualifiers['locus_tag'][0]

    # 获取产物信息
    if 'product' in qualifiers:
        info['product'] = qualifiers['product'][0]

    # 获取注释信息
    if 'note' in qualifiers:
        info['note'] = '; '.join(qualifiers['note'])

    # 检查是否是假基因
    if 'pseudo' in qualifiers or 'pseudogene' in qualifiers:
        info['pseudo'] = True

    # 检查是否是部分基因
    if feature.location.partial:
        info['partial'] = True

    # 获取蛋白质ID
    if 'protein_id' in qualifiers:
        info['protein_id'] = qualifiers['protein_id'][0]

    # 获取EC号
    if 'EC_number' in qualifiers:
        info['ec_number'] = qualifiers['EC_number'][0]

    return info

def load_gbff_data(gbff_file):
    """
    从GBFF文件加载序列和特征信息

    Args:
        gbff_file: GBFF文件路径

    Returns:
        tuple: (参考序列, 特征字典)
    """
    gbff_dict = {}
    reference_genome = None

    try:
        for record in SeqIO.parse(gbff_file, "genbank"):
            # 保存第一个序列作为参考序列
            if reference_genome is None:
                reference_genome = str(record.seq)
                print(f"成功加载参考序列，长度为 {len(reference_genome)} bp")

            # 处理特征
            for feature in record.features:
                if feature.type in ['CDS', 'gene', 'tRNA', 'rRNA']:
                    start = int(feature.location.start) + 1  # 转换为1-based坐标
                    end = int(feature.location.end)
                    strand = '+' if feature.location.strand == 1 else '-'

                    gbff_dict[(start, end)] = {
                        'strand': strand,
                        'feature_type': feature.type,
                        'qualifiers': feature.qualifiers
                    }

        print(f"成功加载 {len(gbff_dict)} 个特征注释")
        return reference_genome, gbff_dict

    except Exception as e:
        print(f"加载GBFF文件时发生错误: {str(e)}")
        return None, None

def is_completely_contained(seq_start, seq_end, feat_start, feat_end):
    """
    判断特征是否完全包含在序列中
    """
    return seq_start <= feat_start and seq_end >= feat_end

def get_feature_name_with_direction_and_info(feat_info):
    """
    从特征信息中提取并格式化完整的注释信息

    Args:
        feat_info: 特征信息字典

    Returns:
        str: 格式化的注释信息字符串
    """
    qualifiers = feat_info['qualifiers']
    parts = []

    # 获取基因名称
    name = None
    if 'gene' in qualifiers:
        name = qualifiers['gene'][0]
    elif 'locus_tag' in qualifiers:
        name = qualifiers['locus_tag'][0]

    if not name:
        return None

    # 基本信息
    name_part = name
    if 'pseudo' in qualifiers:
        name_part += '[pseudo]'
    direction = feat_info['strand']
    parts.append(f"{name_part}({direction})")

    # 产物信息
    if 'product' in qualifiers:
        parts.append(f"Product: {qualifiers['product'][0]}")

    # 蛋白质ID
    if 'protein_id' in qualifiers:
        parts.append(f"Protein ID: {qualifiers['protein_id'][0]}")

    # EC号
    if 'EC_number' in qualifiers:
        parts.append(f"EC: {qualifiers['EC_number'][0]}")

    # 注释信息
    if 'note' in qualifiers:
        parts.append(f"Note: {'; '.join(qualifiers['note'])}")

    return ' | '.join(parts)

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

            if row[0].startswith("Group"):
                in_group_section = True
                group_id = int(row[0].split()[1])
                continue

            if row[0].strip() == "Sequence Comparisons":
                in_group_section = False
                group_id = None
                continue

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

def annotate_sequences(input_csv, reference_genome, gbff_dict):
    """
    使用GBFF信息对序列进行注释
    """
    with open(input_csv, 'r', newline='', encoding='utf-8') as f:
        original_content = list(csv.reader(f))

    output_content = []
    current_group_rows = []
    in_group = False

    for row in original_content:
        if not row:
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gbff_dict, reference_genome)
                output_content.extend(processed_rows)
                current_group_rows = []
                in_group = False
            output_content.append(row)
            continue

        if row[0].startswith("Group"):
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gbff_dict, reference_genome)
                output_content.extend(processed_rows)

            current_group_rows = [row]
            in_group = True
            continue

        if row[0].strip() == "Sequence Comparisons":
            if in_group and current_group_rows:
                processed_rows = process_group_rows(current_group_rows, gbff_dict, reference_genome)
                output_content.extend(processed_rows)
            output_content.append(row)
            in_group = False
            continue

        if in_group:
            current_group_rows.append(row)
        else:
            output_content.append(row)

    if in_group and current_group_rows:
        processed_rows = process_group_rows(current_group_rows, gbff_dict, reference_genome)
        output_content.extend(processed_rows)

    with open(input_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerows(output_content)

    print(f"注释完成，结果已更新至 {input_csv}")

def process_group_rows(group_rows, gbff_dict, reference_genome):
    """
    处理单个group中的所有行，添加GBFF注释信息
    """
    processed_rows = [group_rows[0]]

    for row in group_rows[1:]:
        try:
            start = int(row[0])
            end = int(row[1])

            features = []
            for (feat_start, feat_end), feat_info in gbff_dict.items():
                if is_completely_contained(start, end, feat_start, feat_end):
                    feature_name = get_feature_name_with_direction_and_info(feat_info)
                    if feature_name:
                        features.append(feature_name)

            processed_rows.append(row + features)
        except ValueError:
            processed_rows.append(row)

    return processed_rows

def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description='对序列进行GBFF注释分析')

    parser.add_argument('-i', '--input',
                        required=True,
                        help='输入的CSV文件路径 (包含序列分组信息)')

    parser.add_argument('-g', '--gbff',
                        required=True,
                        help='GBFF文件路径 (包含序列和注释信息)')

    return parser.parse_args()

def main():
    """
    主函数：处理命令行参数并执行注释
    """
    args = parse_arguments()

    try:
        # 从GBFF文件加载序列和特征信息
        reference_genome, gbff_dict = load_gbff_data(args.gbff)
        if reference_genome is None or gbff_dict is None:
            print("错误：无法从GBFF文件加载数据")
            return

        # 执行注释
        annotate_sequences(args.input, reference_genome, gbff_dict)

    except Exception as e:
        print(f"错误：程序执行过程中出现异常")
        print(f"错误信息：{str(e)}")

if __name__ == "__main__":
    main()