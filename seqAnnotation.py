import csv
from collections import defaultdict
import checkSimilarity


def read_gff(gff_file):
    """
    读取GFF文件并返回所有特征类型的注释信息字典
    参数:
        gff_file: GFF格式文件的路径
    返回:
        包含所有特征注释信息的嵌套字典，以位置为键
    """
    gff_dict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            score = parts[5]
            strand = parts[6]
            phase = parts[7]
            attributes = parts[8]

            # 解析属性字段
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=')
                    attr_dict[key.strip()] = value.strip()

            # 使用位置作为键，以便于重叠检查
            gff_dict[(start, end)] = {
                'seqid': seqid,
                'source': source,
                'feature_type': feature_type,
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attr_dict
            }

    return gff_dict

def classify_overlap(seq_start, seq_end, gene_start, gene_end):
    """
    判断重叠类型
    返回: 'complete' 表示完全重叠, 'partial' 表示部分重叠, None 表示无重叠
    """
    if seq_start <= gene_start and seq_end >= gene_end:
        return 'complete'
    elif seq_end < gene_start or seq_start > gene_end:
        return None
    else:
        return 'partial'


def annotate_sequences(input_csv, gff_file, output_csv, reference_genome):
    """
    主函数：读取CSV文件，注释序列，并输出结果
    """
    # 读取GFF文件
    gff_dict = read_gff(gff_file)

    # 使用load_intervals读取分组数据
    group_sequences = checkSimilarity.load_intervals(input_csv, reference_genome)

    # 存储结果
    results = []

    # 处理每个分组
    for group_id in sorted(group_sequences.keys()):
        # 添加组标题
        results.append({
            'Start': f'Group {group_id}',
            'End': 'End',
            'Length': 'Length',
            'Seqid': '',
            'Feature_Type': '',
            'Strand': '',
            'Complete_Overlap_Features': '',
            'Partial_Overlap_Features': '',
            'Attributes': ''
        })

        # 处理组内的每个序列
        for seq_data in group_sequences[group_id]:
            start = seq_data['start']
            end = seq_data['end']
            length = end - start

            # 存储完全重叠和部分重叠的特征
            complete_overlaps = []
            partial_overlaps = []
            seqids = set()
            feature_types = set()
            strands = set()
            all_attributes = set()

            # 检查与每个特征的重叠情况
            for (feat_start, feat_end), feat_info in gff_dict.items():
                overlap_type = classify_overlap(start, end, feat_start, feat_end)

                if overlap_type:  # 如果有重叠（完全或部分）
                    seqids.add(feat_info['seqid'])
                    feature_types.add(feat_info['feature_type'])
                    strands.add(feat_info['strand'])

                    # 格式化特征信息
                    feature_detail = f"{feat_info['feature_type']}:{feat_info['attributes'].get('Name', 'Unknown')}"

                    # 收集属性信息
                    attrs = feat_info['attributes']
                    attr_str = ';'.join(f"{k}={v}" for k, v in attrs.items())
                    all_attributes.add(attr_str)

                    if overlap_type == 'complete':
                        complete_overlaps.append(feature_detail)
                    else:  # partial
                        partial_overlaps.append(feature_detail)

            # 保存该序列的注释结果
            results.append({
                'Start': start,
                'End': end,
                'Length': length,
                'Seqid': ';'.join(seqids),
                'Feature_Type': ';'.join(feature_types),
                'Strand': ';'.join(strands),
                'Complete_Overlap_Features': ';'.join(complete_overlaps) if complete_overlaps else 'None',
                'Partial_Overlap_Features': ';'.join(partial_overlaps) if partial_overlaps else 'None',
                'Attributes': '||'.join(all_attributes) if all_attributes else 'None'
            })

    # 将结果写入CSV文件
    with open(output_csv, 'w', newline='') as f:
        fieldnames = ['Start', 'End', 'Length', 'Seqid', 'Feature_Type',
                      'Strand', 'Complete_Overlap_Features',
                      'Partial_Overlap_Features', 'Attributes']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"注释完成，结果已保存至 {output_csv}")


def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description='对序列进行GFF注释分析')

    # 必需参数
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
    # 解析命令行参数
    args = parse_arguments()

    # 加载参考基因组
    try:
        reference_genome = checkSimilarity.load_genome_sequence(args.fasta)
    except Exception as e:
        print(f"错误：无法加载参考基因组文件 {args.fasta}")
        print(f"错误信息：{str(e)}")
        return

    # 执行注释
    try:
        annotate_sequences(args.input, args.gff, args.output, reference_genome)
        print(f"注释完成！结果已保存至：{args.output}")
    except Exception as e:
        print(f"错误：注释过程中出现异常")
        print(f"错误信息：{str(e)}")


if __name__ == "__main__":
    main()