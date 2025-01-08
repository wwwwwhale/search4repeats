import csv
from collections import defaultdict
import checkSimilarity


def read_gff(gff_file):
    """
    读取GFF文件并返回基因注释信息的字典
    """
    gff_dict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            # 只关注gene类型的特征
            if feature_type == 'gene':
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=')
                        attr_dict[key] = value

                gff_dict[(start, end)] = {
                    'strand': strand,
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




def annotate_sequences(input_csv, gff_file, output_csv,reference_genome):
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
            'End': '',
            'Complete_Overlap_Genes': '',
            'Partial_Overlap_Genes': ''
        })

        # 处理组内的每个序列
        for seq_data in group_sequences[group_id]:
            start = seq_data['start']
            end = seq_data['end']

            # 存储完全重叠和部分重叠的基因
            complete_overlaps = []
            partial_overlaps = []

            # 检查与每个基因的重叠情况
            for (gene_start, gene_end), gene_info in gff_dict.items():
                # 判断重叠类型
                if start <= gene_start and end >= gene_end:
                    # 完全重叠
                    gene_name = gene_info['attributes'].get('Name', 'Unknown')
                    complete_overlaps.append(f"{gene_name}({gene_info['strand']})")
                elif (start <= gene_end and end >= gene_end) or \
                        (start <= gene_start and end >= gene_start):
                    # 部分重叠
                    gene_name = gene_info['attributes'].get('Name', 'Unknown')
                    partial_overlaps.append(f"{gene_name}({gene_info['strand']})")

            # 保存该序列的注释结果
            results.append({
                'Start': start,
                'End': end,
                'Complete_Overlap_Genes': ';'.join(complete_overlaps) if complete_overlaps else 'None',
                'Partial_Overlap_Genes': ';'.join(partial_overlaps) if partial_overlaps else 'None'
            })

    # 将结果写入CSV文件
    with open(output_csv, 'w', newline='') as f:
        fieldnames = ['Start', 'End', 'Complete_Overlap_Genes', 'Partial_Overlap_Genes']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"注释完成，结果已保存至 {output_csv}")



# 使用示例
if __name__ == "__main__":
    input_csv = "results/sequential_grouping_results.csv"  # 输入的CSV文件
    gff_file = "data/GCA_000005845.2_ASM584v2_genomic.gff"  # GFF注释文件
    output_csv = "results/annotated_sequences_1950.csv"  # 输出的CSV文件
    fasta_file = "data/GCA_000005845.2_ASM584v2_genomic.fna"
    reference_genome = checkSimilarity.load_genome_sequence(fasta_file)
    annotate_sequences(input_csv, gff_file, output_csv,reference_genome)