import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO


class GroupAnalyzer:
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.folder_name = folder_path.name

        # 查找对应的文件
        self.grouped_results_file = list(folder_path.glob("*grouped*results.csv"))
        self.summary_file = list(folder_path.glob("*summary.csv"))
        self.gbff_file = list(Path("./data").glob(f"{folder_path.name}*.gbff"))  # 添加gbff文件查找

        # 检查文件是否存在
        self.files_exist = (len(self.grouped_results_file) > 0 and
                            len(self.summary_file) > 0)

        if self.files_exist:
            self.grouped_results_file = self.grouped_results_file[0]
            self.summary_file = self.summary_file[0]

        # 获取gbff文件路径
        self.gbff_path = self.gbff_file[0] if self.gbff_file else None
        self.genome_length = 0  # 初始化基因组长度

    def load_genome_sequence(self):
        """从GBFF文件加载基因组序列并获取长度"""
        if not self.gbff_path:
            print(f"警告: {self.folder_name} 中未找到GBFF文件")
            return 0

        try:
            for record in SeqIO.parse(self.gbff_path, "genbank"):
                print(f"读取基因组序列: {record.id}")
                sequence_length = len(record.seq)
                print(f"序列长度: {sequence_length} bp")

                # 记录额外的序列信息
                if hasattr(record, 'annotations'):
                    source = record.annotations.get('source', 'Unknown source')
                    print(f"来源: {source}")

                return sequence_length

            print(f"警告: {self.gbff_path} 中未找到序列")
            return 0

        except Exception as e:
            print(f"读取基因组序列时出错 {self.gbff_path}: {str(e)}")
            return 0

    def read_grouped_results(self):
        """读取分组结果文件，获取序列长度信息"""
        try:
            with open(self.grouped_results_file, 'r') as file:
                lines = file.readlines()

            groups = {}
            current_group = None
            lengths = []

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('Group'):
                    if current_group is not None and lengths:
                        groups[current_group] = np.mean(lengths)
                    current_group = int(line.split()[1])
                    lengths = []
                    continue

                if line.startswith('Start,') or line.startswith('Seq1'):
                    continue

                parts = line.split(',')
                if len(parts) >= 3:
                    try:
                        length = int(parts[2])
                        if 0 < length < 100000:
                            lengths.append(length)
                    except ValueError:
                        continue

            if current_group is not None and lengths:
                groups[current_group] = np.mean(lengths)

            return groups

        except Exception as e:
            print(f"读取文件 {self.folder_name} 的分组结果时出错: {str(e)}")
            return {}

    def read_summary(self):
        """读取摘要文件，获取每个组的序列数量"""
        try:
            with open(self.summary_file, 'r') as file:
                lines = file.readlines()

            group_counts = {}

            for line in lines:
                line = line.strip()
                if not line or line.startswith('Total') or line.startswith('Interval Count'):
                    continue

                parts = line.split(',')
                if len(parts) >= 3:
                    interval_count = int(parts[0])
                    groups_str = parts[2]
                    groups = groups_str.split(';')
                    for group in groups:
                        group = group.strip()
                        if 'Group' in group:
                            group_num = int(group.split()[1])
                            group_counts[group_num] = interval_count

            return group_counts

        except Exception as e:
            print(f"读取文件 {self.folder_name} 的摘要时出错: {str(e)}")
            return {}

    def analyze(self):
        """分析文件并返回结果"""
        if not self.files_exist:
            print(f"跳过 {self.folder_name}: 缺少必要的文件")
            return pd.DataFrame()

        try:
            # 获取基因组长度
            self.genome_length = self.load_genome_sequence()

            group_lengths = self.read_grouped_results()
            group_counts = self.read_summary()

            # Debug信息
            if group_lengths:
                print(f"\n{self.folder_name} 的长度信息:")
                for group, length in group_lengths.items():
                    print(f"Group {group}: 平均长度 = {length:.2f}")

            results = []
            all_groups = sorted(set(list(group_lengths.keys()) + list(group_counts.keys())))

            for group_num in all_groups:
                results.append({
                    'folder_name': self.folder_name,
                    'group_number': group_num,
                    'sequence_count': group_counts.get(group_num, 0),
                    'average_length': round(group_lengths.get(group_num, 0), 2),
                    'genome_length': self.genome_length  # 添加基因组长度
                })

            return pd.DataFrame(results)
        except Exception as e:
            print(f"分析文件 {self.folder_name} 时出错: {str(e)}")
            return pd.DataFrame()


def process_all_folders(base_path="./results"):
    """处理results目录下的所有文件夹"""
    base_path = Path(base_path)
    all_results = []

    folders = [f for f in base_path.iterdir() if f.is_dir()]
    print(f"开始处理 {len(folders)} 个文件夹...")

    for folder in folders:
        print(f"\n正在处理: {folder.name}")
        analyzer = GroupAnalyzer(folder)
        results = analyzer.analyze()

        if not results.empty:
            all_results.append(results)

    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        return final_results
    else:
        print("没有找到有效的数据")
        return pd.DataFrame()


def main():
    results = process_all_folders()

    if not results.empty:
        output_file = "all_groups_analysis.csv"
        results.to_csv(output_file, index=False)
        print(f"\n结果已保存到: {output_file}")

        print("\n统计摘要:")
        print("-" * 60)
        print(f"总文件夹数: {len(results['folder_name'].unique())}")
        print(f"总组数: {len(results)}")
        print(f"平均序列数量: {results['sequence_count'].mean():.2f}")
        print(f"平均序列长度: {results['average_length'].mean():.2f}")
        print(f"平均基因组长度: {results['genome_length'].mean():.2f}")
    else:
        print("未能生成有效的分析结果")


if __name__ == "__main__":
    main()