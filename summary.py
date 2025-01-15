import csv
from collections import defaultdict
import argparse
import sys

def process_csv_and_group(input_file, output_file):
    """
    读取并处理CSV文件，统计每组的区间数量，并按区间数量分组。

    :param input_file: 输入的CSV文件路径
    :param output_file: 输出的CSV文件路径
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时发生错误：{str(e)}")
        sys.exit(1)

    group_data = {}  # 保存组信息，键为组号，值为区间数量
    current_group = None
    current_group_count = 0
    in_interval_section = False  # 标记是否在区间部分

    for line in lines:
        line = line.strip()

        if line.startswith("Group"):  # 遇到新组的开头
            if current_group is not None:
                group_data[current_group] = current_group_count
            current_group = line  # 更新当前组号
            current_group_count = 0  # 重置区间计数
            in_interval_section = True  # 开始读取区间部分

        elif in_interval_section:
            if line == "":  # 遇到空行，结束当前组的区间部分
                in_interval_section = False
            elif not line.startswith("Start"):  # 不是标题行，则计为区间
                current_group_count += 1

    # 保存最后一组数据
    if current_group is not None:
        group_data[current_group] = current_group_count

    # 按区间数量分组
    grouped_by_intervals = defaultdict(list)
    for group, count in group_data.items():
        grouped_by_intervals[count].append(group)

    # 将结果写入输出文件
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            # 写入总组数
            writer.writerow([f"Total number of groups: {len(group_data)}"])
            writer.writerow([])  # 添加一个空行
            writer.writerow(["Interval Count", "Groups"])  # 写入表头
            for count, groups in sorted(grouped_by_intervals.items()):
                writer.writerow([count, "; ".join(groups)])  # 用分号分隔组号
        print(f"处理完成！结果已保存到 {output_file}")
    except Exception as e:
        print(f"写入文件时发生错误：{str(e)}")
        sys.exit(1)

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='处理CSV文件并按区间数量分组')
    parser.add_argument('-i', '--input', required=True, help='输入CSV文件的路径')
    parser.add_argument('-o', '--output', required=True, help='输出CSV文件的路径')

    # 解析命令行参数
    args = parser.parse_args()

    # 调用处理函数
    process_csv_and_group(args.input, args.output)

if __name__ == "__main__":
    main()