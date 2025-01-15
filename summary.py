import csv
from collections import defaultdict


def process_csv_and_group(input_file, output_file):
    """
    读取并处理CSV文件，统计每组的区间数量，并按区间数量分组。

    :param input_file: 输入的CSV文件路径
    :param output_file: 输出的CSV文件路径
    """
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

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
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Interval Count", "Groups"])  # 写入表头
        for count, groups in sorted(grouped_by_intervals.items()):
            writer.writerow([count, "; ".join(groups)])  # 用分号分隔组号

    print(f"处理完成！结果已保存到 {output_file}")


# 调用函数
input_csv = "./Jeju/Jeju_grouping_results.csv"  # 输入文件路径
output_csv = "./Jeju/Jeju_summary.csv"  # 输出文件路径
process_csv_and_group(input_csv, output_csv)
