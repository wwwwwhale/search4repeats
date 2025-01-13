import csv

def read_intervals_from_csv(csv_file_path):
    """
    读取CSV文件，返回除第一列外的区间和length。

    Args:
        csv_file_path (str): CSV文件的路径

    Returns:
        list: 包含区间和length的元组列表，每个元组的形式为 ((current_start1, current_end1), (current_start2, current_end2), length)
    """
    seqpair = []

    with open(csv_file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # 跳过第一行（标题行）

        for row in reader:
            # 跳过空行
            if not row:
                continue

            # 获取当前行数据，排除第一列
            current_start1, current_end1, current_start2, current_end2, length = map(int, row[1:])

            # 每行保存两个区间和一个length，并添加到results列表
            seqpair.append(((current_start1, current_end1), (current_start2, current_end2), length))

        # 按照length值进行排序
        seqpair.sort(key=lambda x: x[2])

    return seqpair


def calculate_overlap(interval_a, interval_b):
    """
    计算两个区间的重叠比例。
    :param interval_a: 第一个区间 (start, end)
    :param interval_b: 第二个区间 (start, end)
    :return: 如果重叠比例 >= 90%，返回 True，否则返回 False
    """
    start_a, end_a = interval_a
    start_b, end_b = interval_b
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    overlap_length = max(0, overlap_end - overlap_start)  # 重叠长度
    total_length = max(end_a - start_a, end_b - start_b)  # 两区间的最大长度
    return (overlap_length / total_length) >= 0.85



def group_intervals_by_overlap(seqpairs):
    """
    将区间按照重叠逻辑进行分组。
    :param seqpairs: 包含多个区间的列表，每个元素是 ((interval1), (interval2), length)
    :return: 分组结果，列表形式
    """
    grouped_results = []  # 用于存储分组的结果

    for seqpair in seqpairs:
        interval1, interval2, length = seqpair
        added_to_group = False  # 标记是否已加入某个分组

        # 遍历现有分组，尝试将当前区间加入分组
        for group in grouped_results:
            if any(
                calculate_overlap(interval1, g[0]) or calculate_overlap(interval1, g[1]) or
                calculate_overlap(interval2, g[0]) or calculate_overlap(interval2, g[1])
                for g in group
            ):
                group.append(seqpair)
                added_to_group = True
                break

        # 如果没有加入任何分组，创建一个新的分组
        if not added_to_group:
            grouped_results.append([seqpair])

    return grouped_results


def main():
    csv_file_path = "./ecoli/ecoli_similar_seq.csv"
    seqpairs = read_intervals_from_csv(csv_file_path)  # 读取区间对
    grouped_results = group_intervals_by_overlap(seqpairs)  # 分组区间

    # 打印分组结果
    print("分组结果:")
    for i, group in enumerate(grouped_results, start=1):
        print(f"分组 {i}:")
        for seq in group:
            print(seq)

if __name__ == "__main__":
    main()