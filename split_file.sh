#!/bin/bash

# 确保输出目录存在
mkdir -p file_lists

# 获取所有 .gbff 文件并保存到数组
files=(data/*.gbff)
total_files=${#files[@]}

# 计算每份大约应该有多少文件
files_per_part=$((total_files / 3))
remainder=$((total_files % 3))

# 清空或创建新的文件列表
> file_lists/part1.txt
> file_lists/part2.txt
> file_lists/part3.txt

# 分配文件
current=0

# 分配第一部分
end1=$((files_per_part + (remainder > 0 ? 1 : 0)))
for ((i=0; i<end1; i++)); do
    echo "${files[current]}" >> file_lists/part1.txt
    current=$((current + 1))
done

# 分配第二部分
end2=$((current + files_per_part + (remainder > 1 ? 1 : 0)))
for ((i=current; i<end2; i++)); do
    echo "${files[current]}" >> file_lists/part2.txt
    current=$((current + 1))
done

# 分配第三部分
for ((i=current; i<total_files; i++)); do
    echo "${files[current]}" >> file_lists/part3.txt
    current=$((current + 1))
done

# 显示分配结果
echo "Total files: $total_files"
echo "Files in part1: $(wc -l < file_lists/part1.txt)"
echo "Files in part2: $(wc -l < file_lists/part2.txt)"
echo "Files in part3: $(wc -l < file_lists/part3.txt)"

# 显示每个文件的内容预览
echo -e "\nPart 1 preview:"
head -n 3 file_lists/part1.txt

echo -e "\nPart 2 preview:"
head -n 3 file_lists/part2.txt

echo -e "\nPart 3 preview:"
head -n 3 file_lists/part3.txt
