#!/bin/bash
#SBATCH -J repeat_seq
#SBATCH -p cnall
#SBATCH -N 1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=10
#SBATCH -o %j.out
#SBATCH -e %j.err

source /home/zouqin/anaconda3/etc/profile.d/conda.sh
conda activate FindRepeat


# 确保 results 目录存在
mkdir -p results

# 最大并行任务数
MAX_PARALLEL=5

# 从文件读取文件列表
mapfile -t files < file_lists/part1.txt
total_files=${#files[@]}
current=0



# 处理函数
process_file() {
    local file=$1
    local basename=$(basename "$file" .gbff)
    python run.py "$file" "results/$basename"
}

# 等待有空闲槽位
wait_for_slot() {
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
}

# 主处理循环
for file in "${files[@]}"; do
    wait_for_slot  # 等待有空闲槽位
    
    # 显示进度
    current=$((current + 1))
    echo "Processing $current of $total_files: $file"
    
    # 启动处理
    process_file "$file" &
done

# 等待所有任务完成
wait

echo "All tasks completed!"