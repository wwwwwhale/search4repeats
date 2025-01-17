#!/bin/bash

# 创建脚本生成目录
mkdir -p generated_scripts

# 基础脚本模板
generate_script() {
    local part_num=$1
    local output_file="generated_scripts/run_part${part_num}.sh"

    cat > "$output_file" << 'EOF'
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

# 确保 results 和日志目录存在
mkdir -p results
mkdir -p logs

# 最大并行任务数
MAX_PARALLEL=5

# 从文件读取文件列表
EOF

    # 添加特定的part数字
    echo "mapfile -t files < file_lists/part${part_num}.txt" >> "$output_file"

    # 继续添加剩余的脚本内容
    cat >> "$output_file" << 'EOF'
total_files=${#files[@]}
current=0
failed_tasks=()

# 处理函数
process_file() {
    local file=$1
    local basename=$(basename "$file" .gbff)
    local log_file="logs/${basename}.log"
    
    # 运行命令并记录日志
    if python run.py "$file" "results/$basename" > "$log_file" 2>&1; then
        echo "Successfully processed: $file"
    else
        echo "Error processing: $file. Check $log_file for details" >&2
        failed_tasks+=("$file")
    fi
}

# 等待有空闲槽位
wait_for_slot() {
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        # 等待任意一个子进程结束
         sleep 1
        # 不检查返回状态，让新任务继续
    done
}

# 主处理循环
for file in "${files[@]}"; do
    wait_for_slot
    
    # 显示进度
    current=$((current + 1))
    echo "Processing $current of $total_files: $file"
    
    # 启动处理并立即转入后台
    (process_file "$file") &
done

# 等待所有剩余任务完成
wait

# 报告结果
echo "All tasks completed!"
if [ ${#failed_tasks[@]} -eq 0 ]; then
    echo "All tasks succeeded!"
else
    echo "Failed tasks:"
    printf '%s\n' "${failed_tasks[@]}"
    echo "Total failed tasks: ${#failed_tasks[@]}"
fi
EOF

    # 设置执行权限
    chmod +x "$output_file"
    echo "Generated script: $output_file"
}

# 生成10个脚本
for i in {1..10}; do
    generate_script $i
done

echo "All scripts have been generated in the 'generated_scripts' directory!"