import pandas as pd
import re
import os
from pathlib import Path
from datetime import datetime

def extract_numbers_from_filename(filename):
    """从文件名中提取所有精确到小数点后3位的数字"""
    # 匹配小数点后3位的数字
    pattern = r'\d+\.\d{3}'
    numbers = re.findall(pattern, filename)
    return [float(n) for n in numbers]

def process_csv_file(filepath, return_details=False):
    """处理单个CSV文件
    
    Args:
        filepath: CSV文件路径
        return_details: 是否返回详细的spine数据
    
    Returns:
        如果return_details=False，返回汇总数据字典
        如果return_details=True，返回(汇总数据字典, 详细spine数据列表)
    """
    filename = os.path.basename(filepath)
    
    # 提取文件名中的数字
    numbers = extract_numbers_from_filename(filename)
    
    # 判断是否为特殊文件
    is_special = filename in ["2-1 63x 1-3 1 35.987 0.932.csv", "2-1 63x 1-3 2 17.226 0.924.csv"]
    
    # 读取CSV文件
    df = pd.read_csv(filepath)
    
    # 获取Length列的数据（跳过表头）
    length_data = df['Length'].tolist()
    
    if is_special:
        # 特殊文件：从第一行开始每4个就是spine数据
        # dendrite width是文件名中较小的数字，dendrite length是较大的数字
        if len(numbers) >= 2:
            dendrite_length = max(numbers)
            dendrite_width = min(numbers)
        else:
            dendrite_length = numbers[0] if numbers else None
            dendrite_width = None
        
        # 从第一行开始处理
        spine_data = length_data
        spine_count = len(spine_data) // 4
        
    else:
        # 普通文件：第一行是dendrite width，后面每4个是spine数据
        if len(length_data) == 0:
            return None if not return_details else (None, [])
        
        dendrite_width = length_data[0]
        spine_data = length_data[1:]  # 跳过第一行
        spine_count = len(spine_data) // 4
        
        # dendrite length从文件名中提取（应该是最大的那个数字，或者唯一的小数点后3位的数字）
        if numbers:
            dendrite_length = max(numbers) if len(numbers) > 1 else numbers[0]
        else:
            dendrite_length = None
    
    # 计算每个spine的详细数据
    spine_details = []
    spine_lengths = []
    
    for i in range(spine_count):
        idx = i * 4
        if idx + 3 < len(spine_data):  # 确保有完整的4个数据
            shaft_length = spine_data[idx]
            shaft_width = spine_data[idx + 1]
            head_length = spine_data[idx + 2]
            head_width = spine_data[idx + 3]
            spine_length = shaft_length + head_length
            
            spine_lengths.append(spine_length)
            
            if return_details:
                spine_details.append({
                    'filename': filename,
                    'spine_number': i + 1,
                    'shaft_length': shaft_length,
                    'shaft_width': shaft_width,
                    'head_length': head_length,
                    'head_width': head_width,
                    'spine_length': spine_length,
                    'dendrite_length': dendrite_length,
                    'dendrite_width': dendrite_width
                })
    
    # 计算平均spine length
    avg_spine_length = sum(spine_lengths) / len(spine_lengths) if spine_lengths else 0
    
    # 计算总spine length
    total_spine_length = sum(spine_lengths)
    
    summary = {
        'filename': filename,
        'dendrite_length': dendrite_length,
        'dendrite_width': dendrite_width,
        'spine_numbers': spine_count,
        'avg_spine_length': avg_spine_length,
        'total_spine_length': total_spine_length
    }
    
    if return_details:
        return summary, spine_details
    else:
        return summary

def main():
    # 获取当前目录下所有CSV文件
    current_dir = Path('.')
    csv_files = list(current_dir.glob('*.csv'))
    
    # 排除结果文件本身
    csv_files = [f for f in csv_files if f.name not in ['dendrite_analysis_results.csv', 'dendrite_spine_details.csv']]
    
    # 处理所有CSV文件
    results = []
    all_spine_details = []
    
    for csv_file in csv_files:
        try:
            # 获取汇总数据
            result = process_csv_file(csv_file, return_details=False)
            if result:
                results.append(result)
            
            # 获取详细spine数据
            summary, spine_details = process_csv_file(csv_file, return_details=True)
            if spine_details:
                all_spine_details.extend(spine_details)
                
        except Exception as e:
            print(f"处理文件 {csv_file} 时出错: {e}")
    
    # 创建汇总数据表
    df_results = pd.DataFrame(results)
    
    # 创建详细spine数据表
    df_spine_details = pd.DataFrame(all_spine_details)
    
    # 保存结果
    summary_file = 'dendrite_analysis_results.csv'
    details_file = 'dendrite_spine_details.csv'
    
    # 尝试保存，如果文件被占用则使用带时间戳的文件名
    try:
        df_results.to_csv(summary_file, index=False, encoding='utf-8-sig')
        df_spine_details.to_csv(details_file, index=False, encoding='utf-8-sig')
    except PermissionError:
        # 如果文件被占用，使用带时间戳的文件名
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        summary_file = f'dendrite_analysis_results_{timestamp}.csv'
        details_file = f'dendrite_spine_details_{timestamp}.csv'
        df_results.to_csv(summary_file, index=False, encoding='utf-8-sig')
        df_spine_details.to_csv(details_file, index=False, encoding='utf-8-sig')
        print(f"\n注意：原文件被占用，已使用新文件名保存")
    
    print(f"处理完成！共处理 {len(results)} 个文件")
    print(f"共提取 {len(all_spine_details)} 个spine的详细数据")
    print(f"\n汇总结果已保存到: {summary_file}")
    print(f"详细spine数据已保存到: {details_file}")
    print("\n汇总数据预览:")
    print(df_results.to_string())
    print("\n详细spine数据预览（前10行）:")
    print(df_spine_details.head(10).to_string())
    
    return df_results, df_spine_details

if __name__ == '__main__':
    main()

