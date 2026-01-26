#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVL和E核团分析结果可视化脚本
在同一图表中展示两个核团的数据，并添加统计学分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import os
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置绘图风格
sns.set_style("whitegrid")
sns.set_palette("husl")

def perform_statistical_test(data1, data2, param_name):
    """
    执行统计学检验
    返回：统计方法名称，统计量，p值
    """
    # 移除NaN值
    data1_clean = data1.dropna()
    data2_clean = data2.dropna()
    
    if len(data1_clean) < 3 or len(data2_clean) < 3:
        return None, None, None
    
    # 检查正态性（Shapiro-Wilk test，样本量小于50时使用）
    if len(data1_clean) <= 50 and len(data2_clean) <= 50:
        _, p_norm1 = stats.shapiro(data1_clean)
        _, p_norm2 = stats.shapiro(data2_clean)
        is_normal = p_norm1 > 0.05 and p_norm2 > 0.05
    else:
        # 样本量大时假设正态分布
        is_normal = True
    
    # 检查方差齐性
    _, p_var = stats.levene(data1_clean, data2_clean)
    equal_var = p_var > 0.05
    
    # 选择适当的检验方法
    if is_normal:
        # 使用t检验
        if equal_var:
            stat, p_value = stats.ttest_ind(data1_clean, data2_clean)
            method = "独立样本t检验"
        else:
            stat, p_value = stats.ttest_ind(data1_clean, data2_clean, equal_var=False)
            method = "Welch's t检验"
    else:
        # 使用非参数检验
        stat, p_value = stats.mannwhitneyu(data1_clean, data2_clean, alternative='two-sided')
        method = "Mann-Whitney U检验"
    
    return method, stat, p_value

def get_significance_stars(p_value):
    """根据p值返回显著性标记"""
    if p_value is None or np.isnan(p_value):
        return "ns"
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"

def create_comparison_plots(df_mvl, df_e, output_prefix="mvl_e_comparison"):
    """
    创建MVL和E核团的对比可视化图表
    """
    # 获取数值型列（排除原始参数列，只保留分析结果列）
    exclude_cols = ['a', 'b', 'theta', 'xc', 'yc']
    numeric_cols = [col for col in df_mvl.columns 
                   if col not in exclude_cols and df_mvl[col].dtype in [np.number]]
    
    if len(numeric_cols) == 0:
        print("未找到可分析的数值型列")
        return
    
    print(f"找到 {len(numeric_cols)} 个分析参数")
    print(f"参数列表: {numeric_cols}")
    
    # 准备数据
    plot_data = []
    stats_results = []
    
    for col in numeric_cols:
        mvl_data = df_mvl[col].dropna()
        e_data = df_e[col].dropna()
        
        # 添加数据到绘图列表
        for val in mvl_data:
            plot_data.append({'Parameter': col, 'Value': val, 'Region': 'MVL'})
        for val in e_data:
            plot_data.append({'Parameter': col, 'Value': val, 'Region': 'E'})
        
        # 执行统计学检验
        method, stat, p_value = perform_statistical_test(df_mvl[col], df_e[col], col)
        stars = get_significance_stars(p_value)
        
        stats_results.append({
            'Parameter': col,
            'Method': method,
            'Statistic': stat,
            'P_value': p_value,
            'Significance': stars,
            'MVL_mean': mvl_data.mean(),
            'MVL_std': mvl_data.std(),
            'MVL_median': mvl_data.median(),
            'E_mean': e_data.mean(),
            'E_std': e_data.std(),
            'E_median': e_data.median(),
            'MVL_n': len(mvl_data),
            'E_n': len(e_data)
        })
    
    plot_df = pd.DataFrame(plot_data)
    stats_df = pd.DataFrame(stats_results)
    
    # 保存统计结果
    stats_df.to_csv(f"{output_prefix}_statistics.csv", index=False, encoding='utf-8-sig')
    print(f"\n统计结果已保存到: {output_prefix}_statistics.csv")
    
    # 创建多个子图
    n_params = len(numeric_cols)
    n_cols = 3
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_params == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, col in enumerate(numeric_cols):
        ax = axes[idx]
        
        # 准备数据
        mvl_vals = df_mvl[col].dropna()
        e_vals = df_e[col].dropna()
        
        # 获取统计结果
        stat_row = stats_df[stats_df['Parameter'] == col].iloc[0]
        p_val = stat_row['P_value']
        stars = stat_row['Significance']
        method = stat_row['Method']
        
        # 创建箱线图和小提琴图组合
        positions = [1, 2]
        data_to_plot = [mvl_vals.values, e_vals.values]
        
        # 小提琴图
        parts = ax.violinplot(data_to_plot, positions=positions, widths=0.6, 
                             showmeans=True, showmedians=True)
        
        # 设置颜色（使用用户指定的颜色）
        # MVL: #4d7e54, E: #2c4ca0
        for pc, color in zip(parts['bodies'], ['#4d7e54', '#2c4ca0']):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans', 'cmedians'):
            if partname in parts:
                parts[partname].set_edgecolor('black')
                parts[partname].set_linewidth(1.5)
        
        # 添加散点图显示数据点
        for i, (pos, data) in enumerate(zip(positions, data_to_plot)):
            x_jitter = np.random.normal(pos, 0.05, len(data))
            ax.scatter(x_jitter, data, alpha=0.5, s=30, color='black', zorder=3)
        
        # 设置x轴标签
        ax.set_xticks(positions)
        ax.set_xticklabels(['MVL', 'E'], fontsize=12, fontweight='bold')
        ax.set_ylabel(col, fontsize=11, fontweight='bold')
        p_val_str = f'{p_val:.4f}' if p_val and not np.isnan(p_val) else 'N/A'
        ax.set_title(f'{col}\n{stars} (p={p_val_str})', 
                    fontsize=12, fontweight='bold', pad=10)
        
        # 均值标注已删除（根据用户要求）
        
        # 添加显著性标记线
        if p_val and p_val < 0.05:
            y_max = max(mvl_vals.max(), e_vals.max())
            y_min = min(mvl_vals.min(), e_vals.min())
            y_range = y_max - y_min
            y_pos = y_max + y_range * 0.1
            
            ax.plot([1, 2], [y_pos, y_pos], 'k-', linewidth=1.5)
            ax.text(1.5, y_pos + y_range * 0.02, stars, 
                   ha='center', fontsize=14, fontweight='bold')
        
        ax.grid(True, alpha=0.3, axis='y')
    
    # 隐藏多余的子图
    for idx in range(n_params, len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_all_parameters.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_all_parameters.pdf", bbox_inches='tight')
    print(f"所有参数对比图已保存: {output_prefix}_all_parameters.png/pdf")
    
    # 创建关键参数的详细对比图
    key_params = ['Area', 'Aspect_Ratio', 'Eccentricity', 'Direction_Cos', 'Direction_Sin']
    available_key_params = [p for p in key_params if p in numeric_cols]
    
    if len(available_key_params) > 0:
        fig, axes = plt.subplots(1, len(available_key_params), figsize=(6*len(available_key_params), 6))
        if len(available_key_params) == 1:
            axes = [axes]
        
        for idx, col in enumerate(available_key_params):
            ax = axes[idx]
            
            mvl_vals = df_mvl[col].dropna()
            e_vals = df_e[col].dropna()
            stat_row = stats_df[stats_df['Parameter'] == col].iloc[0]
            p_val = stat_row['P_value']
            stars = stat_row['Significance']
            
            # 箱线图
            bp = ax.boxplot([mvl_vals.values, e_vals.values], 
                           labels=['MVL', 'E'],
                           patch_artist=True,
                           widths=0.6,
                           showmeans=True)
            
            # 设置颜色（使用用户指定的颜色）
            # MVL: #4d7e54, E: #2c4ca0
            bp['boxes'][0].set_facecolor('#4d7e54')
            bp['boxes'][0].set_alpha(0.7)
            bp['boxes'][1].set_facecolor('#2c4ca0')
            bp['boxes'][1].set_alpha(0.7)
            
            # 添加散点
            x1 = np.random.normal(1, 0.05, len(mvl_vals))
            x2 = np.random.normal(2, 0.05, len(e_vals))
            ax.scatter(x1, mvl_vals.values, alpha=0.6, s=50, color='#4d7e54', zorder=3)
            ax.scatter(x2, e_vals.values, alpha=0.6, s=50, color='#2c4ca0', zorder=3)
            
            ax.set_ylabel(col, fontsize=13, fontweight='bold')
            p_val_str = f'{p_val:.4f}' if p_val and not np.isnan(p_val) else 'N/A'
            ax.set_title(f'{col}\n{stars} (p={p_val_str})', 
                        fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
            
            # 添加统计信息
            textstr = f'MVL: n={len(mvl_vals)}, μ={mvl_vals.mean():.3f}±{mvl_vals.std():.3f}\n'
            textstr += f'E: n={len(e_vals)}, μ={e_vals.mean():.3f}±{e_vals.std():.3f}'
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, 
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_key_parameters.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_key_parameters.pdf", bbox_inches='tight')
        print(f"关键参数对比图已保存: {output_prefix}_key_parameters.png/pdf")
    
    # 打印统计摘要
    print("\n" + "="*80)
    print("统计学分析摘要")
    print("="*80)
    print(stats_df.to_string(index=False))
    print("="*80)
    print("\n显著性标记说明: *** p<0.001, ** p<0.01, * p<0.05, ns: 不显著")
    
    return stats_df

def main():
    """主函数"""
    print("="*80)
    print("MVL和E核团分析结果可视化")
    print("="*80)
    
    # --- 路径设置 ---
    # 脚本将自动查找数据文件和设置输出目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    excel_file = os.path.join(script_dir, "RF_MVL_E.xlsx")

    print(f"\n读取Excel文件: {excel_file}")
    
    try:
        # 读取分析结果sheet
        df_mvl = pd.read_excel(excel_file, sheet_name='MVL_分析结果')
        df_e = pd.read_excel(excel_file, sheet_name='E_分析结果')
        
        print(f"MVL数据: {df_mvl.shape[0]} 个样本, {df_mvl.shape[1]} 个参数")
        print(f"E数据: {df_e.shape[0]} 个样本, {df_e.shape[1]} 个参数")
        
        # --- 创建可视化 ---
        # 定义输出文件的前缀，并确保它们保存在正确的输出目录中
        output_prefix = os.path.join(script_dir, "mvl_e_comparison")
        stats_df = create_comparison_plots(df_mvl, df_e, output_prefix=output_prefix)
        
        print("\n可视化完成！")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

