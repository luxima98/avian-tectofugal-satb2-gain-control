import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# 获取脚本所在目录
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except:
    script_dir = os.getcwd()
os.chdir(script_dir)
print(f"工作目录: {os.getcwd()}")

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

print("=" * 60)
print("Spine Density Analysis Script")
print("=" * 60)

# ==================== 第一部分：数据读取和分类 ====================
print("\n[步骤1] 正在读取数据...")
e_mvl = pd.read_csv('E_MVL/spine_analysis_detail.csv')
mvl1 = pd.read_csv('MVL/spine_data_detail.csv')
mvl2 = pd.read_csv('MVL/dendrite_spine_details_20251111_200858.csv')
mvl3 = pd.read_csv('MVL/detail.csv')

# 统一列名
e_mvl.columns = ['File', 'Dendrite_Length', 'Dendrite_Width', 'Spine_Number', 
                 'Shaft_Length', 'Shaft_Width', 'Head_Length', 'Head_Width', 'Spine_Length']
mvl1.columns = ['File', 'Dendrite_Length', 'Dendrite_Width', 'Spine_Number',
                'Shaft_Length', 'Shaft_Width', 'Head_Length', 'Head_Width', 'Spine_Length']
mvl2.columns = ['File', 'Dendrite_Length', 'Dendrite_Width', 'Spine_Number',
                'Shaft_Length', 'Shaft_Width', 'Head_Length', 'Head_Width', 'Spine_Length']
mvl3.columns = ['File', 'Dendrite_Length', 'Dendrite_Width', 'Spine_Number',
                'Shaft_Length', 'Shaft_Width', 'Head_Length', 'Head_Width', 'Spine_Length']

# 添加组别标识
e_mvl['Group'] = 'E_MVL'
mvl1['Group'] = 'MVL'
mvl2['Group'] = 'MVL'
mvl3['Group'] = 'MVL'

# 合并所有数据
all_data = pd.concat([e_mvl, mvl1, mvl2, mvl3], ignore_index=True)
print(f"总共有 {len(all_data)} 个spine数据")
print(f"E_MVL: {len(e_mvl)} 个spine")
print(f"MVL: {len(mvl1) + len(mvl2) + len(mvl3)} 个spine")

# 计算关键比值
all_data['Head_Shaft_Ratio'] = all_data['Head_Width'] / (all_data['Shaft_Width'] + 1e-6)

# 基于规则的初步分类
def classify_by_rules(row):
    hw_sw_ratio = row['Head_Shaft_Ratio']
    shaft_len = row['Shaft_Length']
    
    # Stubby: shaft length几乎为0 且 head width/shaft width ≈ 1
    if shaft_len < 0.15 and 0.8 <= hw_sw_ratio <= 1.2:
        return 'Stubby'
    
    # Filopodia: head width/shaft width ≈ 1
    elif 0.8 <= hw_sw_ratio <= 1.2:
        return 'Filopodia'
    
    # Mushroom 和 Thin: head width/shaft width > 1
    elif hw_sw_ratio > 1.2:
        if shaft_len < 0.5:
            return 'Mushroom_Candidate'
        else:
            return 'Thin_Candidate'
    
    else:
        return 'Unclassified'

all_data['Rule_Classification'] = all_data.apply(classify_by_rules, axis=1)

# 对Mushroom和Thin候选进行统计分类
mushroom_thin_candidates = all_data[all_data['Rule_Classification'].isin(['Mushroom_Candidate', 'Thin_Candidate'])].copy()

if len(mushroom_thin_candidates) > 0:
    features = ['Head_Shaft_Ratio', 'Shaft_Length', 'Head_Width', 'Shaft_Width', 'Head_Length']
    X = mushroom_thin_candidates[features].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)
    
    cluster_0_mean_ratio = mushroom_thin_candidates.iloc[clusters == 0]['Head_Shaft_Ratio'].mean()
    cluster_0_mean_shaft = mushroom_thin_candidates.iloc[clusters == 0]['Shaft_Length'].mean()
    cluster_1_mean_ratio = mushroom_thin_candidates.iloc[clusters == 1]['Head_Shaft_Ratio'].mean()
    cluster_1_mean_shaft = mushroom_thin_candidates.iloc[clusters == 1]['Shaft_Length'].mean()
    
    if cluster_0_mean_ratio > cluster_1_mean_ratio:
        cluster_labels = {0: 'Mushroom', 1: 'Thin'}
    else:
        cluster_labels = {1: 'Mushroom', 0: 'Thin'}
    
    mushroom_thin_candidates['Final_Classification'] = [cluster_labels[c] for c in clusters]
    all_data.loc[mushroom_thin_candidates.index, 'Final_Classification'] = mushroom_thin_candidates['Final_Classification']

# 更新最终分类
all_data['Final_Classification'] = all_data['Final_Classification'].fillna(all_data['Rule_Classification'])
all_data['Final_Classification'] = all_data['Final_Classification'].replace({
    'Mushroom_Candidate': 'Mushroom',
    'Thin_Candidate': 'Thin',
    'Unclassified': 'Other'
})

print("\n[步骤2] 分类完成")
print("分类结果统计:")
print(all_data['Final_Classification'].value_counts())

# ==================== 第二部分：计算Density ====================
print("\n[步骤3] 正在计算Spine Density...")

# 按File（dendrite）分组计算
def calculate_density_by_dendrite(group_data):
    """计算每条dendrite的spine density"""
    dendrite_length = group_data['Dendrite_Length'].iloc[0]  # 每条dendrite的长度相同
    
    if dendrite_length <= 0:
        return pd.Series({
            'Dendrite_Length': dendrite_length,
            'Total_Spine_Count': len(group_data),
            'Total_Spine_Density': 0,
            'Mushroom_Count': 0,
            'Mushroom_Density': 0,
            'Thin_Count': 0,
            'Thin_Density': 0,
            'Filopodia_Count': 0,
            'Filopodia_Density': 0,
            'Stubby_Count': 0,
            'Stubby_Density': 0,
            'Mushroom_Proportion': 0,
            'Thin_Proportion': 0,
            'Filopodia_Proportion': 0,
            'Stubby_Proportion': 0
        })
    
    total_count = len(group_data)
    total_density = total_count / dendrite_length  # number/μm
    
    # 计算各类型数量
    mushroom_count = len(group_data[group_data['Final_Classification'] == 'Mushroom'])
    thin_count = len(group_data[group_data['Final_Classification'] == 'Thin'])
    filopodia_count = len(group_data[group_data['Final_Classification'] == 'Filopodia'])
    stubby_count = len(group_data[group_data['Final_Classification'] == 'Stubby'])
    
    # 计算各类型density
    mushroom_density = mushroom_count / dendrite_length
    thin_density = thin_count / dendrite_length
    filopodia_density = filopodia_count / dendrite_length
    stubby_density = stubby_count / dendrite_length
    
    # 计算占比
    mushroom_prop = mushroom_count / total_count if total_count > 0 else 0
    thin_prop = thin_count / total_count if total_count > 0 else 0
    filopodia_prop = filopodia_count / total_count if total_count > 0 else 0
    stubby_prop = stubby_count / total_count if total_count > 0 else 0
    
    return pd.Series({
        'Dendrite_Length': dendrite_length,
        'Total_Spine_Count': total_count,
        'Total_Spine_Density': total_density,
        'Mushroom_Count': mushroom_count,
        'Mushroom_Density': mushroom_density,
        'Thin_Count': thin_count,
        'Thin_Density': thin_density,
        'Filopodia_Count': filopodia_count,
        'Filopodia_Density': filopodia_density,
        'Stubby_Count': stubby_count,
        'Stubby_Density': stubby_density,
        'Mushroom_Proportion': mushroom_prop,
        'Thin_Proportion': thin_prop,
        'Filopodia_Proportion': filopodia_prop,
        'Stubby_Proportion': stubby_prop
    })

# 按File分组计算
density_summary = all_data.groupby(['File', 'Group']).apply(calculate_density_by_dendrite).reset_index()
# 移除多余的列（如果有）
if 'level_2' in density_summary.columns:
    density_summary = density_summary.drop('level_2', axis=1)

# 分离E_MVL和MVL数据
e_mvl_density = density_summary[density_summary['Group'] == 'E_MVL'].copy()
mvl_density = density_summary[density_summary['Group'] == 'MVL'].copy()

print(f"\nE_MVL: {len(e_mvl_density)} 条dendrite")
print(f"MVL: {len(mvl_density)} 条dendrite")

# 保存数据表
e_mvl_density.to_csv('E_MVL_density_analysis.csv', index=False)
mvl_density.to_csv('MVL_density_analysis.csv', index=False)
print("\n[步骤4] Density数据表已保存:")
print("  - E_MVL_density_analysis.csv")
print("  - MVL_density_analysis.csv")

# ==================== 第三部分：统计学分析 ====================
print("\n[步骤5] 正在进行统计学分析...")

# 准备分析数据
metrics_to_analyze = {
    'Total_Spine_Density': 'Total Spine Density (number/μm)',
    'Mushroom_Density': 'Mushroom Density (number/μm)',
    'Thin_Density': 'Thin Density (number/μm)',
    'Filopodia_Density': 'Filopodia Density (number/μm)',
    'Stubby_Density': 'Stubby Density (number/μm)',
    'Mushroom_Proportion': 'Mushroom Proportion',
    'Thin_Proportion': 'Thin Proportion',
    'Filopodia_Proportion': 'Filopodia Proportion',
    'Stubby_Proportion': 'Stubby Proportion'
}

# 进行统计检验
statistical_results = []

for metric, metric_name in metrics_to_analyze.items():
    e_mvl_values = e_mvl_density[metric].values
    mvl_values = mvl_density[metric].values
    
    # 移除NaN和无效值
    e_mvl_values = e_mvl_values[~np.isnan(e_mvl_values)]
    mvl_values = mvl_values[~np.isnan(mvl_values)]
    
    if len(e_mvl_values) == 0 or len(mvl_values) == 0:
        continue
    
    # 描述性统计
    e_mvl_mean = np.mean(e_mvl_values)
    e_mvl_std = np.std(e_mvl_values, ddof=1)
    e_mvl_sem = stats.sem(e_mvl_values)
    
    mvl_mean = np.mean(mvl_values)
    mvl_std = np.std(mvl_values, ddof=1)
    mvl_sem = stats.sem(mvl_values)
    
    # 正态性检验
    _, e_mvl_normality = stats.shapiro(e_mvl_values) if len(e_mvl_values) <= 5000 else (None, 1.0)
    _, mvl_normality = stats.shapiro(mvl_values) if len(mvl_values) <= 5000 else (None, 1.0)
    
    # 选择统计检验方法
    if e_mvl_normality > 0.05 and mvl_normality > 0.05:
        # 正态分布，使用t检验
        statistic, p_value = stats.ttest_ind(e_mvl_values, mvl_values)
        test_name = 'Independent t-test'
    else:
        # 非正态分布，使用Mann-Whitney U检验
        statistic, p_value = stats.mannwhitneyu(e_mvl_values, mvl_values, alternative='two-sided')
        test_name = 'Mann-Whitney U test'
    
    # 效应量计算 (Cohen's d)
    pooled_std = np.sqrt(((len(e_mvl_values) - 1) * e_mvl_std**2 + (len(mvl_values) - 1) * mvl_std**2) / 
                        (len(e_mvl_values) + len(mvl_values) - 2))
    cohens_d = (e_mvl_mean - mvl_mean) / pooled_std if pooled_std > 0 else 0
    
    statistical_results.append({
        'Metric': metric_name,
        'E_MVL_Mean': e_mvl_mean,
        'E_MVL_STD': e_mvl_std,
        'E_MVL_SEM': e_mvl_sem,
        'MVL_Mean': mvl_mean,
        'MVL_STD': mvl_std,
        'MVL_SEM': mvl_sem,
        'Test': test_name,
        'Statistic': statistic,
        'P_Value': p_value,
        'Significant': 'Yes' if p_value < 0.05 else 'No',
        'Cohens_d': cohens_d
    })
    
    print(f"\n{metric_name}:")
    print(f"  E_MVL: {e_mvl_mean:.4f} ± {e_mvl_sem:.4f} (SEM)")
    print(f"  MVL: {mvl_mean:.4f} ± {mvl_sem:.4f} (SEM)")
    print(f"  {test_name}: statistic={statistic:.4f}, p={p_value:.4f}")
    print(f"  {'***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'}")

# 对配比进行统计检验
print("\n[步骤6.2] 对Spine类型配比进行统计检验...")
proportion_metrics = {
    'Mushroom_Proportion': 'Mushroom Proportion',
    'Thin_Proportion': 'Thin Proportion',
    'Filopodia_Proportion': 'Filopodia Proportion',
    'Stubby_Proportion': 'Stubby Proportion'
}

for metric, metric_name in proportion_metrics.items():
    e_mvl_values = e_mvl_density[metric].dropna() * 100  # 转换为百分比
    mvl_values = mvl_density[metric].dropna() * 100
    
    if len(e_mvl_values) == 0 or len(mvl_values) == 0:
        continue
    
    e_mvl_mean = np.mean(e_mvl_values)
    e_mvl_std = np.std(e_mvl_values, ddof=1)
    e_mvl_sem = stats.sem(e_mvl_values)
    mvl_mean = np.mean(mvl_values)
    mvl_std = np.std(mvl_values, ddof=1)
    mvl_sem = stats.sem(mvl_values)
    
    # 正态性检验
    _, p_norm_e = stats.shapiro(e_mvl_values) if len(e_mvl_values) <= 5000 else (None, 0.01)
    _, p_norm_m = stats.shapiro(mvl_values) if len(mvl_values) <= 5000 else (None, 0.01)
    
    # 选择统计检验方法
    if p_norm_e > 0.05 and p_norm_m > 0.05:
        statistic, p_value = stats.ttest_ind(e_mvl_values, mvl_values)
        test_name = 't-test'
    else:
        statistic, p_value = stats.mannwhitneyu(e_mvl_values, mvl_values, alternative='two-sided')
        test_name = 'Mann-Whitney U test'
    
    # 效应量计算
    pooled_std = np.sqrt(((len(e_mvl_values) - 1) * e_mvl_std**2 + (len(mvl_values) - 1) * mvl_std**2) / 
                        (len(e_mvl_values) + len(mvl_values) - 2))
    cohens_d = (e_mvl_mean - mvl_mean) / pooled_std if pooled_std > 0 else 0
    
    statistical_results.append({
        'Metric': metric_name,
        'E_MVL_Mean': e_mvl_mean,
        'E_MVL_STD': e_mvl_std,
        'E_MVL_SEM': e_mvl_sem,
        'MVL_Mean': mvl_mean,
        'MVL_STD': mvl_std,
        'MVL_SEM': mvl_sem,
        'Test': test_name,
        'Statistic': statistic,
        'P_Value': p_value,
        'Significant': 'Yes' if p_value < 0.05 else 'No',
        'Cohens_d': cohens_d
    })
    
    print(f"\n{metric_name}:")
    print(f"  E_MVL: {e_mvl_mean:.2f}% ± {e_mvl_sem:.2f}% (SEM)")
    print(f"  MVL: {mvl_mean:.2f}% ± {mvl_sem:.2f}% (SEM)")
    print(f"  {test_name}: statistic={statistic:.4f}, p={p_value:.4f}")
    print(f"  {'***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'}")

# 保存统计结果
stats_df = pd.DataFrame(statistical_results)
stats_df.to_csv('statistical_analysis_results.csv', index=False)
print("\n[步骤6] 统计学分析结果已保存到 statistical_analysis_results.csv")

# ==================== 配比分析 ====================
print("\n[步骤6.5] 正在进行配比分析...")

# 计算每个组的主要spine类型
def get_dominant_spine_type(row):
    """确定每条dendrite的主要spine类型"""
    props = {
        'Mushroom': row['Mushroom_Proportion'],
        'Thin': row['Thin_Proportion'],
        'Filopodia': row['Filopodia_Proportion'],
        'Stubby': row['Stubby_Proportion']
    }
    dominant = max(props, key=props.get)
    return dominant

# 添加主要类型列
e_mvl_density['Dominant_Spine_Type'] = e_mvl_density.apply(get_dominant_spine_type, axis=1)
mvl_density['Dominant_Spine_Type'] = mvl_density.apply(get_dominant_spine_type, axis=1)

# 统计每个组的主要spine类型分布
print("\nE_MVL组的主要Spine类型分布:")
e_mvl_dominant = e_mvl_density['Dominant_Spine_Type'].value_counts()
print(e_mvl_dominant)
e_mvl_dominant_pct = e_mvl_density['Dominant_Spine_Type'].value_counts(normalize=True) * 100
print("\nE_MVL组的主要Spine类型占比:")
for spine_type, pct in e_mvl_dominant_pct.items():
    print(f"  {spine_type}: {pct:.2f}%")

print("\nMVL组的主要Spine类型分布:")
mvl_dominant = mvl_density['Dominant_Spine_Type'].value_counts()
print(mvl_dominant)
mvl_dominant_pct = mvl_density['Dominant_Spine_Type'].value_counts(normalize=True) * 100
print("\nMVL组的主要Spine类型占比:")
for spine_type, pct in mvl_dominant_pct.items():
    print(f"  {spine_type}: {pct:.2f}%")

# 确定每个组的总体主要类型
e_mvl_overall_dominant = e_mvl_dominant.index[0] if len(e_mvl_dominant) > 0 else 'Unknown'
mvl_overall_dominant = mvl_dominant.index[0] if len(mvl_dominant) > 0 else 'Unknown'

print(f"\nE_MVL组的主要Spine类型: {e_mvl_overall_dominant} ({e_mvl_dominant_pct[e_mvl_overall_dominant]:.2f}%)")
print(f"MVL组的主要Spine类型: {mvl_overall_dominant} ({mvl_dominant_pct[mvl_overall_dominant]:.2f}%)")

# 计算平均配比
print("\n各组的平均Spine类型配比:")
print("E_MVL组:")
print(f"  Mushroom: {e_mvl_density['Mushroom_Proportion'].mean()*100:.2f}%")
print(f"  Thin: {e_mvl_density['Thin_Proportion'].mean()*100:.2f}%")
print(f"  Filopodia: {e_mvl_density['Filopodia_Proportion'].mean()*100:.2f}%")
print(f"  Stubby: {e_mvl_density['Stubby_Proportion'].mean()*100:.2f}%")

print("\nMVL组:")
print(f"  Mushroom: {mvl_density['Mushroom_Proportion'].mean()*100:.2f}%")
print(f"  Thin: {mvl_density['Thin_Proportion'].mean()*100:.2f}%")
print(f"  Filopodia: {mvl_density['Filopodia_Proportion'].mean()*100:.2f}%")
print(f"  Stubby: {mvl_density['Stubby_Proportion'].mean()*100:.2f}%")

# 保存更新后的数据表（包含主要类型）
e_mvl_density.to_csv('E_MVL_density_analysis.csv', index=False)
mvl_density.to_csv('MVL_density_analysis.csv', index=False)

# ==================== 第四部分：可视化 ====================
print("\n[步骤7] 正在生成可视化图表...")

# 创建可视化
fig, axes = plt.subplots(3, 3, figsize=(20, 18))
fig.suptitle('Spine Density Analysis: E_MVL vs MVL', fontsize=20, fontweight='bold', y=0.995)

# 准备数据用于绘图
plot_data = []
for metric, metric_name in metrics_to_analyze.items():
    e_mvl_vals = e_mvl_density[metric].dropna()
    mvl_vals = mvl_density[metric].dropna()
    
    for val in e_mvl_vals:
        plot_data.append({'Group': 'E_MVL', 'Metric': metric_name, 'Value': val})
    for val in mvl_vals:
        plot_data.append({'Group': 'MVL', 'Metric': metric_name, 'Value': val})

plot_df = pd.DataFrame(plot_data)

# 绘制每个指标的箱线图和统计检验
group_palette = {'E_MVL': '#2c4ca0', 'MVL': '#4d7e54'}
print(f"\n[调色信息] 箱线图颜色: E_MVL={group_palette['E_MVL']} | MVL={group_palette['MVL']}")
metrics_keys = list(metrics_to_analyze.keys())
for idx, metric_key in enumerate(metrics_keys):
    metric_name = metrics_to_analyze[metric_key]
    row, col = divmod(idx, 3)
    ax = axes[row, col]
    
    data_e = e_mvl_density[metric_key].dropna()
    data_m = mvl_density[metric_key].dropna()
    
    bp = ax.boxplot([data_e, data_m], labels=['E_MVL', 'MVL'],
                    patch_artist=True, widths=0.6)
    
    box_colors = [group_palette['E_MVL'], group_palette['MVL']]
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)
        patch.set_alpha(1.0)
    for median, color in zip(bp['medians'], box_colors):
        median.set_color(color)
        median.set_linewidth(2)
    whisker_colors = box_colors[0:1]*2 + box_colors[1:2]*2
    for element in ['whiskers', 'caps', 'fliers']:
        for artist, color in zip(bp[element], whisker_colors):
            if element == 'fliers':
                artist.set(markerfacecolor=color, markeredgecolor=color, alpha=0.8)
            else:
                artist.set_color(color)
                artist.set_linewidth(1.5)
    
    x_e = np.random.normal(1, 0.04, len(data_e))
    x_m = np.random.normal(2, 0.04, len(data_m))
    ax.scatter(x_e, data_e, alpha=0.5, s=20, color=group_palette['E_MVL'])
    ax.scatter(x_m, data_m, alpha=0.5, s=20, color=group_palette['MVL'])
    
    mean_e, sem_e = np.mean(data_e), stats.sem(data_e) if len(data_e) > 1 else 0
    mean_m, sem_m = np.mean(data_m), stats.sem(data_m) if len(data_m) > 1 else 0
    
    ax.errorbar([1, 2], [mean_e, mean_m],
               yerr=[sem_e, sem_m],
               fmt='o', color='black', markersize=8,
               capsize=5, capthick=2, label='Mean ± SEM')
    
    result = stats_df[stats_df['Metric'] == metric_name].iloc[0]
    p_val = result['P_Value']
    
    y_max = max(data_e.max(), data_m.max())
    y_min = min(data_e.min(), data_m.min())
    y_range = max(y_max - y_min, 1e-6)
    
    sig_text = '***' if p_val < 0.001 else \
               '**' if p_val < 0.01 else \
               '*' if p_val < 0.05 else 'ns'
    
    line_y = y_max + y_range * 0.1
    sig_y = line_y + y_range * 0.05
    p_y = sig_y + y_range * 0.05
    ax.plot([1, 2], [line_y, line_y], 'k-', linewidth=1)
    ax.text(1.5, sig_y, sig_text, ha='center', fontsize=12, fontweight='bold')
    ax.text(1.5, p_y, f'p={p_val:.3f}', ha='center', fontsize=9)
    ax.set_ylim(top=max(ax.get_ylim()[1], y_max + y_range * 0.3))
    
    ax.set_ylabel(metric_name, fontsize=11)
    ax.set_title(metric_name, fontsize=12, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, axis='y')
    if idx == 0:
        ax.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig('spine_density_statistical_analysis.png', dpi=300, bbox_inches='tight')
print("可视化图表已保存到 spine_density_statistical_analysis.png")

# ==================== 配比可视化 ====================
print("\n[步骤7.5] 正在生成配比可视化图表...")

# 创建配比分析图
fig2, axes2 = plt.subplots(2, 2, figsize=(16, 14))
fig2.suptitle('Spine Type Proportion Analysis: E_MVL vs MVL', fontsize=18, fontweight='bold', y=0.995)

# 1. 堆叠柱状图 - 平均配比
ax1 = axes2[0, 0]
spine_types_ordered = ['Mushroom', 'Thin', 'Filopodia', 'Stubby']
colors_prop = {
    'Mushroom': '#326db6',
    'Thin': '#d16d5b',
    'Filopodia': '#669877',
    'Stubby': '#b25f79'
}

e_mvl_props = [
    e_mvl_density['Mushroom_Proportion'].mean() * 100,
    e_mvl_density['Thin_Proportion'].mean() * 100,
    e_mvl_density['Filopodia_Proportion'].mean() * 100,
    e_mvl_density['Stubby_Proportion'].mean() * 100
]

mvl_props = [
    mvl_density['Mushroom_Proportion'].mean() * 100,
    mvl_density['Thin_Proportion'].mean() * 100,
    mvl_density['Filopodia_Proportion'].mean() * 100,
    mvl_density['Stubby_Proportion'].mean() * 100
]

x = np.arange(2)
width = 0.6
bottom_e = 0
bottom_m = 0

bars = []
for i, spine_type in enumerate(spine_types_ordered):
    bar1 = ax1.bar(0, e_mvl_props[i], width, bottom=bottom_e, 
                   color=colors_prop[spine_type], alpha=0.8)
    ax1.bar(1, mvl_props[i], width, bottom=bottom_m, 
           color=colors_prop[spine_type], alpha=0.8)
    bars.append(bar1)
    
    # 添加数值标签
    if e_mvl_props[i] > 2:
        ax1.text(0, bottom_e + e_mvl_props[i]/2, f'{e_mvl_props[i]:.1f}%', 
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    if mvl_props[i] > 2:
        ax1.text(1, bottom_m + mvl_props[i]/2, f'{mvl_props[i]:.1f}%', 
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    bottom_e += e_mvl_props[i]
    bottom_m += mvl_props[i]

ax1.set_ylabel('Proportion (%)', fontsize=12)
ax1.set_title('Average Spine Type Proportion', fontsize=14, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(['E_MVL', 'MVL'])
ax1.set_ylim(0, 100)
ax1.legend(bars, spine_types_ordered, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')

# 2. 主要类型分布 - 饼图
ax2 = axes2[0, 1]
e_mvl_dominant_counts = e_mvl_density['Dominant_Spine_Type'].value_counts()
e_mvl_dominant_counts = e_mvl_dominant_counts.reindex(spine_types_ordered, fill_value=0)
colors_pie = [colors_prop[st] for st in e_mvl_dominant_counts.index]
wedges, texts, autotexts = ax2.pie(e_mvl_dominant_counts.values, 
                                   labels=e_mvl_dominant_counts.index,
                                   autopct='%1.1f%%',
                                   colors=colors_pie,
                                   startangle=90)
ax2.set_title('E_MVL: Dominant Spine Type Distribution', fontsize=14, fontweight='bold')

# 3. MVL主要类型分布 - 饼图
ax3 = axes2[1, 0]
mvl_dominant_counts = mvl_density['Dominant_Spine_Type'].value_counts()
mvl_dominant_counts = mvl_dominant_counts.reindex(spine_types_ordered, fill_value=0)
colors_pie_m = [colors_prop[st] for st in mvl_dominant_counts.index]
wedges, texts, autotexts = ax3.pie(mvl_dominant_counts.values, 
                                   labels=mvl_dominant_counts.index,
                                   autopct='%1.1f%%',
                                   colors=colors_pie_m,
                                   startangle=90)
ax3.set_title('MVL: Dominant Spine Type Distribution', fontsize=14, fontweight='bold')

# 4. 配比对比箱线图
ax4 = axes2[1, 1]
prop_data_to_plot = []
prop_labels = []
for spine_type in spine_types_ordered:
    prop_data_to_plot.append(e_mvl_density[f'{spine_type}_Proportion'].dropna() * 100)
    prop_data_to_plot.append(mvl_density[f'{spine_type}_Proportion'].dropna() * 100)
    prop_labels.extend([f'E_MVL\n{spine_type}', f'MVL\n{spine_type}'])

bp = ax4.boxplot(prop_data_to_plot, labels=prop_labels, patch_artist=True, widths=0.6)
for i, patch in enumerate(bp['boxes']):
    spine_idx = i // 2
    group_idx = i % 2
    patch.set_facecolor(colors_prop[spine_types_ordered[spine_idx]])
    patch.set_alpha(0.7 if group_idx == 0 else 0.5)

ax4.set_ylabel('Proportion (%)', fontsize=12)
ax4.set_title('Spine Type Proportion Comparison', fontsize=14, fontweight='bold')
ax4.tick_params(axis='x', rotation=45)
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('spine_proportion_analysis.png', dpi=300, bbox_inches='tight')
print("配比可视化图表已保存到 spine_proportion_analysis.png")

# ==================== 第五部分：生成汇总报告 ====================
print("\n[步骤8] 正在生成汇总报告...")

report = f"""
================================================================================
                    SPINE DENSITY ANALYSIS REPORT
================================================================================

1. DATA SUMMARY
   - Total spines analyzed: {len(all_data)}
   - E_MVL dendrites: {len(e_mvl_density)}
   - MVL dendrites: {len(mvl_density)}

2. CLASSIFICATION RESULTS
{all_data['Final_Classification'].value_counts().to_string()}

3. STATISTICAL COMPARISON (E_MVL vs MVL)
"""

for result in statistical_results:
    report += f"""
   {result['Metric']}:
   - E_MVL: {result['E_MVL_Mean']:.4f} ± {result['E_MVL_SEM']:.4f} (SEM)
   - MVL: {result['MVL_Mean']:.4f} ± {result['MVL_SEM']:.4f} (SEM)
   - {result['Test']}: p = {result['P_Value']:.4f} {'***' if result['P_Value'] < 0.001 else '**' if result['P_Value'] < 0.01 else '*' if result['P_Value'] < 0.05 else '(ns)'}
   - Cohen's d = {result['Cohens_d']:.4f}

"""

report += f"""
4. PROPORTION ANALYSIS

   E_MVL Group:
   - Average Mushroom Proportion: {e_mvl_density['Mushroom_Proportion'].mean()*100:.2f}%
   - Average Thin Proportion: {e_mvl_density['Thin_Proportion'].mean()*100:.2f}%
   - Average Filopodia Proportion: {e_mvl_density['Filopodia_Proportion'].mean()*100:.2f}%
   - Average Stubby Proportion: {e_mvl_density['Stubby_Proportion'].mean()*100:.2f}%
   
   MVL Group:
   - Average Mushroom Proportion: {mvl_density['Mushroom_Proportion'].mean()*100:.2f}%
   - Average Thin Proportion: {mvl_density['Thin_Proportion'].mean()*100:.2f}%
   - Average Filopodia Proportion: {mvl_density['Filopodia_Proportion'].mean()*100:.2f}%
   - Average Stubby Proportion: {mvl_density['Stubby_Proportion'].mean()*100:.2f}%

5. DOMINANT SPINE TYPE ANALYSIS

   E_MVL Group:
   - Dominant Spine Type: {e_mvl_overall_dominant}
   - Distribution of Dominant Types:
"""
for spine_type, count in e_mvl_dominant.items():
    pct = e_mvl_dominant_pct[spine_type]
    report += f"     * {spine_type}: {count} dendrites ({pct:.2f}%)\n"

report += f"""
   MVL Group:
   - Dominant Spine Type: {mvl_overall_dominant}
   - Distribution of Dominant Types:
"""
for spine_type, count in mvl_dominant.items():
    pct = mvl_dominant_pct[spine_type]
    report += f"     * {spine_type}: {count} dendrites ({pct:.2f}%)\n"

report += """
6. FILES GENERATED
   - E_MVL_density_analysis.csv: Density data for E_MVL group
   - MVL_density_analysis.csv: Density data for MVL group
   - statistical_analysis_results.csv: Statistical test results
   - spine_density_statistical_analysis.png: Visualization of comparisons
   - spine_proportion_analysis.png: Proportion analysis visualization

================================================================================
"""

print(report)

# 保存报告
with open('spine_density_analysis_report.txt', 'w', encoding='utf-8') as f:
    f.write(report)

print("\n" + "=" * 60)
print("所有分析完成！")
print("=" * 60)
print("\n生成的文件:")
print("  1. E_MVL_density_analysis.csv - E_MVL组density数据表（包含配比和主要类型）")
print("  2. MVL_density_analysis.csv - MVL组density数据表（包含配比和主要类型）")
print("  3. statistical_analysis_results.csv - 统计学分析结果（包含配比统计）")
print("  4. spine_density_statistical_analysis.png - Density对比可视化图表")
print("  5. spine_proportion_analysis.png - 配比分析可视化图表")
print("  6. spine_density_analysis_report.txt - 完整分析报告")
print("\n主要发现:")
print(f"  - E_MVL组的主要Spine类型: {e_mvl_overall_dominant}")
print(f"  - MVL组的主要Spine类型: {mvl_overall_dominant}")

