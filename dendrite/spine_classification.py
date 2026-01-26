import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.patches as patches
from matplotlib.patches import Rectangle, FancyBboxPatch
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def identify_stubby_spines(df):
    """识别stubby类型的spine
    
    条件：
    1. head length < dendrite length
    2. shaft width 和 head width 大致相似（差异小于20%）
    3. 或者 shaft length 几乎为0（< 0.1）
    """
    df = df.copy()
    
    # 计算width相似度（差异百分比）
    width_diff_ratio = abs(df['shaft_width'] - df['head_width']) / (df['shaft_width'] + df['head_width'] + 1e-10)
    
    # 判断条件
    condition1 = df['head_length'] < df['dendrite_length']  # head length < dendrite length
    condition2 = width_diff_ratio < 0.2  # shaft width 和 head width 大致相似（差异<20%）
    condition3 = df['shaft_length'] < 0.1  # shaft length 几乎为0
    
    # stubby条件：满足(条件1 AND 条件2) OR 条件3
    is_stubby = (condition1 & condition2) | condition3
    
    # 对于stubby类型，length = head length, width = head width
    df.loc[is_stubby, 'spine_type'] = 'stubby'
    df.loc[is_stubby, 'spine_length_stubby'] = df.loc[is_stubby, 'head_length']
    df.loc[is_stubby, 'spine_width_stubby'] = df.loc[is_stubby, 'head_width']
    
    return df, is_stubby

def cluster_spines(df_non_stubby):
    """对非stubby的spine进行聚类分析
    
    使用特征：
    - shaft_width
    - head_width / shaft_width (ratio)
    - shaft_length
    - shaft_length + head_length (total length)
    """
    # 准备特征
    features = df_non_stubby.copy()
    
    # 计算 head_width / shaft_width 比率（避免除零）
    features['head_shaft_width_ratio'] = features['head_width'] / (features['shaft_width'] + 1e-10)
    
    # 选择用于聚类的特征
    X = features[['shaft_width', 'head_shaft_width_ratio', 'shaft_length', 'spine_length']].values
    
    # 标准化特征
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # K-means聚类，分成3类
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)
    
    # 将聚类结果和比率添加到数据框
    df_non_stubby = df_non_stubby.copy()
    df_non_stubby['cluster'] = clusters
    df_non_stubby['head_shaft_width_ratio'] = features['head_shaft_width_ratio'].values
    
    # 根据特征值给聚类命名（thin, mushroom, filopodia等）
    # 计算每个聚类的平均特征值
    cluster_stats = df_non_stubby.groupby('cluster').agg({
        'shaft_width': 'mean',
        'head_width': 'mean',
        'shaft_length': 'mean',
        'head_length': 'mean',
        'spine_length': 'mean',
        'head_shaft_width_ratio': 'mean'
    })
    
    # 根据特征命名（可以根据实际情况调整）
    # 通常：thin (细长), mushroom (蘑菇型), filopodia (丝状)
    cluster_names = {}
    
    # 按shaft_length排序，最长的可能是thin
    # 按head_width/shaft_width排序，比率最大的可能是mushroom
    # 按spine_length排序
    
    # 计算每个聚类的特征值用于排序
    cluster_features = []
    for cluster_id in range(3):
        stats = cluster_stats.loc[cluster_id]
        cluster_features.append({
            'id': cluster_id,
            'shaft_length': stats['shaft_length'],
            'head_shaft_ratio': stats['head_shaft_width_ratio'],
            'head_width': stats['head_width'],
            'shaft_width': stats['shaft_width'],
            'spine_length': stats['spine_length']
        })
    
    # 根据特征值命名
    # 1. 找出shaft_length最长的 -> thin
    # 2. 找出head_width/shaft_width比率最大的 -> mushroom  
    # 3. 剩下的 -> filopodia
    
    sorted_by_shaft_length = sorted(cluster_features, key=lambda x: x['shaft_length'], reverse=True)
    sorted_by_ratio = sorted(cluster_features, key=lambda x: x['head_shaft_ratio'], reverse=True)
    
    # 如果shaft_length差异明显，最长的为thin
    if sorted_by_shaft_length[0]['shaft_length'] > sorted_by_shaft_length[1]['shaft_length'] * 1.3:
        cluster_names[sorted_by_shaft_length[0]['id']] = 'thin'
        remaining = [c for c in cluster_features if c['id'] != sorted_by_shaft_length[0]['id']]
    else:
        remaining = cluster_features
        cluster_names[sorted_by_shaft_length[0]['id']] = 'filopodia'  # 默认
    
    # 在剩余的中，找出head_width/shaft_width比率最大的 -> mushroom
    if len(remaining) > 0:
        sorted_remaining = sorted(remaining, key=lambda x: x['head_shaft_ratio'], reverse=True)
        if sorted_remaining[0]['head_shaft_ratio'] > 1.3:
            cluster_names[sorted_remaining[0]['id']] = 'mushroom'
            final_remaining = [c for c in remaining if c['id'] != sorted_remaining[0]['id']]
        else:
            final_remaining = remaining
            cluster_names[sorted_remaining[0]['id']] = 'filopodia'
        
        # 最后一个
        for c in final_remaining:
            if c['id'] not in cluster_names:
                cluster_names[c['id']] = 'filopodia'
    
    # 确保所有聚类都有名称
    for cluster_id in range(3):
        if cluster_id not in cluster_names:
            cluster_names[cluster_id] = 'filopodia'
    
    # 添加类型标签
    df_non_stubby['spine_type'] = df_non_stubby['cluster'].map(cluster_names)
    
    return df_non_stubby, cluster_stats, cluster_names, X_scaled, scaler, kmeans

def draw_spine_shape(ax, spine_type, stats, x_offset=0, scale=1.0):
    """绘制单个spine的形状
    
    spine形状：shaft (细长部分) + head (头部)
    """
    if spine_type == 'stubby':
        # stubby类型：只有head，没有明显的shaft
        length = stats['head_length'] * scale
        width = stats['head_width'] * scale
        
        # 绘制一个椭圆或圆形（短粗型）
        ellipse = patches.Ellipse((x_offset, length/2), width, length, 
                                  angle=0, color='skyblue', alpha=0.7, edgecolor='black', linewidth=2)
        ax.add_patch(ellipse)
        
        # 添加dendrite连接线（底部）
        ax.plot([x_offset, x_offset], [0, length/4], 'k-', linewidth=2)
        
        # 添加标签
        ax.text(x_offset, -length*0.2, f'Stubby\nLength={stats["head_length"]:.3f}\nWidth={stats["head_width"]:.3f}', 
                ha='center', va='top', fontsize=9, fontweight='bold')
        
        # 设置坐标轴范围
        max_dim = max(length, width)
        ax.set_xlim(x_offset - max_dim*0.7, x_offset + max_dim*0.7)
        ax.set_ylim(-0.3, max_dim*1.2)
        
    else:
        # 其他类型：shaft + head
        shaft_length = stats['shaft_length'] * scale
        shaft_width = stats['shaft_width'] * scale
        head_length = stats['head_length'] * scale
        head_width = stats['head_width'] * scale
        
        # 绘制shaft（细长部分，从底部开始）
        shaft_rect = Rectangle((x_offset - shaft_width/2, 0), 
                              shaft_width, shaft_length,
                              facecolor='lightblue', edgecolor='black', linewidth=2)
        ax.add_patch(shaft_rect)
        
        # 绘制head（头部，在shaft上方）
        head_rect = Rectangle((x_offset - head_width/2, shaft_length), 
                             head_width, head_length,
                             facecolor='lightcoral', edgecolor='black', linewidth=2)
        ax.add_patch(head_rect)
        
        # 添加标签
        total_length = shaft_length + head_length
        label_text = f'{spine_type.capitalize()}\nShaft L={stats["shaft_length"]:.3f}\nShaft W={stats["shaft_width"]:.3f}\nHead L={stats["head_length"]:.3f}\nHead W={stats["head_width"]:.3f}'
        ax.text(x_offset, -total_length*0.15, label_text, 
                ha='center', va='top', fontsize=8, fontweight='bold')
        
        # 设置坐标轴范围
        max_width = max(shaft_width, head_width)
        ax.set_xlim(x_offset - max_width*0.7, x_offset + max_width*0.7)
        ax.set_ylim(-0.3, total_length*1.1)
    
    ax.set_aspect('equal')
    ax.axis('off')

def visualize_kmeans_clustering(df_non_stubby, X_scaled, scaler, kmeans, cluster_names):
    """可视化K-means聚类过程"""
    # 使用PCA降维到2D进行可视化
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    # 获取聚类标签
    clusters = df_non_stubby['cluster'].values
    spine_types = df_non_stubby['spine_type'].values
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图：PCA降维后的聚类结果
    ax1 = axes[0]
    colors_map = {'thin': 'green', 'mushroom': 'red', 'filopodia': 'blue'}
    markers_map = {'thin': 'o', 'mushroom': 's', 'filopodia': '^'}
    
    for spine_type in ['thin', 'mushroom', 'filopodia']:
        mask = spine_types == spine_type
        if mask.sum() > 0:
            ax1.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                       c=colors_map[spine_type], marker=markers_map[spine_type],
                       label=spine_type.capitalize(), s=50, alpha=0.6, edgecolors='black', linewidth=0.5)
    
    # 绘制聚类中心
    centers_pca = pca.transform(kmeans.cluster_centers_)
    ax1.scatter(centers_pca[:, 0], centers_pca[:, 1], 
               c='yellow', marker='*', s=500, edgecolors='black', linewidth=2,
               label='Cluster Centers', zorder=5)
    
    # 添加聚类中心标签
    for i, (x, y) in enumerate(centers_pca):
        cluster_type = cluster_names.get(i, 'unknown')
        ax1.annotate(f'C{i}\n({cluster_type})', 
                    (x, y), xytext=(5, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', 
                    facecolor='yellow', alpha=0.7))
    
    ax1.set_xlabel(f'PC1 (解释方差: {pca.explained_variance_ratio_[0]:.2%})', fontsize=12)
    ax1.set_ylabel(f'PC2 (解释方差: {pca.explained_variance_ratio_[1]:.2%})', fontsize=12)
    ax1.set_title('K-means Clustering Results (PCA Visualization)', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 右图：特征分布对比（使用原始值，非标准化）
    ax2 = axes[1]
    
    # 准备特征数据用于可视化（使用原始值）
    feature_data_raw = {
        'Shaft Width': df_non_stubby['shaft_width'].values,
        'Head/Shaft Ratio': df_non_stubby['head_shaft_width_ratio'].values,
        'Shaft Length': df_non_stubby['shaft_length'].values,
        'Total Length': df_non_stubby['spine_length'].values
    }
    
    # 按类型分组数据
    plot_data = {spine_type: {} for spine_type in ['thin', 'mushroom', 'filopodia']}
    for spine_type in plot_data.keys():
        mask = spine_types == spine_type
        if mask.sum() > 0:
            for feat_name, feat_values in feature_data_raw.items():
                plot_data[spine_type][feat_name] = feat_values[mask]
    
    # 绘制分组箱线图
    feature_names = list(feature_data_raw.keys())
    n_features = len(feature_names)
    n_types = len([t for t in plot_data.keys() if len(plot_data[t]) > 0])
    
    positions = np.arange(n_features)
    width = 0.25
    
    box_plots = []
    for idx, spine_type in enumerate(['thin', 'mushroom', 'filopodia']):
        if len(plot_data[spine_type]) == 0:
            continue
        data_to_plot = [plot_data[spine_type][feat] for feat in feature_names]
        
        bp = ax2.boxplot(data_to_plot, positions=positions + idx*width, 
                       widths=width*0.8, patch_artist=True)
        
        for patch in bp['boxes']:
            patch.set_facecolor(colors_map[spine_type])
            patch.set_alpha(0.7)
        
        box_plots.append(bp)
    
    ax2.set_xticks(positions + width)
    ax2.set_xticklabels(feature_names, fontsize=10, rotation=15, ha='right')
    ax2.set_ylabel('Feature Value', fontsize=12)
    ax2.set_title('Feature Distribution by Spine Type', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors_map[t], alpha=0.7, label=t.capitalize()) 
                      for t in ['thin', 'mushroom', 'filopodia'] if len(plot_data.get(t, {})) > 0]
    ax2.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('kmeans_clustering_visualization.png', dpi=300, bbox_inches='tight')
    print("K-means聚类可视化已保存到: kmeans_clustering_visualization.png")
    
    return fig

def print_classification_logic():
    """打印详细的分类逻辑代码"""
    print("\n" + "="*80)
    print("详细的Spine分类逻辑代码")
    print("="*80)
    
    print("\n【1. Stubby类型识别逻辑】")
    print("""
def identify_stubby_spines(df):
    # 计算width相似度（差异百分比）
    width_diff_ratio = abs(df['shaft_width'] - df['head_width']) / (df['shaft_width'] + df['head_width'] + 1e-10)
    
    # 判断条件
    condition1 = df['head_length'] < df['dendrite_length']  # head length < dendrite length
    condition2 = width_diff_ratio < 0.2  # shaft width 和 head width 大致相似（差异<20%）
    condition3 = df['shaft_length'] < 0.1  # shaft length 几乎为0
    
    # stubby条件：满足(条件1 AND 条件2) OR 条件3
    is_stubby = (condition1 & condition2) | condition3
    
    # 对于stubby类型，length = head length, width = head width
    df.loc[is_stubby, 'spine_length_stubby'] = df.loc[is_stubby, 'head_length']
    df.loc[is_stubby, 'spine_width_stubby'] = df.loc[is_stubby, 'head_width']
    
    return df, is_stubby
    """)
    
    print("\n【2. K-means聚类逻辑】")
    print("""
def cluster_spines(df_non_stubby):
    # 准备特征
    features = df_non_stubby.copy()
    features['head_shaft_width_ratio'] = features['head_width'] / (features['shaft_width'] + 1e-10)
    
    # 选择用于聚类的特征
    X = features[['shaft_width', 'head_shaft_width_ratio', 'shaft_length', 'spine_length']].values
    
    # 标准化特征
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # K-means聚类，分成3类
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)
    
    return clusters, kmeans, X_scaled
    """)
    
    print("\n【3. 聚类命名逻辑】")
    print("""
def name_clusters(cluster_stats):
    # 计算每个聚类的特征值
    cluster_features = []
    for cluster_id in range(3):
        stats = cluster_stats.loc[cluster_id]
        cluster_features.append({
            'id': cluster_id,
            'shaft_length': stats['shaft_length'],
            'head_shaft_ratio': stats['head_shaft_width_ratio'],
            'head_width': stats['head_width'],
            'shaft_width': stats['shaft_width'],
            'spine_length': stats['spine_length']
        })
    
    # 根据特征值命名
    sorted_by_shaft_length = sorted(cluster_features, key=lambda x: x['shaft_length'], reverse=True)
    
    # 1. 如果shaft_length差异明显（>1.3倍），最长的为thin
    if sorted_by_shaft_length[0]['shaft_length'] > sorted_by_shaft_length[1]['shaft_length'] * 1.3:
        cluster_names[sorted_by_shaft_length[0]['id']] = 'thin'
        remaining = [c for c in cluster_features if c['id'] != sorted_by_shaft_length[0]['id']]
    else:
        remaining = cluster_features
        cluster_names[sorted_by_shaft_length[0]['id']] = 'filopodia'
    
    # 2. 在剩余的中，找出head_width/shaft_width比率最大的 -> mushroom
    sorted_remaining = sorted(remaining, key=lambda x: x['head_shaft_ratio'], reverse=True)
    if sorted_remaining[0]['head_shaft_ratio'] > 1.3:
        cluster_names[sorted_remaining[0]['id']] = 'mushroom'
        final_remaining = [c for c in remaining if c['id'] != sorted_remaining[0]['id']]
    else:
        final_remaining = remaining
        cluster_names[sorted_remaining[0]['id']] = 'filopodia'
    
    # 3. 最后一个默认为filopodia
    for c in final_remaining:
        if c['id'] not in cluster_names:
            cluster_names[c['id']] = 'filopodia'
    
    return cluster_names
    """)
    
    print("\n【4. 分类总结】")
    print("""
总共4种类型：
1. Stubby（短粗型）:
   - 条件: (head_length < dendrite_length AND shaft_width ≈ head_width) OR shaft_length < 0.1
   - 特征: length = head_length, width = head_width

2. Thin（细长型）:
   - 特征: shaft_length 最长（比其他类型长30%以上）
   - 特点: 细长的shaft，中等大小的head

3. Mushroom（蘑菇型）:
   - 特征: head_width/shaft_width 比率最大（>1.3）
   - 特点: 小的shaft，大的head

4. Filopodia（丝状）:
   - 特征: 其他情况
   - 特点: 中等大小的shaft和head
    """)
    
    print("="*80 + "\n")

def visualize_spine_types(df_all):
    """可视化所有spine类型及其平均数值"""
    # 按类型分组统计
    type_stats = df_all.groupby('spine_type').agg({
        'shaft_length': 'mean',
        'shaft_width': 'mean',
        'head_length': 'mean',
        'head_width': 'mean',
        'spine_length': 'mean'
    }).round(4)
    
    # 对于stubby类型，使用head_length和head_width作为length和width
    if 'stubby' in type_stats.index:
        type_stats.loc['stubby', 'spine_length'] = type_stats.loc['stubby', 'head_length']
    
    # 创建图形
    fig = plt.figure(figsize=(18, 12))
    
    # 第一部分：绘制每种类型的spine形状
    spine_types = sorted(df_all['spine_type'].unique())
    n_types = len(spine_types)
    
    # 创建子图布局：第一行显示spine形状，第二行显示统计表格
    gs = fig.add_gridspec(3, max(4, n_types), hspace=0.4, wspace=0.3)
    
    colors = {'stubby': 'skyblue', 'thin': 'lightgreen', 'mushroom': 'lightcoral', 'filopodia': 'lightyellow'}
    
    for idx, spine_type in enumerate(spine_types):
        ax = fig.add_subplot(gs[0, idx])
        stats = type_stats.loc[spine_type]
        draw_spine_shape(ax, spine_type, stats, x_offset=0, scale=0.6)
        count = len(df_all[df_all["spine_type"]==spine_type])
        ax.set_title(f'{spine_type.capitalize()} Type\n(n={count})', 
                    fontsize=12, fontweight='bold')
    
    # 第二部分：显示统计表格
    ax_table = fig.add_subplot(gs[1:, :])
    ax_table.axis('off')
    
    # 创建表格数据
    table_data = []
    for spine_type in spine_types:
        stats = type_stats.loc[spine_type]
        count = len(df_all[df_all['spine_type'] == spine_type])
        
        if spine_type == 'stubby':
            # stubby类型：length = head_length, width = head_width
            table_data.append([
                spine_type.capitalize(),
                count,
                'N/A (stubby)',
                'N/A (stubby)',
                f"{stats['head_length']:.4f}",
                f"{stats['head_width']:.4f}",
                f"{stats['head_length']:.4f}",
                f"{stats['head_width']:.4f}"
            ])
        else:
            table_data.append([
                spine_type.capitalize(),
                count,
                f"{stats['shaft_length']:.4f}",
                f"{stats['shaft_width']:.4f}",
                f"{stats['head_length']:.4f}",
                f"{stats['head_width']:.4f}",
                f"{stats['spine_length']:.4f}",
                'N/A'
            ])
    
    table = ax_table.table(cellText=table_data,
                          colLabels=['Spine Type', 'Count', 'Shaft Length', 'Shaft Width', 
                                    'Head Length', 'Head Width', 'Total Length', 'Width (stubby)'],
                          cellLoc='center',
                          loc='center',
                          bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.5)
    
    # 设置表头样式
    for i in range(len(table_data[0])):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # 设置不同行的颜色
    for i in range(1, len(table_data) + 1):
        for j in range(len(table_data[0])):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')
            else:
                table[(i, j)].set_facecolor('#ffffff')
    
    plt.suptitle('Spine Classification and Statistics', fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig('spine_classification_results.png', dpi=300, bbox_inches='tight')
    print("可视化结果已保存到: spine_classification_results.png")
    
    return type_stats

def main():
    # 读取详细spine数据
    print("正在读取spine详细数据...")
    detail_files = list(Path('.').glob('dendrite_spine_details*.csv'))
    if not detail_files:
        print("错误：未找到spine详细数据文件！")
        return
    
    # 使用最新的文件
    latest_file = max(detail_files, key=lambda p: p.stat().st_mtime)
    print(f"读取文件: {latest_file}")
    
    df = pd.read_csv(latest_file)
    print(f"总共 {len(df)} 个spine")
    
    # 识别stubby类型
    print("\n正在识别stubby类型spine...")
    df, is_stubby = identify_stubby_spines(df)
    stubby_count = is_stubby.sum()
    print(f"识别出 {stubby_count} 个stubby类型spine")
    
    # 对非stubby进行聚类
    print("\n正在对非stubby spine进行聚类分析...")
    df_non_stubby = df[~is_stubby].copy()
    print(f"剩余 {len(df_non_stubby)} 个spine用于聚类")
    
    if len(df_non_stubby) > 3:
        df_non_stubby, cluster_stats, cluster_names, X_scaled, scaler, kmeans = cluster_spines(df_non_stubby)
        
        # 更新原始数据框
        df.loc[~is_stubby, 'spine_type'] = df_non_stubby['spine_type'].values
        
        print("\n聚类结果:")
        for cluster_id, cluster_name in cluster_names.items():
            count = len(df_non_stubby[df_non_stubby['cluster'] == cluster_id])
            print(f"  {cluster_name.capitalize()}: {count} 个spine")
        
        # 生成K-means聚类可视化
        print("\n正在生成K-means聚类可视化...")
        visualize_kmeans_clustering(df_non_stubby, X_scaled, scaler, kmeans, cluster_names)
    else:
        print("警告：非stubby spine数量太少，无法进行聚类分析")
        df.loc[~is_stubby, 'spine_type'] = 'other'
    
    # 保存分类结果
    output_file = 'spine_classification_results.csv'
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n分类结果已保存到: {output_file}")
    
    # 可视化
    print("\n正在生成可视化图表...")
    type_stats = visualize_spine_types(df)
    
    # 打印详细的分类逻辑
    print_classification_logic()
    
    # 打印统计摘要
    print("\n" + "="*60)
    print("分类统计摘要:")
    print("="*60)
    print(type_stats)
    print("\n各类型spine数量:")
    print(df['spine_type'].value_counts())
    
    return df, type_stats

if __name__ == '__main__':
    main()

