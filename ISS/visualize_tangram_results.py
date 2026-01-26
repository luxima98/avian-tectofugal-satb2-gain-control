#!/usr/bin/env python3
"""
Tangram结果可视化脚本
展示anno_level_3_merged分类的细胞空间分布，使用真实坐标
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import tangram as tg
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_tangram_results():
    """加载Tangram映射结果"""
    base_dir = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data"
    spatial_dir = os.path.join(base_dir, "spatial", "spatial_genesymbol_converted")
    results_dir = os.path.join(spatial_dir, "V", "official_tangram_results")
    tangram_matrices_dir = os.path.join(spatial_dir, "V", "baysor_tangram_output", "tangram_matrices")
    
    # 加载数据
    sc_file = os.path.join(base_dir, "snRNAseq", "seurat_obj_merged.h5ad")
    spatial_file = os.path.join(tangram_matrices_dir, "spatial_cells.h5ad")
    cells_result_file = os.path.join(results_dir, "tangram_cells_result.h5ad")
    clusters_result_file = os.path.join(results_dir, "tangram_clusters_result.h5ad")
    
    print("加载数据...")
    ad_sc = sc.read_h5ad(sc_file)
    
    # 加载空间数据（用于获取坐标信息）
    if os.path.exists(spatial_file):
        ad_spatial = sc.read_h5ad(spatial_file)
        print("✓ 加载空间数据用于坐标信息")
    else:
        print("⚠️ 未找到空间数据文件，将无法显示空间坐标")
        ad_spatial = None
    
    # 检查cells结果文件是否存在
    if os.path.exists(cells_result_file):
        ad_map_cells = sc.read_h5ad(cells_result_file)
        print("✓ 找到cells级别映射结果")
    else:
        print("⚠️ 未找到tangram_cells_result.h5ad文件，将使用clusters结果")
        ad_map_cells = None
    
    ad_map_clusters = sc.read_h5ad(clusters_result_file)
    
    print(f"✓ 单细胞数据: {ad_sc.shape}")
    if ad_spatial is not None:
        print(f"✓ 空间数据: {ad_spatial.shape}")
    else:
        print("✓ 空间数据: 未找到")
    if ad_map_cells is not None:
        print(f"✓ Cells映射结果: {ad_map_cells.shape}")
    else:
        print("✓ Cells映射结果: 未找到")
    print(f"✓ Clusters映射结果: {ad_map_clusters.shape}")
    
    return ad_sc, ad_spatial, ad_map_cells, ad_map_clusters

def load_filtered_spatial_data():
    """加载筛选后的空间数据"""
    base_dir = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data"
    spatial_dir = os.path.join(base_dir, "spatial", "spatial_genesymbol_converted")
    assignment_dir = os.path.join(spatial_dir, "V", "cell_type_assignment")
    
    # 加载筛选后的空间数据
    filtered_spatial_file = os.path.join(assignment_dir, "spatial_cells_with_assignments.h5ad")
    
    if os.path.exists(filtered_spatial_file):
        ad_spatial_filtered = sc.read_h5ad(filtered_spatial_file)
        print("✓ 加载筛选后的空间数据")
        print(f"✓ 筛选后空间数据: {ad_spatial_filtered.shape}")
        print(f"✓ 细胞类型分配列: {list(ad_spatial_filtered.obs.columns)}")
        return ad_spatial_filtered
    else:
        print("⚠️ 未找到筛选后的空间数据文件")
        return None

def plot_cell_type_distribution(ad_map_clusters, ad_spatial, ad_sc, output_dir):
    """绘制细胞类型在空间中的分布"""
    print("\n1. 绘制细胞类型空间分布...")
    
    # 获取细胞类型索引和真实名称的映射
    cell_type_indices = ad_map_clusters.obs_names  # 这些是数字索引
    print(f"   细胞类型索引: {list(cell_type_indices)}")
    
    # 从Tangram映射结果中获取真实的细胞类型名称（这是正确的映射）
    if 'anno_level_3_merged' in ad_map_clusters.obs.columns:
        tangram_cell_types = ad_map_clusters.obs['anno_level_3_merged'].tolist()
        print(f"   Tangram中的细胞类型数: {len(tangram_cell_types)}")
        print(f"   Tangram中的细胞类型: {tangram_cell_types}")
        
        # 创建索引到真实名称的映射（使用Tangram保存的正确映射）
        index_to_name = {}
        for idx, cell_type in zip(cell_type_indices, tangram_cell_types):
            index_to_name[idx] = cell_type
    else:
        # 如果没有保存的细胞类型信息，使用单细胞数据的顺序（不推荐）
        unique_cell_types = ad_sc.obs['anno_level_3_merged'].unique()
        print(f"   ⚠️ 使用单细胞数据顺序（可能不正确）")
        print(f"   单细胞细胞类型数: {len(unique_cell_types)}")
        print(f"   单细胞细胞类型: {list(unique_cell_types)}")
        
        index_to_name = {}
        for i, idx in enumerate(cell_type_indices):
            if i < len(unique_cell_types):
                index_to_name[idx] = unique_cell_types[i]
            else:
                index_to_name[idx] = f"Type_{idx}"
    
    # 获取空间坐标
    if 'spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['spatial']
    elif 'X_spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['X_spatial']
    else:
        # 尝试从obs中获取坐标
        if 'x' in ad_spatial.obs and 'y' in ad_spatial.obs:
            coords = np.column_stack([ad_spatial.obs['x'], ad_spatial.obs['y']])
        else:
            print("   ❌ 未找到空间坐标信息")
            return None
    
    print(f"   空间坐标形状: {coords.shape}")
    
    # 创建子图
    n_types = len(cell_type_indices)
    n_cols = 4
    n_rows = (n_types + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 为每个细胞类型绘制分布
    for i, cell_type_idx in enumerate(cell_type_indices):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # 获取该细胞类型的映射概率
        cell_type_probs = ad_map_clusters.X[i, :]
        
        # 创建散点图，使用反转的颜色映射（概率越高颜色越深）
        # 使用百分位数来改善颜色对比度
        vmin = np.percentile(cell_type_probs, 1)  # 1%分位数
        vmax = np.percentile(cell_type_probs, 99)  # 99%分位数
        
        scatter = ax.scatter(coords[:, 0], 
                           coords[:, 1],
                           c=cell_type_probs, 
                           cmap='viridis_r',  # 反转颜色映射
                           s=1, 
                           alpha=0.8,
                           vmin=vmin,
                           vmax=vmax)
        
        # 使用真实的细胞类型名称
        real_name = index_to_name[cell_type_idx]
        ax.set_title(f'{real_name}', fontsize=10)
        ax.set_xlabel('X Coordinate (pixel)')
        ax.set_ylabel('Y Coordinate (pixel)')
        ax.set_aspect('equal')
        
        # 添加颜色条
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label('Mapping Probability', fontsize=8)
    
    # 隐藏多余的子图
    for i in range(n_types, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, "cell_type_spatial_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ 细胞类型空间分布图已保存: {output_file}")
    
    # 保存PDF版本
    output_file_pdf = os.path.join(output_dir, "cell_type_spatial_distribution.pdf")
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ 细胞类型空间分布图PDF已保存: {output_file_pdf}")
    
    plt.close()
    
    return output_file

def plot_dominant_cell_types(ad_map_clusters, ad_spatial, ad_sc, output_dir, top_n=10):
    """绘制主要细胞类型的分布"""
    print(f"\n2. 绘制前{top_n}个主要细胞类型分布...")
    
    # 获取空间坐标
    if 'spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['spatial']
    elif 'X_spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['X_spatial']
    else:
        if 'x' in ad_spatial.obs and 'y' in ad_spatial.obs:
            coords = np.column_stack([ad_spatial.obs['x'], ad_spatial.obs['y']])
        else:
            print("   ❌ 未找到空间坐标信息")
            return None
    
    # 使用Tangram保存的正确细胞类型信息
    if 'anno_level_3_merged' in ad_map_clusters.obs.columns:
        tangram_cell_types = ad_map_clusters.obs['anno_level_3_merged'].tolist()
        cell_type_indices = ad_map_clusters.obs_names
        
        # 创建索引到真实名称的映射（使用Tangram保存的正确映射）
        index_to_name = {}
        for idx, cell_type in zip(cell_type_indices, tangram_cell_types):
            index_to_name[idx] = cell_type
    else:
        # 如果没有保存的细胞类型信息，使用单细胞数据的顺序（不推荐）
        unique_cell_types = ad_sc.obs['anno_level_3_merged'].unique()
        cell_type_indices = ad_map_clusters.obs_names
        
        index_to_name = {}
        for i, idx in enumerate(cell_type_indices):
            if i < len(unique_cell_types):
                index_to_name[idx] = unique_cell_types[i]
            else:
                index_to_name[idx] = f"Type_{idx}"
    
    # 计算每个细胞类型的平均映射概率（而不是总和，因为每行总和都是1）
    cell_type_means = np.mean(ad_map_clusters.X, axis=1)
    
    # 获取前N个细胞类型
    top_indices = np.argsort(cell_type_means)[-top_n:][::-1]
    top_cell_type_indices = [cell_type_indices[i] for i in top_indices]
    top_means = [cell_type_means[i] for i in top_indices]
    
    print(f"   前{top_n}个细胞类型:")
    for i, (cell_type_idx, mean_prob) in enumerate(zip(top_cell_type_indices, top_means)):
        real_name = index_to_name[cell_type_idx]
        print(f"   {i+1}. {real_name} (索引{cell_type_idx}): {mean_prob:.4f}")
    
    # 创建子图
    n_cols = 3
    n_rows = (top_n + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 为每个主要细胞类型绘制分布
    for i, cell_type_idx in enumerate(top_cell_type_indices):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # 获取该细胞类型的映射概率
        cell_type_probs = ad_map_clusters.X[top_indices[i], :]
        
        # 创建散点图，使用反转的颜色映射
        # 使用百分位数改善颜色对比度
        vmin = np.percentile(cell_type_probs, 1)  # 1%分位数
        vmax = np.percentile(cell_type_probs, 99)  # 99%分位数
        
        scatter = ax.scatter(coords[:, 0], 
                           coords[:, 1],
                           c=cell_type_probs, 
                           cmap='plasma_r',  # 反转颜色映射
                           s=2, 
                           alpha=0.8,
                           vmin=vmin,
                           vmax=vmax)
        
        real_name = index_to_name[cell_type_idx]
        ax.set_title(f'{real_name}\n(Mean Probability: {top_means[i]:.4f})', fontsize=10)
        ax.set_xlabel('X Coordinate (pixel)')
        ax.set_ylabel('Y Coordinate (pixel)')
        ax.set_aspect('equal')
        
        # 添加颜色条
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label('Mapping Probability', fontsize=8)
    
    # 隐藏多余的子图
    for i in range(top_n, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, f"top_{top_n}_cell_types_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ 前{top_n}个细胞类型分布图已保存: {output_file}")
    
    # 保存PDF版本
    output_file_pdf = os.path.join(output_dir, f"top_{top_n}_cell_types_distribution.pdf")
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ 前{top_n}个细胞类型分布图PDF已保存: {output_file_pdf}")
    
    plt.close()
    
    return output_file

def plot_cell_type_abundance(ad_map_clusters, ad_sc, output_dir):
    """绘制细胞类型丰度图"""
    print("\n3. 绘制细胞类型丰度图...")
    
    # 使用Tangram保存的正确细胞类型信息
    if 'anno_level_3_merged' in ad_map_clusters.obs.columns:
        tangram_cell_types = ad_map_clusters.obs['anno_level_3_merged'].tolist()
        cell_type_indices = ad_map_clusters.obs_names
        
        # 创建索引到真实名称的映射（使用Tangram保存的正确映射）
        index_to_name = {}
        for idx, cell_type in zip(cell_type_indices, tangram_cell_types):
            index_to_name[idx] = cell_type
    else:
        # 如果没有保存的细胞类型信息，使用单细胞数据的顺序（不推荐）
        unique_cell_types = ad_sc.obs['anno_level_3_merged'].unique()
        cell_type_indices = ad_map_clusters.obs_names
        
        index_to_name = {}
        for i, idx in enumerate(cell_type_indices):
            if i < len(unique_cell_types):
                index_to_name[idx] = unique_cell_types[i]
            else:
                index_to_name[idx] = f"Type_{idx}"
    
    # 计算每个细胞类型的平均映射概率
    cell_type_means = np.mean(ad_map_clusters.X, axis=1)
    
    # 创建DataFrame，使用真实细胞类型名称
    real_names = [index_to_name[idx] for idx in cell_type_indices]
    df = pd.DataFrame({
        'Cell_Type': real_names,
        'Cell_Type_Index': cell_type_indices,
        'Mean_Probability': cell_type_means
    }).sort_values('Mean_Probability', ascending=True)
    
    # 创建水平条形图
    plt.figure(figsize=(14, 8))
    bars = plt.barh(range(len(df)), df['Mean_Probability'], color='skyblue', alpha=0.7)
    
    plt.yticks(range(len(df)), df['Cell_Type'])
    plt.xlabel('Mean Mapping Probability')
    plt.title('Cell Type Abundance Distribution (anno_level_3_merged)')
    plt.grid(axis='x', alpha=0.3)
    
    # 添加数值标签
    for i, (bar, value) in enumerate(zip(bars, df['Mean_Probability'])):
        plt.text(value + 0.0001, bar.get_y() + bar.get_height()/2, 
                f'{value:.4f}', va='center', fontsize=8)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, "cell_type_abundance.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ 细胞类型丰度图已保存: {output_file}")
    
    # 保存PDF版本
    output_file_pdf = os.path.join(output_dir, "cell_type_abundance.pdf")
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ 细胞类型丰度图PDF已保存: {output_file_pdf}")
    
    plt.close()
    
    return output_file

def plot_spatial_heatmap(ad_map_clusters, ad_spatial, output_dir):
    """绘制空间热图"""
    print("\n4. 绘制空间热图...")
    
    # 获取空间坐标
    if 'spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['spatial']
    elif 'X_spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['X_spatial']
    else:
        if 'x' in ad_spatial.obs and 'y' in ad_spatial.obs:
            coords = np.column_stack([ad_spatial.obs['x'], ad_spatial.obs['y']])
        else:
            print("   ❌ 未找到空间坐标信息")
            return None
    
    # 创建网格
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    
    # 创建网格点
    grid_size = 50
    x_grid = np.linspace(x_min, x_max, grid_size)
    y_grid = np.linspace(y_min, y_max, grid_size)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
    
    # 为每个细胞类型创建热图
    cell_types = ad_map_clusters.obs_names
    n_types = len(cell_types)
    n_cols = 4
    n_rows = (n_types + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, cell_type in enumerate(cell_types):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # 获取该细胞类型的映射概率
        cell_type_probs = ad_map_clusters.X[i, :]
        
        # 插值到网格
        from scipy.interpolate import griddata
        Z = griddata(coords, cell_type_probs, (X_grid, Y_grid), method='cubic', fill_value=0)
        
        # 绘制热图
        im = ax.contourf(X_grid, Y_grid, Z, levels=20, cmap='viridis', alpha=0.8)
        ax.contour(X_grid, Y_grid, Z, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        
        ax.set_title(f'{cell_type}', fontsize=10)
        ax.set_xlabel('X Coordinate (pixel)')
        ax.set_ylabel('Y Coordinate (pixel)')
        ax.set_aspect('equal')
        
        # 添加颜色条
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    # 隐藏多余的子图
    for i in range(n_types, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, "spatial_heatmap.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ 空间热图已保存: {output_file}")
    
    # 保存PDF版本
    output_file_pdf = os.path.join(output_dir, "spatial_heatmap.pdf")
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ 空间热图PDF已保存: {output_file_pdf}")
    
    plt.close()
    
    return output_file

def create_summary_report(ad_map_clusters, ad_spatial, output_dir):
    """创建总结报告"""
    print("\n5. 创建总结报告...")
    
    # 获取空间坐标
    if 'spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['spatial']
    elif 'X_spatial' in ad_spatial.obsm:
        coords = ad_spatial.obsm['X_spatial']
    else:
        if 'x' in ad_spatial.obs and 'y' in ad_spatial.obs:
            coords = np.column_stack([ad_spatial.obs['x'], ad_spatial.obs['y']])
        else:
            coords = None
    
    # 计算统计信息
    cell_type_means = np.mean(ad_map_clusters.X, axis=1)
    cell_type_names = ad_map_clusters.obs_names
    
    # 创建DataFrame
    df = pd.DataFrame({
        'Cell_Type_Index': cell_type_names,
        'Mean_Probability': cell_type_means,
        'Max_Probability': np.max(ad_map_clusters.X, axis=1),
        'Std_Probability': np.std(ad_map_clusters.X, axis=1)
    }).sort_values('Mean_Probability', ascending=False)
    
    # 保存统计表
    stats_file = os.path.join(output_dir, "cell_type_statistics.csv")
    df.to_csv(stats_file, index=False, encoding='utf-8-sig')
    print(f"   ✓ 细胞类型统计表已保存: {stats_file}")
    
    # 创建文本报告
    report_file = os.path.join(output_dir, "visualization_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Tangram可视化结果报告\n")
        f.write("=" * 40 + "\n\n")
        
        f.write("数据概览:\n")
        f.write(f"  细胞类型数: {len(cell_type_names)}\n")
        f.write(f"  空间位置数: {ad_map_clusters.shape[1]}\n")
        if coords is not None:
            f.write(f"  坐标范围: X({coords[:, 0].min():.2f}, {coords[:, 0].max():.2f})\n")
            f.write(f"            Y({coords[:, 1].min():.2f}, {coords[:, 1].max():.2f})\n\n")
        else:
            f.write("  坐标信息: 未找到\n\n")
        
        f.write("前10个主要细胞类型:\n")
        for i, (_, row) in enumerate(df.head(10).iterrows()):
            f.write(f"  {i+1}. 索引{row['Cell_Type_Index']}: 平均概率={row['Mean_Probability']:.4f}, 最大概率={row['Max_Probability']:.4f}\n")
        
        f.write("\n生成的可视化文件:\n")
        f.write("  1. cell_type_spatial_distribution.png - 所有细胞类型的空间分布\n")
        f.write("  2. top_10_cell_types_distribution.png - 前10个主要细胞类型分布\n")
        f.write("  3. cell_type_abundance.png - 细胞类型丰度条形图\n")
        f.write("  4. spatial_heatmap.png - 空间热图\n")
        f.write("  5. cell_type_statistics.csv - 细胞类型统计表\n")
        
        f.write("\n可视化特点:\n")
        f.write("  - 使用anno_level_3_merged的真实细胞类型名称\n")
        f.write("  - 反转颜色映射：概率越高颜色越深\n")
        f.write("  - 基于真实空间坐标绘制\n")
    
    print(f"   ✓ 可视化报告已保存: {report_file}")
    
    return report_file

def plot_filtered_cell_type_assignment(ad_spatial_filtered, output_dir):
    """绘制筛选后的细胞类型分配结果"""
    print("\n6. 绘制筛选后的细胞类型分配...")
    
    if ad_spatial_filtered is None:
        print("   ⚠️ 没有筛选后的数据，跳过此步骤")
        return None
    
    # 获取空间坐标
    if 'x' in ad_spatial_filtered.obs and 'y' in ad_spatial_filtered.obs:
        coords = np.column_stack([ad_spatial_filtered.obs['x'], ad_spatial_filtered.obs['y']])
    else:
        print("   ❌ 未找到空间坐标信息")
        return None
    
    # 获取细胞类型分配信息
    assigned_cell_types = ad_spatial_filtered.obs['assigned_cell_type']
    assignment_probabilities = ad_spatial_filtered.obs['assignment_probability']
    is_top_75 = ad_spatial_filtered.obs['is_top_75_percent']
    
    print(f"   空间位置数: {len(coords)}")
    print(f"   细胞类型数: {assigned_cell_types.nunique()}")
    print(f"   前75%位置数: {is_top_75.sum()}")
    
    # 创建细胞类型颜色映射
    unique_cell_types = assigned_cell_types.unique()
    n_types = len(unique_cell_types)
    
    # 使用tab20颜色映射
    colors = plt.cm.tab20(np.linspace(0, 1, n_types))
    cell_type_colors = dict(zip(unique_cell_types, colors))
    
    # 创建子图
    n_cols = 3
    n_rows = (n_types + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 6*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 为每个细胞类型绘制分布
    for i, cell_type in enumerate(unique_cell_types):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # 获取该细胞类型的位置
        cell_type_mask = assigned_cell_types == cell_type
        cell_type_coords = coords[cell_type_mask]
        cell_type_probs = assignment_probabilities[cell_type_mask]
        cell_type_top_75 = is_top_75[cell_type_mask]
        
        # 绘制所有该细胞类型的位置（浅色）
        ax.scatter(cell_type_coords[:, 0], 
                  cell_type_coords[:, 1],
                  c=cell_type_colors[cell_type], 
                  s=1, 
                  alpha=0.3,
                  label=f'All {cell_type}')
        
        # 绘制前75%的位置（深色）
        if cell_type_top_75.sum() > 0:
            top_75_coords = cell_type_coords[cell_type_top_75]
            top_75_probs = cell_type_probs[cell_type_top_75]
            
            scatter = ax.scatter(top_75_coords[:, 0], 
                               top_75_coords[:, 1],
                               c=top_75_probs, 
                               cmap='viridis', 
                               s=2, 
                               alpha=0.8,
                               vmin=top_75_probs.min(),
                               vmax=top_75_probs.max())
            
            # 添加颜色条
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label('Assignment Probability', fontsize=8)
        
        ax.set_title(f'{cell_type}\n(Total: {cell_type_mask.sum()}, Top75%: {cell_type_top_75.sum()})', fontsize=10)
        ax.set_xlabel('X Coordinate (pixel)')
        ax.set_ylabel('Y Coordinate (pixel)')
        ax.set_aspect('equal')
    
    # 隐藏多余的子图
    for i in range(n_types, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, "filtered_cell_type_assignment.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✓ 筛选后细胞类型分配图已保存: {output_file}")
    
    return output_file

def plot_cell_type_assignment_overview(ad_spatial_filtered, output_dir):
    """绘制细胞类型分配概览"""
    print("\n7. 绘制细胞类型分配概览...")
    
    if ad_spatial_filtered is None:
        print("   ⚠️ 没有筛选后的数据，跳过此步骤")
        return None
    
    # 获取空间坐标
    if 'x' in ad_spatial_filtered.obs and 'y' in ad_spatial_filtered.obs:
        coords = np.column_stack([ad_spatial_filtered.obs['x'], ad_spatial_filtered.obs['y']])
    else:
        print("   ❌ 未找到空间坐标信息")
        return None
    
    # 获取细胞类型分配信息
    assigned_cell_types = ad_spatial_filtered.obs['assigned_cell_type']
    is_top_75 = ad_spatial_filtered.obs['is_top_75_percent']
    
    # 创建细胞类型颜色映射
    unique_cell_types = assigned_cell_types.unique()
    n_types = len(unique_cell_types)
    colors = plt.cm.tab20(np.linspace(0, 1, n_types))
    cell_type_colors = dict(zip(unique_cell_types, colors))
    
    # 创建两个子图：所有分配 vs 前75%分配
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 左图：所有细胞类型分配
    for cell_type in unique_cell_types:
        cell_type_mask = assigned_cell_types == cell_type
        cell_type_coords = coords[cell_type_mask]
        
        ax1.scatter(cell_type_coords[:, 0], 
                   cell_type_coords[:, 1],
                   c=[cell_type_colors[cell_type]], 
                   s=1, 
                   alpha=0.7,
                   label=cell_type)
    
    ax1.set_title('All Cell Type Assignments', fontsize=14)
    ax1.set_xlabel('X Coordinate (pixel)')
    ax1.set_ylabel('Y Coordinate (pixel)')
    ax1.set_aspect('equal')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # 右图：前75%细胞类型分配
    for cell_type in unique_cell_types:
        cell_type_mask = (assigned_cell_types == cell_type) & is_top_75
        cell_type_coords = coords[cell_type_mask]
        
        if len(cell_type_coords) > 0:
            ax2.scatter(cell_type_coords[:, 0], 
                       cell_type_coords[:, 1],
                       c=[cell_type_colors[cell_type]], 
                       s=2, 
                       alpha=0.8,
                       label=cell_type)
    
    ax2.set_title('Top 75% Cell Type Assignments', fontsize=14)
    ax2.set_xlabel('X Coordinate (pixel)')
    ax2.set_ylabel('Y Coordinate (pixel)')
    ax2.set_aspect('equal')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = os.path.join(output_dir, "cell_type_assignment_overview.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✓ 细胞类型分配概览图已保存: {output_file}")
    
    return output_file

def main():
    """主函数"""
    print("=" * 80)
    print("🎨 Tangram结果可视化")
    print("=" * 80)
    
    # 设置输出目录
    base_dir = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data"
    spatial_dir = os.path.join(base_dir, "spatial", "spatial_genesymbol_converted")
    output_dir = os.path.join(spatial_dir, "V", "tangram_visualization_filtered")
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 加载数据
        ad_sc, ad_spatial, ad_map_cells, ad_map_clusters = load_tangram_results()
        
        # 加载筛选后的空间数据
        ad_spatial_filtered = load_filtered_spatial_data()
        
        # 生成可视化
        plot_cell_type_distribution(ad_map_clusters, ad_spatial, ad_sc, output_dir)
        plot_dominant_cell_types(ad_map_clusters, ad_spatial, ad_sc, output_dir, top_n=10)
        plot_cell_type_abundance(ad_map_clusters, ad_sc, output_dir)
        plot_spatial_heatmap(ad_map_clusters, ad_spatial, output_dir)
        create_summary_report(ad_map_clusters, ad_spatial, output_dir)
        
        # 生成筛选后的可视化
        if ad_spatial_filtered is not None:
            plot_filtered_cell_type_assignment(ad_spatial_filtered, output_dir)
            plot_cell_type_assignment_overview(ad_spatial_filtered, output_dir)
        
        print("\n" + "=" * 80)
        print("🎯 Tangram可视化完成!")
        print("=" * 80)
        print(f"📁 输出目录: {output_dir}")
        print("📊 生成的文件:")
        print("  - cell_type_spatial_distribution.png")
        print("  - top_10_cell_types_distribution.png") 
        print("  - cell_type_abundance.png")
        print("  - spatial_heatmap.png")
        print("  - cell_type_statistics.csv")
        print("  - visualization_report.txt")
        
    except Exception as e:
        print(f"❌ 可视化失败: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
