#!/usr/bin/env python3
"""
按照官方Tangram notebook流程运行映射
与tutorial_tangram_with_squidpy.ipynb和tutorial_tangram_without_squidpy.ipynb保持一致
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

def convert_ensembl_to_symbol(spatial_data, mapping_file):
    """将ENSEMBL基因ID转换为Gene Symbol"""
    print("\n0. 转换ENSEMBL基因ID为Gene Symbol...")
    
    # 读取基因映射文件
    mapping_df = pd.read_excel(mapping_file)
    print(f"   ✓ 读取基因映射文件: {mapping_df.shape}")
    print(f"   ✓ 映射文件列名: {list(mapping_df.columns)}")
    
    # 创建ENSEMBL到Symbol的映射字典
    ensembl_to_symbol = {}
    for _, row in mapping_df.iterrows():
        ensembl_id = row['ENSEMBL GeneID ']  # 注意列名末尾有空格
        gene_symbol = row['Gene Symbol']
        if pd.notna(ensembl_id) and pd.notna(gene_symbol):
            ensembl_to_symbol[ensembl_id] = gene_symbol
    
    print(f"   ✓ 创建了 {len(ensembl_to_symbol)} 个基因映射")
    
    # 转换空间数据的基因名
    original_genes = spatial_data.var_names.tolist()
    converted_genes = []
    
    for gene in original_genes:
        if gene in ensembl_to_symbol:
            converted_genes.append(ensembl_to_symbol[gene])
        else:
            converted_genes.append(gene)  # 如果找不到映射，保持原名
    
    # 更新基因名
    spatial_data.var_names = converted_genes
    spatial_data.var.index = converted_genes
    
    print(f"   ✓ 转换前基因数: {len(original_genes)}")
    print(f"   ✓ 转换后基因数: {len(converted_genes)}")
    print(f"   ✓ 成功转换的基因数: {sum(1 for orig, conv in zip(original_genes, converted_genes) if orig != conv)}")
    
    return spatial_data

def run_tangram_official_flow():
    """按照官方notebook流程运行Tangram映射"""
    
    print("=" * 80)
    print("🚀 按照官方Tangram notebook流程运行映射")
    print("=" * 80)
    
    # 设置路径
    base_dir = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data"
    spatial_dir = os.path.join(base_dir, "spatial", "spatial_genesymbol_converted")
    baysor_output_dir = os.path.join(spatial_dir, "IV_BM_CR_25_baysor_tangram_output")
    tangram_matrices_dir = os.path.join(baysor_output_dir, "tangram_matrices")
    output_dir = os.path.join(spatial_dir, "IV_BM_CR_25_official_tangram_results")
    os.makedirs(output_dir, exist_ok=True)
    
    # 基因映射文件路径
    mapping_file = r"D:\r\data\chicken\science.adp5182_data_s2.xlsx"
    
    # 1. 加载数据 (对应官方notebook的Cell 5-12)
    print("\n1. 加载数据...")
    
    # 使用您提供的单细胞数据路径
    sc_file = os.path.join(base_dir, "snRNAseq", "seurat_obj_merged.h5ad")
    # 使用Baysor生成的空间细胞数据
    spatial_file = os.path.join(tangram_matrices_dir, "spatial_cells.h5ad")
    
    # 检查文件是否存在
    if not os.path.exists(sc_file):
        print(f"   ❌ 单细胞数据文件不存在: {sc_file}")
        return None
    
    if not os.path.exists(spatial_file):
        print(f"   ❌ 空间数据文件不存在: {spatial_file}")
        return None
    
    # 加载数据
    ad_sc = sc.read_h5ad(sc_file)
    ad_spatial = sc.read_h5ad(spatial_file)
    
    print(f"   ✓ 单细胞数据: {ad_sc.shape}")
    print(f"   ✓ 空间数据: {ad_spatial.shape}")
    print(f"   ✓ 细胞类型数: {ad_sc.obs['anno_level_3_merged'].nunique()}")
    
    # 转换ENSEMBL基因ID为Gene Symbol
    ad_spatial = convert_ensembl_to_symbol(ad_spatial, mapping_file)
    
    # 2. 数据预处理 (对应官方notebook的Cell 16)
    print("\n2. 数据预处理...")
    
    # 标准化 - 与官方notebook一致
    sc.pp.normalize_total(ad_sc, target_sum=1e4)
    sc.pp.log1p(ad_sc)
    
    print(f"   ✓ 单细胞数据预处理完成")
    
    # 3. 选择训练基因 (对应官方notebook的Cell 12-13)
    print("\n3. 选择训练基因...")
    
    # 大小写不敏感的基因匹配
    sc_genes_lower = [gene.lower() for gene in ad_sc.var_names]
    spatial_genes_lower = [gene.lower() for gene in ad_spatial.var_names]
    
    # 创建大小写不敏感的基因映射
    sc_gene_mapping = {gene.lower(): gene for gene in ad_sc.var_names}
    spatial_gene_mapping = {gene.lower(): gene for gene in ad_spatial.var_names}
    
    # 找到重叠基因（大小写不敏感）
    overlap_genes_lower = list(set(sc_genes_lower) & set(spatial_genes_lower))
    print(f"   ✓ 大小写不敏感重叠基因数: {len(overlap_genes_lower)}")
    
    # 转换为原始基因名
    overlap_genes = [sc_gene_mapping[gene_lower] for gene_lower in overlap_genes_lower]
    print(f"   ✓ 原始基因名重叠基因数: {len(overlap_genes)}")
    
    # 使用所有重叠基因作为训练基因
    training_genes = overlap_genes
    print(f"   ✓ 训练基因数: {len(training_genes)}")
    
    # 显示一些重叠基因示例
    if len(overlap_genes) > 0:
        print(f"   ✓ 重叠基因示例: {overlap_genes[:10]}")
    else:
        print("   ⚠ 警告: 没有找到重叠基因!")
        return None
    
    # 4. 准备数据 (对应官方notebook的Cell 14)
    print("\n4. 准备数据...")
    
    # 使用官方推荐的pp_adatas函数
    tg.pp_adatas(ad_sc, ad_spatial, genes=training_genes)
    
    print(f"   ✓ 数据准备完成")
    print(f"   ✓ 训练基因数: {len(ad_sc.uns['training_genes'])}")
    print(f"   ✓ 重叠基因数: {len(ad_sc.uns['overlap_genes'])}")
    
    # 5. 运行映射 - Cells模式 (对应官方notebook的Cell 18)
    print("\n5. 运行映射 - Cells模式...")
    
    try:
        ad_map_cells = tg.map_cells_to_space(
            ad_sc, 
            ad_spatial,
            mode="cells",
            density_prior='rna_count_based',
            num_epochs=500,
            device='cpu',  # 使用CPU，与官方notebook一致
            verbose=True
        )
        
        print(f"   ✅ Cells模式映射成功: {ad_map_cells.shape}")
        
        # 保存结果
        cells_file = os.path.join(output_dir, "tangram_cells_result.h5ad")
        ad_map_cells.write_h5ad(cells_file)
        print(f"   ✓ Cells模式结果已保存: {cells_file}")
        
    except Exception as e:
        print(f"   ❌ Cells模式映射失败: {e}")
        ad_map_cells = None
    
    # 6. 运行映射 - Clusters模式 (对应官方notebook的Cell 18注释部分)
    print("\n6. 运行映射 - Clusters模式...")
    
    try:
        ad_map_clusters = tg.map_cells_to_space(
            ad_sc, 
            ad_spatial,
            mode="clusters",
            cluster_label='anno_level_3_merged',
            density_prior='rna_count_based',
            num_epochs=500,
            device='cpu',
            verbose=True
        )
        
        print(f"   ✅ Clusters模式映射成功: {ad_map_clusters.shape}")
        
        # 手动添加细胞类型信息到obs中（确保细胞类型信息正确保存）
        if 'anno_level_3_merged' not in ad_map_clusters.obs.columns:
            print("   🔧 手动添加细胞类型信息...")
            # 从单细胞数据获取细胞类型列表（按字母顺序排序）
            sc_cell_types = sorted(ad_sc.obs['anno_level_3_merged'].unique().tolist())
            print(f"   ✓ 单细胞细胞类型数: {len(sc_cell_types)}")
            
            # 确保细胞类型数量与映射结果行数一致
            if len(sc_cell_types) == ad_map_clusters.shape[0]:
                ad_map_clusters.obs['anno_level_3_merged'] = sc_cell_types
                print(f"   ✅ 细胞类型信息已添加到obs中")
            else:
                print(f"   ⚠️ 细胞类型数量不匹配: 单细胞{len(sc_cell_types)} vs 映射结果{ad_map_clusters.shape[0]}")
        else:
            print("   ✅ 细胞类型信息已存在")
        
        # 保存结果
        clusters_file = os.path.join(output_dir, "tangram_clusters_result.h5ad")
        ad_map_clusters.write_h5ad(clusters_file)
        print(f"   ✓ Clusters模式结果已保存: {clusters_file}")
        
    except Exception as e:
        print(f"   ❌ Clusters模式映射失败: {e}")
        ad_map_clusters = None
    
    # 7. 结果分析 (对应官方notebook的Cell 19-28)
    print("\n7. 结果分析...")
    
    results = {}
    
    if ad_map_cells is not None:
        print(f"   ✓ Cells模式结果分析:")
        print(f"     形状: {ad_map_cells.shape}")
        
        # 安全地获取训练基因数
        train_genes_count = len(ad_map_cells.uns.get('training_genes', []))
        print(f"     训练基因数: {train_genes_count}")
        
        # 分析训练分数 (对应官方notebook的Cell 25)
        if 'train_genes_df' in ad_map_cells.uns:
            train_scores = ad_map_cells.uns['train_genes_df']['train_score']
            print(f"     平均训练分数: {train_scores.mean():.4f}")
            print(f"     训练分数范围: {train_scores.min():.4f} - {train_scores.max():.4f}")
        
        results['cells'] = {
            'shape': ad_map_cells.shape,
            'train_genes': train_genes_count
        }
    
    if ad_map_clusters is not None:
        print(f"   ✓ Clusters模式结果分析:")
        print(f"     形状: {ad_map_clusters.shape}")
        
        # 安全地获取训练基因数
        train_genes_count = len(ad_map_clusters.uns.get('training_genes', []))
        print(f"     训练基因数: {train_genes_count}")
        
        # 分析训练分数
        if 'train_genes_df' in ad_map_clusters.uns:
            train_scores = ad_map_clusters.uns['train_genes_df']['train_score']
            print(f"     平均训练分数: {train_scores.mean():.4f}")
            print(f"     训练分数范围: {train_scores.min():.4f} - {train_scores.max():.4f}")
        
        results['clusters'] = {
            'shape': ad_map_clusters.shape,
            'train_genes': train_genes_count
        }
    
    # 8. 生成可视化 (对应官方notebook的Cell 25)
    print("\n8. 生成可视化...")
    
    try:
        if ad_map_cells is not None:
            # 训练分数可视化
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            
            if 'train_genes_df' in ad_map_cells.uns:
                train_scores = ad_map_cells.uns['train_genes_df']['train_score']
                ax.hist(train_scores, bins=20, alpha=0.7, color='skyblue')
                ax.set_title('Training Scores Distribution - Cells Mode')
                ax.set_xlabel('Training Score')
                ax.set_ylabel('Frequency')
                ax.axvline(train_scores.mean(), color='red', linestyle='--', 
                           label=f'Mean: {train_scores.mean():.4f}')
                ax.legend()
            
            plt.tight_layout()
            
            # 保存可视化
            plot_file = os.path.join(output_dir, "training_scores_cells.png")
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   ✓ Cells模式训练分数可视化已保存: {plot_file}")
        
        if ad_map_clusters is not None:
            # 训练分数可视化
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            
            if 'train_genes_df' in ad_map_clusters.uns:
                train_scores = ad_map_clusters.uns['train_genes_df']['train_score']
                ax.hist(train_scores, bins=20, alpha=0.7, color='lightgreen')
                ax.set_title('Training Scores Distribution - Clusters Mode')
                ax.set_xlabel('Training Score')
                ax.set_ylabel('Frequency')
                ax.axvline(train_scores.mean(), color='red', linestyle='--', 
                           label=f'Mean: {train_scores.mean():.4f}')
                ax.legend()
            
            plt.tight_layout()
            
            # 保存可视化
            plot_file = os.path.join(output_dir, "training_scores_clusters.png")
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   ✓ Clusters模式训练分数可视化已保存: {plot_file}")
        
    except Exception as e:
        print(f"   ⚠ 可视化生成失败: {e}")
    
    # 9. 生成报告
    print("\n9. 生成报告...")
    
    report_file = os.path.join(output_dir, "official_tangram_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Tangram映射结果报告 (官方流程)\n")
        f.write("=" * 40 + "\n\n")
        
        f.write("项目概述:\n")
        f.write("按照官方Tangram notebook流程运行映射。\n\n")
        
        f.write("数据信息:\n")
        f.write(f"  单细胞数据: {ad_sc.shape}\n")
        f.write(f"  空间数据: {ad_spatial.shape}\n")
        f.write(f"  训练基因数: {len(training_genes)}\n")
        f.write(f"  细胞类型数: {ad_sc.obs['anno_level_3_merged'].nunique()}\n")
        f.write(f"  单细胞数据路径: {sc_file}\n")
        f.write(f"  空间数据路径: {spatial_file}\n")
        f.write(f"  使用Baysor生成的空间细胞数据\n")
        f.write(f"  基因匹配: 大小写不敏感\n\n")
        
        f.write("流程步骤:\n")
        f.write("1. 数据加载\n")
        f.write("2. 数据预处理 (标准化)\n")
        f.write("3. 训练基因选择\n")
        f.write("4. 数据准备 (pp_adatas)\n")
        f.write("5. 映射 (map_cells_to_space)\n")
        f.write("6. 结果分析\n\n")
        
        if 'cells' in results:
            f.write(f"Cells模式结果:\n")
            f.write(f"  形状: {results['cells']['shape']}\n")
            f.write(f"  训练基因数: {results['cells']['train_genes']}\n\n")
        
        if 'clusters' in results:
            f.write(f"Clusters模式结果:\n")
            f.write(f"  形状: {results['clusters']['shape']}\n")
            f.write(f"  训练基因数: {results['clusters']['train_genes']}\n\n")
        
        f.write("与官方notebook的一致性:\n")
        f.write("1. 使用相同的预处理步骤\n")
        f.write("2. 使用pp_adatas准备数据\n")
        f.write("3. 使用相同的映射参数\n")
        f.write("4. 使用相同的分析流程\n")
        f.write("5. 生成相同的可视化\n\n")
        
        f.write("文件输出:\n")
        if 'cells' in results:
            f.write(f"  Cells模式结果: {cells_file}\n")
        if 'clusters' in results:
            f.write(f"  Clusters模式结果: {clusters_file}\n")
        f.write(f"  综合报告: {report_file}\n")
    
    print(f"   ✓ 综合报告已保存: {report_file}")
    
    print("\n" + "=" * 80)
    print("🎯 官方流程Tangram映射完成!")
    print("=" * 80)
    
    return results

# 如果在Jupyter notebook中运行
if __name__ == "__main__":
    results = run_tangram_official_flow()
    if results:
        print(f"\n📊 映射结果:")
        for mode, result in results.items():
            print(f"   {mode.upper()}模式: {result['shape']}, 训练基因: {result['train_genes']}")
    else:
        print("\n❌ 映射失败")
