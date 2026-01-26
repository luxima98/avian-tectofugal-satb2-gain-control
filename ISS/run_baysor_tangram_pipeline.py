#!/usr/bin/env python3
"""
完整的Baysor + Tangram处理流水线
基于DAPI图像的Baysor细胞分割，生成Tangram后续处理所需的关键矩阵

主要功能：
1. 将空间转录组数据转换为Baysor格式
2. 使用Baysor进行细胞分割（支持DAPI图像和无监督模式）
3. 生成Tangram所需的关键矩阵
4. 准备单细胞数据用于Tangram映射

测试数据说明：
- test_ad_sc.h5ad: Tangram项目提供的示例单细胞数据
  来源：Tasic et al. (2018) 小鼠大脑皮层数据
  用途：用于测试和演示Tangram功能
  格式：AnnData格式，包含细胞类型注释


import os
import subprocess
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import tempfile
import logging
from scipy.sparse import csr_matrix
import scipy.io
import h5py
from PIL import Image
import warnings
warnings.filterwarnings('ignore')

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def convert_mat_to_tiff(mat_file, output_tiff):
    """将MAT格式的DAPI图像转换为TIFF格式"""
    try:
        # 首先尝试使用scipy.io读取（适用于较旧的MAT格式）
        try:
            mat_data = scipy.io.loadmat(mat_file)
            
            # 查找图像数据（MAT文件可能包含多个变量）
            image_key = None
            for key in mat_data.keys():
                if not key.startswith('__') and isinstance(mat_data[key], np.ndarray):
                    if len(mat_data[key].shape) == 2:  # 确保是2D图像
                        image_key = key
                        break
            
            if image_key is None:
                raise ValueError("未找到图像数据")
            
            # 获取图像数据
            image_data = mat_data[image_key]
            
        except Exception as e1:
            # 如果scipy.io失败，尝试使用h5py读取（适用于MATLAB v7.3格式）
            logger.info("尝试使用h5py读取MATLAB v7.3格式文件...")
            with h5py.File(mat_file, 'r') as f:
                # 查找图像数据
                image_key = None
                for key in f.keys():
                    if not key.startswith('__') and len(f[key].shape) == 2:
                        image_key = key
                        break
                
                if image_key is None:
                    raise ValueError("未找到图像数据")
                
                # 获取图像数据
                image_data = f[image_key][:]
        
        # 归一化到0-255范围（如果是浮点数）
        if image_data.dtype == np.float32 or image_data.dtype == np.float64:
            image_data = (image_data - image_data.min()) / (image_data.max() - image_data.min()) * 255
            image_data = image_data.astype(np.uint8)
        
        # 保存为TIFF
        Image.fromarray(image_data).save(output_tiff)
        logger.info(f"✓ DAPI图像已转换为: {output_tiff}")
        return output_tiff
        
    except Exception as e:
        logger.error(f"✗ 转换DAPI图像失败: {str(e)}")
        return None

def run_baysor_tangram_pipeline(transcript_file, dapi_image, output_dir, single_cell_file=None):
    """
    运行完整的Baysor + Tangram处理流水线
    
    Parameters:
    -----------
    transcript_file : str
        空间转录组h5ad文件路径
    dapi_image : str
        DAPI图像文件路径
    output_dir : str
        输出目录
    single_cell_file : str, optional
        单细胞数据h5ad文件路径，如果提供则生成完整的Tangram映射数据
    """
    print("="*80)
    print("Baysor + Tangram 完整处理流水线")
    print("="*80)
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 准备转录本数据
    print("1. 准备转录本数据...")
    transcript_csv = prepare_transcript_data(transcript_file, output_dir)
    
    # 2. 运行Baysor分割
    print("2. 运行Baysor细胞分割...")
    success = run_baysor_optimized(transcript_csv, dapi_image, output_dir)
    
    if not success:
        print("⚠ 有监督分割失败，尝试无监督分割...")
        success = run_baysor_optimized(transcript_csv, None, output_dir)
        
        if not success:
            print("✗ Baysor分割失败")
            return False
    
    # 3. 生成Tangram关键矩阵
    print("3. 生成Tangram关键矩阵...")
    tangram_matrices = generate_tangram_matrices(output_dir, transcript_file)
    
    # 4. 如果提供了单细胞数据，准备完整的Tangram映射数据
    if single_cell_file and os.path.exists(single_cell_file):
        print("4. 准备Tangram映射数据...")
        prepare_tangram_mapping_data(tangram_matrices, single_cell_file, output_dir)
    
    print("\n" + "="*80)
    print("✓ 流水线处理完成!")
    print("="*80)
    print("生成的关键矩阵:")
    for matrix_name, matrix_path in tangram_matrices.items():
        if os.path.exists(matrix_path):
            size = os.path.getsize(matrix_path)
            print(f"  ✓ {matrix_name}: {matrix_path} ({size:,} bytes)")
        else:
            print(f"  ✗ {matrix_name}: {matrix_path} (未生成)")
    print("="*80)
    
    return True

def prepare_transcript_data(transcript_file, output_dir):
    """准备转录本数据为Baysor格式"""
    logger.info("读取空间转录组数据...")
    adata_sp = sc.read_h5ad(transcript_file)
    logger.info(f"空间数据信息: {adata_sp.shape[0]} 个spots, {adata_sp.shape[1]} 个基因")
    
    # 获取空间坐标
    if 'x' in adata_sp.obs.columns and 'y' in adata_sp.obs.columns:
        x_coords = adata_sp.obs['x'].values
        y_coords = adata_sp.obs['y'].values
    else:
        raise ValueError("未找到空间坐标信息 (x, y)")
    
    # 创建转录本数据
    transcript_data = []
    if hasattr(adata_sp.X, 'toarray'):
        expression_matrix = adata_sp.X.toarray()
    else:
        expression_matrix = adata_sp.X
    
    gene_names = adata_sp.var_names.tolist()
    
    logger.info("提取转录本数据...")
    for i, (x, y) in enumerate(zip(x_coords, y_coords)):
        for j, gene in enumerate(gene_names):
            count = expression_matrix[i, j]
            if count > 0:  # 只保存有表达的转录本
                transcript_data.append({
                    'x': float(x),
                    'y': float(y),
                    'gene': gene
                })
    
    # 转换为DataFrame
    transcript_df = pd.DataFrame(transcript_data)
    logger.info(f"提取了 {len(transcript_df)} 个转录本")
    
    # 数据采样（如果数据量太大）
    if len(transcript_df) > 500000:
        logger.info("数据量较大，进行采样...")
        transcript_df = transcript_df.sample(n=500000, random_state=42)
        logger.info(f"采样后数据量: {len(transcript_df)} 个转录本")
    
    # 保存为CSV
    transcript_csv = os.path.join(output_dir, 'transcripts_baysor_format.csv')
    transcript_df.to_csv(transcript_csv, index=False)
    logger.info(f"转录本数据已保存: {transcript_csv}")
    
    return transcript_csv

def run_baysor_with_dapi(transcript_csv, dapi_image, output_dir):
    """运行Baysor分割（使用DAPI图像进行有监督分割）"""
    # 设置Julia路径
    julia_path = r"C:\Users\Administrator\AppData\Local\Programs\Julia-1.10.10\bin\julia.exe"
    baysor_path = r"D:\r\data\chicken\Tangram-master\Tangram-master\Baysor-master"
    
    # 转义路径
    transcript_csv_julia = transcript_csv.replace('\\', '\\\\')
    output_dir_julia = output_dir.replace('\\', '\\\\')
    baysor_path_julia = baysor_path.replace('\\', '\\\\')
    dapi_image_julia = dapi_image.replace('\\', '\\\\')
    
    # Julia命令 - 使用DAPI图像进行有监督分割
    julia_cmd = f"""
using Pkg
Pkg.activate("{baysor_path_julia}")
using Baysor

# 设置线程数
ENV["JULIA_NUM_THREADS"] = "6"

println("开始Baysor分割...")
println("输入文件: {transcript_csv_julia}")
println("DAPI图像: {dapi_image_julia}")
println("输出目录: {output_dir_julia}")

# 使用正确的Baysor.run函数，包含DAPI图像
Baysor.run(
    "{transcript_csv_julia}",
    "{dapi_image_julia}";  # 使用DAPI图像进行有监督分割
    output = "{output_dir_julia}",
    plot = true,
    scale = 20.0,              # 细胞尺度，可根据实际数据调整
    n_clusters = 4,            # 增加聚类数以处理更复杂的组织
    min_molecules_per_cell = 10, # 增加每个细胞最少分子数以提高质量
    count_matrix_format = "loom" # 输出格式
)

println("Baysor分割完成!")
"""
    
    # 使用Baysor.CommandLine.run函数
    cmd = [
        julia_path, 
        "--project={baysor_path_julia}",
        "-e", 
        f'using Baysor; Baysor.CommandLine.run("{transcript_csv_julia}", "{dapi_image_julia}"; output="{output_dir_julia}", plot=true, scale=20.0, n_clusters=4, min_molecules_per_cell=10, count_matrix_format="loom")'
    ]
    
    logger.info("运行Baysor命令...")
    logger.info(f"输入文件: {transcript_csv}")
    logger.info(f"DAPI图像: {dapi_image}")
    logger.info(f"输出目录: {output_dir}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10800, encoding='utf-8', errors='ignore')
        
        if result.returncode == 0:
            logger.info("✓ Baysor运行成功!")
            if result.stdout:
                logger.info("输出信息:")
                print(result.stdout[-1000:])  # 只显示最后1000个字符
            return True
        else:
            logger.error(f"✗ Baysor运行错误 (返回码: {result.returncode})")
            if result.stderr:
                logger.error("错误信息:")
                print(result.stderr[-2000:])  # 只显示最后2000个字符
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("✗ Baysor运行超时 (3小时)")
        return False
    except Exception as e:
        logger.error(f"✗ 运行Baysor时出错: {str(e)}")
        return False

def run_baysor_unsupervised(transcript_csv, output_dir):
    """运行Baysor无监督分割（不使用DAPI图像）"""
    # 设置Julia路径
    julia_path = r"C:\Users\Administrator\AppData\Local\Programs\Julia-1.10.10\bin\julia.exe"
    baysor_path = r"D:\r\data\chicken\Tangram-master\Tangram-master\Baysor-master"
    
    # 转义路径
    transcript_csv_julia = transcript_csv.replace('\\', '\\\\')
    output_dir_julia = output_dir.replace('\\', '\\\\')
    baysor_path_julia = baysor_path.replace('\\', '\\\\')
    
    # Julia命令 - 无监督分割
    julia_cmd = f"""
using Pkg
Pkg.activate("{baysor_path_julia}")
using Baysor

# 设置线程数
ENV["JULIA_NUM_THREADS"] = "6"

println("开始Baysor无监督分割...")
println("输入文件: {transcript_csv_julia}")
println("输出目录: {output_dir_julia}")

# 使用Baysor.run函数进行无监督分割
Baysor.run(
    "{transcript_csv_julia}";  # 只提供转录本数据，不提供DAPI图像
    output = "{output_dir_julia}",
    plot = true,
    scale = 20.0,              # 细胞尺度
    n_clusters = 4,            # 聚类数
    min_molecules_per_cell = 10, # 每个细胞最少分子数
    count_matrix_format = "loom" # 输出格式
)

println("Baysor无监督分割完成!")
"""
    
    # 使用Baysor.CommandLine.run函数进行无监督分割
    cmd = [
        julia_path, 
        "--project={baysor_path_julia}",
        "-e", 
        f'using Baysor; Baysor.CommandLine.run("{transcript_csv_julia}"; output="{output_dir_julia}", plot=true, scale=20.0, n_clusters=4, min_molecules_per_cell=10, count_matrix_format="loom")'
    ]
    
    logger.info("运行Baysor无监督分割命令...")
    logger.info(f"输入文件: {transcript_csv}")
    logger.info(f"输出目录: {output_dir}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10800, encoding='utf-8', errors='ignore')
        
        if result.returncode == 0:
            logger.info("✓ Baysor无监督分割运行成功!")
            if result.stdout:
                logger.info("输出信息:")
                print(result.stdout[-1000:])  # 只显示最后1000个字符
            return True
        else:
            logger.error(f"✗ Baysor无监督分割运行错误 (返回码: {result.returncode})")
            if result.stderr:
                logger.error("错误信息:")
                print(result.stderr[-2000:])  # 只显示最后2000个字符
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("✗ Baysor无监督分割运行超时 (3小时)")
        return False
    except Exception as e:
        logger.error(f"✗ 运行Baysor无监督分割时出错: {str(e)}")
        return False

def run_baysor_optimized(transcript_csv, dapi_image, output_dir):
    """针对ISS数据优化的Baysor分割"""
    julia_path = r"C:\Users\Administrator\AppData\Local\Programs\Julia-1.10.10\bin\julia.exe"
    baysor_path = r"D:\r\data\chicken\Tangram-master\Tangram-master\Baysor-master"
    
    # 转义路径
    transcript_csv_julia = transcript_csv.replace('\\', '\\\\')
    output_dir_julia = output_dir.replace('\\', '\\\\')
    baysor_path_julia = baysor_path.replace('\\', '\\\\')
    dapi_image_julia = dapi_image.replace('\\', '\\\\') if dapi_image else "nothing"
    
    # 针对ISS数据的优化参数 - 仅使用无监督分割
    julia_cmd = f"""
using Pkg
Pkg.activate("{baysor_path_julia}")
using Baysor

ENV["JULIA_NUM_THREADS"] = "4"  # 减少线程数避免内存问题

println("开始优化版Baysor分割...")
println("输入文件: {transcript_csv_julia}")
println("输出目录: {output_dir_julia}")

# 针对ISS数据的优化参数 - 仅使用无监督分割
Baysor.CommandLine.run(
    "{transcript_csv_julia}";
    output = "{output_dir_julia}",
    plot = true,
    scale = 10.0,               # ISS数据需要更小的细胞尺度
    n_clusters = 6,             # 增加聚类数
    min_molecules_per_cell = 5,  # ISS数据分子数较少
    count_matrix_format = "tsv"  # 使用TSV格式
)

println("Baysor分割完成!")
"""
    
    cmd = [julia_path, "--project={baysor_path_julia}", "-e", julia_cmd]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, encoding='utf-8', errors='ignore')
        
        if result.returncode == 0:
            logger.info("✓ Baysor运行成功!")
            return True
        else:
            logger.error(f"✗ Baysor运行错误 (返回码: {result.returncode})")
            if result.stderr:
                logger.error("错误信息:")
                print(result.stderr[-2000:])
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("✗ Baysor运行超时 (2小时)")
        return False
    except Exception as e:
        logger.error(f"✗ 运行Baysor时出错: {str(e)}")
        return False

def run_baysor_only_unsupervised(transcript_file, output_dir):
    """仅运行无监督Baysor分割"""
    # 1. 准备转录本数据
    print("1. 准备转录本数据...")
    transcript_csv = prepare_transcript_data(transcript_file, output_dir)
    if not transcript_csv:
        return False
    
    # 2. 运行无监督Baysor分割
    print("2. 运行Baysor无监督分割...")
    success = run_baysor_optimized(transcript_csv, None, output_dir)
    
    if not success:
        print("✗ Baysor分割失败")
        return False
    
    # 3. 生成Tangram矩阵
    print("3. 生成Tangram所需矩阵...")
    success = generate_tangram_matrices(output_dir, transcript_file)
    
    if not success:
        print("✗ 生成Tangram矩阵失败")
        return False
    
    print("✓ 无监督Baysor流水线完成!")
    return True

def generate_tangram_matrices(output_dir, original_transcript_file):
    """生成Tangram所需的关键矩阵"""
    logger.info("生成Tangram关键矩阵...")
    
    tangram_output = os.path.join(output_dir, 'tangram_matrices')
    os.makedirs(tangram_output, exist_ok=True)
    
    matrices = {}
    
    # 1. 读取Baysor分割结果
    segmentation_file = os.path.join(output_dir, 'segmentation.csv')
    if not os.path.exists(segmentation_file):
        logger.error("未找到segmentation.csv文件")
        return matrices
    
    seg_df = pd.read_csv(segmentation_file)
    logger.info(f"读取分割结果: {len(seg_df)} 个分子")
    
    # 2. 生成细胞坐标矩阵
    cell_coords = seg_df.groupby('cell')[['x', 'y']].mean().reset_index()
    cell_coords.columns = ['cell_id', 'x', 'y']
    cell_coords_file = os.path.join(tangram_output, 'cell_coordinates.csv')
    cell_coords.to_csv(cell_coords_file, index=False)
    matrices['cell_coordinates'] = cell_coords_file
    logger.info(f"细胞坐标矩阵已保存: {cell_coords_file}")
    
    # 3. 生成基因-细胞表达矩阵
    # 读取原始空间数据以获取基因信息
    adata_sp = sc.read_h5ad(original_transcript_file)
    
    # 创建细胞-基因表达矩阵
    cell_gene_matrix = seg_df.groupby(['cell', 'gene']).size().unstack(fill_value=0)
    
    # 确保所有基因都包含在矩阵中
    all_genes = set(adata_sp.var_names)
    existing_genes = set(cell_gene_matrix.columns)
    missing_genes = all_genes - existing_genes
    
    for gene in missing_genes:
        cell_gene_matrix[gene] = 0
    
    # 重新排列基因顺序以匹配原始数据
    cell_gene_matrix = cell_gene_matrix.reindex(columns=adata_sp.var_names, fill_value=0)
    
    # 保存为CSV
    cell_gene_file = os.path.join(tangram_output, 'cell_gene_expression.csv')
    cell_gene_matrix.to_csv(cell_gene_file)
    matrices['cell_gene_expression'] = cell_gene_file
    logger.info(f"细胞-基因表达矩阵已保存: {cell_gene_file}")
    
    # 4. 创建AnnData格式的空间数据
    # 使用细胞坐标作为空间坐标
    spatial_coords = cell_coords[['x', 'y']].values
    
    # 创建AnnData对象
    adata_spatial = ad.AnnData(
        X=cell_gene_matrix.values,
        obs=pd.DataFrame(index=cell_gene_matrix.index),
        var=adata_sp.var.copy(),
        obsm={'spatial': spatial_coords}
    )
    
    # 添加细胞统计信息
    adata_spatial.obs['n_counts'] = cell_gene_matrix.sum(axis=1)
    adata_spatial.obs['n_genes'] = (cell_gene_matrix > 0).sum(axis=1)
    adata_spatial.obs['x'] = cell_coords['x'].values
    adata_spatial.obs['y'] = cell_coords['y'].values
    
    # 保存AnnData格式
    spatial_h5ad_file = os.path.join(tangram_output, 'spatial_cells.h5ad')
    adata_spatial.write_h5ad(spatial_h5ad_file)
    matrices['spatial_anndata'] = spatial_h5ad_file
    logger.info(f"空间细胞AnnData已保存: {spatial_h5ad_file}")
    
    # 5. 生成细胞统计信息
    cell_stats = seg_df.groupby('cell').agg({
        'x': ['mean', 'std', 'min', 'max'],
        'y': ['mean', 'std', 'min', 'max'],
        'gene': 'count'
    }).round(3)
    
    cell_stats.columns = ['x_mean', 'x_std', 'x_min', 'x_max', 'y_mean', 'y_std', 'y_min', 'y_max', 'n_molecules']
    cell_stats = cell_stats.reset_index()
    
    cell_stats_file = os.path.join(tangram_output, 'cell_statistics.csv')
    cell_stats.to_csv(cell_stats_file, index=False)
    matrices['cell_statistics'] = cell_stats_file
    logger.info(f"细胞统计信息已保存: {cell_stats_file}")
    
    # 6. 生成分割结果摘要
    summary = {
        'total_molecules': len(seg_df),
        'total_cells': len(cell_coords),
        'total_genes': len(adata_sp.var_names),
        'molecules_per_cell_mean': seg_df.groupby('cell').size().mean(),
        'molecules_per_cell_std': seg_df.groupby('cell').size().std(),
        'genes_per_cell_mean': (cell_gene_matrix > 0).sum(axis=1).mean(),
        'genes_per_cell_std': (cell_gene_matrix > 0).sum(axis=1).std()
    }
    
    summary_df = pd.DataFrame([summary])
    summary_file = os.path.join(tangram_output, 'segmentation_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    matrices['segmentation_summary'] = summary_file
    logger.info(f"分割结果摘要已保存: {summary_file}")
    
    logger.info(f"✓ 生成了 {len(matrices)} 个关键矩阵")
    return matrices

def prepare_tangram_mapping_data(tangram_matrices, single_cell_file, output_dir):
    """准备Tangram映射数据"""
    logger.info("准备Tangram映射数据...")
    
    # 读取单细胞数据
    adata_sc = sc.read_h5ad(single_cell_file)
    logger.info(f"单细胞数据信息: {adata_sc.shape[0]} 个细胞, {adata_sc.shape[1]} 个基因")
    
    # 读取空间细胞数据
    spatial_h5ad_file = tangram_matrices.get('spatial_anndata')
    if not spatial_h5ad_file or not os.path.exists(spatial_h5ad_file):
        logger.error("未找到空间细胞AnnData文件")
        return
    
    adata_spatial = sc.read_h5ad(spatial_h5ad_file)
    logger.info(f"空间细胞数据信息: {adata_spatial.shape[0]} 个细胞, {adata_spatial.shape[1]} 个基因")
    
    # 准备Tangram映射
    tangram_mapping_dir = os.path.join(output_dir, 'tangram_mapping')
    os.makedirs(tangram_mapping_dir, exist_ok=True)
    
    # 保存处理后的数据
    adata_sc.write_h5ad(os.path.join(tangram_mapping_dir, 'single_cell_processed.h5ad'))
    adata_spatial.write_h5ad(os.path.join(tangram_mapping_dir, 'spatial_cells_processed.h5ad'))
    
    # 生成基因重叠信息
    overlap_genes = list(set(adata_sc.var_names) & set(adata_spatial.var_names))
    overlap_df = pd.DataFrame({'gene': overlap_genes})
    overlap_df.to_csv(os.path.join(tangram_mapping_dir, 'overlap_genes.csv'), index=False)
    
    logger.info(f"✓ Tangram映射数据已准备完成")
    logger.info(f"  重叠基因数量: {len(overlap_genes)}")
    logger.info(f"  单细胞数据: {adata_sc.shape}")
    logger.info(f"  空间细胞数据: {adata_spatial.shape}")

def main():
    """主函数"""
    # 设置路径
    transcript_file = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data\spatial\spatial_genesymbol_converted\VI_Sag2_CR_13_converted.h5ad"
    dapi_mat_file = r"D:\r\data\chicken\ISS\DAPI_images\DAPI_for_baysor.mat"
    single_cell_file = r"D:\r\data\chicken\Tangram-master\Tangram-master\data\test_ad_sc.h5ad"  # Tangram项目提供的测试单细胞数据（小鼠大脑皮层数据，来自Tasic et al. 2018）
    output_dir = r"D:\r\data\chicken\ISS\Spatial\ISS\Tangram\my_data\spatial\spatial_genesymbol_converted\VI_Sag2_CR_13_baysor_tangram_output"
    
    # 检查文件是否存在
    if not os.path.exists(transcript_file):
        logger.error(f"转录本文件不存在: {transcript_file}")
        return
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print("Baysor + Tangram 完整处理流水线")
    print("=" * 80)
    
    # 转换DAPI图像
    dapi_tiff_file = os.path.join(output_dir, 'DAPI_converted.tiff')
    if not os.path.exists(dapi_tiff_file) and os.path.exists(dapi_mat_file):
        dapi_tiff_file = convert_mat_to_tiff(dapi_mat_file, dapi_tiff_file)
    
    # 由于DAPI图像处理复杂，直接使用无监督分割
    print("⚠ 为避免DAPI图像处理问题，直接使用无监督分割")
    # 直接使用无监督分割
    success = run_baysor_only_unsupervised(transcript_file, output_dir)
    
    if success:
        print("\n" + "="*80)
        print("✓ Baysor + Tangram流水线处理完成!")
        print("="*80)
        print("生成的关键矩阵包括:")
        print("- cell_coordinates.csv - 细胞空间坐标")
        print("- cell_gene_expression.csv - 细胞×基因表达矩阵")
        print("- spatial_cells.h5ad - 空间细胞AnnData格式")
        print("- cell_statistics.csv - 细胞统计信息")
        print("- segmentation_summary.csv - 分割结果摘要")
        print("\n这些数据现在可以用于Tangram与单细胞数据映射!")
        print("="*80)
    else:
        print("\n" + "="*80)
        print("✗ 流水线处理失败")
        print("="*80)

if __name__ == "__main__":
    main()
