import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines
import numpy as np
import os
from typing import Dict, Optional

# ----------------------------
# Global style (academic)
# ----------------------------
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
# 放大整体字体（学术风格）
sns.set_context("talk", font_scale=1.2)

# 可统一调整的字号
TITLE_SIZE = 18
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 14
LEGEND_FONT_SIZE = 13

# 默认单位映射（可按需修改/补充）。
# 对于 RF_MVL_E.xlsx 的参数若未知单位，先用 a.u. (arbitrary units)。
DEFAULT_UNITS_MAP: Dict[str, str] = {
    # 常见几何量（若存在）
    'xc': 'px',
    'yc': 'px',
    'a': 'px',
    'b': 'px',
    'theta_deg': 'deg',
}


def _pretty_label(name: str, unit_map: Dict[str, str]) -> str:
    """将列名转换为学术英文轴标签: name [unit]。若无单位则使用 a.u."""
    unit = unit_map.get(name, 'a.u.')
    return f"{name} [{unit}]"


def _apply_axis_styles(ax):
    ax.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(TICK_LABEL_SIZE)


def visualize_data(units_map: Optional[Dict[str, str]] = None):
    """
    加载 Excel（RF_MVL_E.xlsx）并创建：
    - 每列分布图（仅数值列）
    - 基础 pairplot（仅数值列）
    所有坐标轴放大学术字体，并附单位（学术英文）。
    """
    if units_map is None:
        units_map = DEFAULT_UNITS_MAP

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, 'RF_MVL_E.xlsx')

        df = pd.read_excel(file_path)
        print(f"成功读取: {file_path}")
        print("列名:", df.columns.tolist())
    except FileNotFoundError:
        print(f"错误: 未在 {file_path} 找到文件")
        return
    except Exception as e:
        print(f"读取 Excel 发生异常: {e}")
        return

    # 仅数值列
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # --- 可视化 1: 每列分布 ---
    print("正在生成各列分布图（数值列）...")
    for col in numeric_cols:
        plt.figure(figsize=(6.5, 4.8))
        sns.histplot(df[col], kde=True, color="#4a6fa5")
        xlabel = _pretty_label(col, units_map)
        plt.xlabel(xlabel, fontsize=AXIS_LABEL_SIZE)
        plt.ylabel('Count', fontsize=AXIS_LABEL_SIZE)
        _apply_axis_styles(plt.gca())
        plt.title(f'Distribution of {col}', fontsize=TITLE_SIZE)
        plt.tight_layout()
        plt.savefig(os.path.join(script_dir, f'dist_{col}.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # --- 可视化 2: 基础 Pair Plot ---
    print("正在生成参数两两关系图 (pairplot)...")
    df_plot = df[numeric_cols].copy()
    # 用带单位的列名仅用于可视化
    rename_map = {c: _pretty_label(c, units_map) for c in df_plot.columns}
    df_plot_renamed = df_plot.rename(columns=rename_map)

    pair_plot = sns.pairplot(
        df_plot_renamed,
        diag_kind='kde',
        plot_kws={'alpha': 0.6, 's': 35, 'edgecolor': 'k', 'linewidth': 0.3},
        diag_kws={'fill': True}
    )
    pair_plot.fig.suptitle('Pairwise Relationships of Parameters', y=1.02, fontsize=TITLE_SIZE)
    # 放大刻度
    for ax in pair_plot.axes.flatten():
        if ax is not None:
            _apply_axis_styles(ax)
    pair_plot.savefig(os.path.join(script_dir, 'parameter_pair_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n基础可视化已生成: {script_dir}")


def generate_mvl_e_comparison_all_parameters(
    units_map: Optional[Dict[str, str]] = None,
    hue_col: str = 'Region',
    palette: Optional[Dict[str, str]] = None,
):
    """
    生成 MVL vs E 的“所有参数两两比较”图（即学术风格 pairplot）。
    - 数据来源: RF_MVL_E.xlsx
    - 若存在 hue_col (默认 'Region')，则以其区分 E 与 MVL
    - 输出: mvl_e_comparison_all_parameters.png / .pdf
    - 坐标轴: 学术英文 + 单位
    """
    if units_map is None:
        units_map = DEFAULT_UNITS_MAP
    if palette is None:
        palette = {'E': '#2c4ca0', 'MVL': '#4d7e54'}

    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, 'RF_MVL_E.xlsx')

    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"读取 {file_path} 失败: {e}")
        return

    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    if not numeric_cols:
        print("未找到数值列，无法生成比较图。")
        return

    # 构造用于绘图的数据框（只对数值列重命名为 带单位 标签）
    df_features = df[numeric_cols].copy()
    rename_map = {c: _pretty_label(c, units_map) for c in df_features.columns}
    df_features = df_features.rename(columns=rename_map)

    # 合并 hue 列（若存在）
    hue_available = (hue_col in df.columns)
    if hue_available:
        hue_series = df[hue_col].astype(str).str.strip().str.upper()
        df_plot = pd.concat([df_features, hue_series.rename(hue_col)], axis=1)
    else:
        df_plot = df_features.copy()

    # 生成 pairplot
    pair_kws = dict(
        diag_kind='kde',
        plot_kws={'alpha': 0.65, 's': 28, 'edgecolor': 'k', 'linewidth': 0.25},
        diag_kws={'fill': True},
    )

    if hue_available:
        g = sns.pairplot(df_plot, hue=hue_col, hue_order=['E','MVL'], palette=palette, **pair_kws)
        # 调整图例样式与位置（放到图外右侧，避免遮挡）
        if g._legend is not None:
            for text in g._legend.texts:
                text.set_fontsize(LEGEND_FONT_SIZE)
            g._legend.set_title(hue_col)
            g._legend.get_title().set_fontsize(LEGEND_FONT_SIZE)
            g._legend.set_bbox_to_anchor((1.02, 0.5))
            g._legend.set_loc('center left')
    else:
        g = sns.pairplot(df_plot, **pair_kws)

    # 为右侧图例留白
    g.fig.subplots_adjust(right=0.85)
    g.fig.suptitle('MVL vs E: Pairwise Comparison of All Parameters', y=1.02, fontsize=TITLE_SIZE)

    # 放大刻度
    for ax in g.axes.flatten():
        if ax is not None:
            _apply_axis_styles(ax)

    out_png = os.path.join(script_dir, 'mvl_e_comparison_all_parameters.png')
    out_pdf = os.path.join(script_dir, 'mvl_e_comparison_all_parameters.pdf')
    g.savefig(out_png, dpi=300, bbox_inches='tight')
    g.savefig(out_pdf, bbox_inches='tight')
    plt.close()

    print(f"已保存: {out_png} / {out_pdf}")


def generate_ellipse_overlays(
    csv_path: str | None = None,
    save_prefix: str = 'ellipse_overlays',
    draw_ellipses: bool = True,
    axis_units: str = 'px'
):
    """
    基于 ellipse_metrics.csv 叠加绘制椭圆轮廓，并以不同形状标注不同区域(E/MVL)的中心点。
    - X/Y 轴明确标注为中心坐标: X center(xc), Y center(yc)
    - E 中心: 圆点(o)，MVL 中心: 方形(s)
    - 颜色: E 蓝色(#2c4ca0)，MVL 绿色(#4d7e54)
    - a/b 视为半轴，width=2a, height=2b，角度使用 theta_deg（度）
    - 学术英文标签与放大字体
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if csv_path is None:
        csv_path = os.path.join(script_dir, 'ellipse_metrics.csv')

    if not os.path.exists(csv_path):
        print(f"未找到数据文件: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    required_cols = {'a', 'b', 'theta_deg', 'xc', 'yc'}
    if not required_cols.issubset(df.columns):
        print(f"列缺失，需包含: {required_cols}，实际列: {list(df.columns)}")
        return

    # 颜色与形状映射
    color_map = {'E': '#2c4ca0', 'MVL': '#4d7e54'}
    marker_map = {'E': 'o', 'MVL': 's'}

    fig, ax = plt.subplots(figsize=(8.5, 7.2))

    # 绘制椭圆轮廓（可选）
    if draw_ellipses:
        for _, row in df.iterrows():
            region = row.get('Region', 'E')
            ec = color_map.get(region, '#666666')
            try:
                e = Ellipse((row['xc'], row['yc']), width=2*row['a'], height=2*row['b'], angle=row['theta_deg'],
                            edgecolor=ec, facecolor='none', lw=1.3, alpha=0.9)
                ax.add_patch(e)
            except Exception as ex:
                print(f"跳过异常椭圆: {ex}")
                continue

    # 绘制中心点（不同形状）
    region_series = df.get('Region')
    if region_series is not None:
        for region, sub in df.groupby('Region'):
            ax.scatter(sub['xc'], sub['yc'], s=60, marker=marker_map.get(region, 'o'),
                       facecolors=color_map.get(region, '#666666'), edgecolors='k', alpha=0.9,
                       label=f"{region} center")
    else:
        ax.scatter(df['xc'], df['yc'], s=60, marker='o', facecolors='#777777', edgecolors='k', alpha=0.9, label='Center')

    ax.set_xlabel(f'X center [{axis_units}]', fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(f'Y center [{axis_units}]', fontsize=AXIS_LABEL_SIZE)
    ax.set_title('Ellipse Overlays (from ellipse_metrics.csv)', fontsize=TITLE_SIZE)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    _apply_axis_styles(ax)

    # 图例：中心点形状
    handles = []
    if 'E' in color_map:
        e_center = mlines.Line2D([], [], color='k', marker=marker_map['E'], linestyle='None',
                                 markerfacecolor=color_map['E'], markersize=9, label='E center')
        handles.append(e_center)
    if 'MVL' in color_map:
        mvl_center = mlines.Line2D([], [], color='k', marker=marker_map['MVL'], linestyle='None',
                                   markerfacecolor=color_map['MVL'], markersize=9, label='MVL center')
        handles.append(mvl_center)

    # 图例放到图外右侧，避免遮挡
    leg = ax.legend(handles=handles, bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=True)
    if leg is not None:
        for text in leg.get_texts():
            text.set_fontsize(LEGEND_FONT_SIZE)

    # 为右侧图外图例留白
    fig.subplots_adjust(right=0.82)

    png_path = os.path.join(script_dir, f'{save_prefix}.png')
    pdf_path = os.path.join(script_dir, f'{save_prefix}.pdf')
    plt.tight_layout()
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)

    print(f"已保存: {png_path} / {pdf_path}")


def generate_axial_rose_plot(
    angle_col: Optional[str] = None,
    bins: int = 18,
    normalize: bool = True,
    zero_at: str = 'N',  # 'N'|'E'|'S'|'W'
    clockwise: bool = True,
    output_prefix: str = 'axial_rose_plot'
):
    """
    生成轴向玫瑰图（axial rose plot）。
    - 数据来源: RF_MVL_E.xlsx
    - 角度单位: deg（常见为 0-180 或 0-360）。轴向数据按 180° 等价处理。
    - 数字刻度字体放大一倍（相对于全局 TICK_LABEL_SIZE）。
    - 输出: axial_rose_plot.png / .pdf
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, 'RF_MVL_E.xlsx')

    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"读取 {file_path} 失败: {e}")
        return

    # 自动识别角度列
    col = None
    if angle_col and angle_col in df.columns:
        col = angle_col
    elif 'theta_deg' in df.columns:
        col = 'theta_deg'
    elif 'theta' in df.columns:
        col = 'theta'

    if col is None:
        print("未找到角度列（期望 'theta' 或 'theta_deg'），无法生成轴向玫瑰图。")
        return

    angles_deg = pd.to_numeric(df[col], errors='coerce').dropna().to_numpy()
    if angles_deg.size == 0:
        print("角度列为空或无法解析。")
        return

    # 将角度归一到 [0, 180)
    angles_axial = np.mod(angles_deg, 180.0)

    # 统计 0-180 的直方分布
    bin_edges_deg = np.linspace(0.0, 180.0, bins + 1)
    counts, _ = np.histogram(angles_axial, bins=bin_edges_deg)

    if normalize and counts.sum() > 0:
        counts = counts / counts.sum()

    # 复制到 180-360 以体现轴向对称
    counts_full = np.concatenate([counts, counts])
    bin_edges_full_deg = np.concatenate([bin_edges_deg, bin_edges_deg[1:] + 180.0])

    # 转为弧度
    bin_edges_full_rad = np.deg2rad(bin_edges_full_deg)
    bin_widths = np.diff(bin_edges_full_rad)
    bin_centers = bin_edges_full_rad[:-1]

    # 极坐标绘制
    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, polar=True)

    # 朝向与零角
    zero_map = {'N': 'N', 'E': 'E', 'S': 'S', 'W': 'W'}
    ax.set_theta_zero_location(zero_map.get(zero_at, 'N'))
    ax.set_theta_direction(-1 if clockwise else 1)

    bars = ax.bar(bin_centers, counts_full, width=bin_widths, bottom=0.0,
                  align='edge', edgecolor='k', linewidth=0.6, color='#4a6fa5', alpha=0.85)

    # 角度刻度：每 30°
    thetaticks_deg = np.arange(0, 360, 30)
    ax.set_thetagrids(thetaticks_deg, labels=[f"{d}°" for d in thetaticks_deg])

    # 径向刻度：自动，且显示百分比（若归一化）
    if normalize:
        rticks = ax.get_yticks()
        ax.set_yticks(rticks)
        ax.set_yticklabels([f"{int(r*100)}%" for r in rticks])
        rlabel = 'Relative frequency'
    else:
        rlabel = 'Frequency'

    # 标题
    title_txt = f"Axial Rose Plot of {col}"
    ax.set_title(title_txt, va='bottom', fontsize=TITLE_SIZE)

    # 数字刻度字体放大一倍
    big_size = TICK_LABEL_SIZE * 2
    ax.tick_params(axis='both', labelsize=big_size)
    # thetagrid labels字体
    for txt in ax.get_xticklabels():
        txt.set_fontsize(big_size)
    for txt in ax.get_yticklabels():
        txt.set_fontsize(big_size)

    # 径向轴标签
    ax.set_rlabel_position(225)  # 放到左下象限，避免遮挡
    ax.text(np.deg2rad(225), ax.get_rmax()*1.05, rlabel, fontsize=AXIS_LABEL_SIZE, rotation=45,
            ha='center', va='center')

    # 保存
    out_png = os.path.join(script_dir, f'{output_prefix}.png')
    out_pdf = os.path.join(script_dir, f'{output_prefix}.pdf')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)

    print(f"已保存: {out_png} / {out_pdf}")


if __name__ == '__main__':
    # 生成 Excel 的基础图（学术风格 + 单位）
    try:
        visualize_data()
    except Exception as e:
        print(f"基础可视化跳过: {e}")

    # 生成 MVL vs E 两两比较总图
    try:
        generate_mvl_e_comparison_all_parameters()
    except Exception as e:
        print(f"MVL/E 比较图跳过: {e}")

    # 生成椭圆叠加图（带明确坐标轴与不同形状中心）
    try:
        generate_ellipse_overlays()
    except Exception as e:
        print(f"椭圆叠加图跳过: {e}")

    # 生成轴向玫瑰图（角度列自动识别）
    try:
        generate_axial_rose_plot()
    except Exception as e:
        print(f"轴向玫瑰图跳过: {e}")
