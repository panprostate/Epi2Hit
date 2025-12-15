import numpy as np
import pandas as pd
import sys
import os
import warnings
import subprocess
import tempfile
import io
warnings.filterwarnings('ignore')
import pickle
import pyBigWig
import statistics
import statistics as st
import seaborn as sns
from pathlib import Path
sns.set_style("whitegrid")
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow, Arc
from collections import defaultdict
from sklearn.mixture import GaussianMixture
from statannot import add_stat_annotation
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import gaussian_kde
from joblib import Parallel, delayed
import plotnine as p9
from scipy import signal
import pybedtools

index_count = defaultdict(int)


def read_dataset(file, sep='\t', header=0, index_col=0, comment=None):
    """
    See pandas.read_csv documentation - this is just a wrapper
    :param file:
    :param sep:
    :param header:
    :param index_col:
    :param comment:
    :return:
    """
    return pd.read_csv(file, sep=sep, header=header, index_col=index_col,
                       na_values=['Na', 'NA', 'NAN'], comment=comment)


def write_dataset(data, path, sep='\t', header=True, index=True, compression=None):
    """
    See pandas.to_csv documentation - this is just a wrapper
    :param data: data to save
    :param path: str or Path object, path to save
    :param sep:
    :param header:
    :param index:
    :param compression:
    :return:
    """
    from pathlib import Path

    if isinstance(path, str):
        path = Path(path)

    output_dir = path.parents[0]
    output_dir.mkdir(parents=True, exist_ok=True)

    return data.to_csv(path, sep=sep, header=header, index=index, compression=compression)

def median_scale(data, clip=None, exclude=None, axis=0):
    
    """
    Scale using median and median absolute deviation (MAD) across an axis.

    :param data: pd.DataFrame or pd.Series, input data to be scaled
    :param clip: float, optional, symmetrically clips scaled data to this value
    :param exclude: pd.Series, optional, samples to exclude while calculating median and MAD
    :param axis: int, default=0, axis to apply scaling on (0 for columns, 1 for rows)
    :return: pd.DataFrame with scaled data
    """
    from statsmodels.robust.scale import mad

    if exclude is not None:
        data_filtered = data.reindex(data.index & exclude[~exclude].index)
    else:
        data_filtered = data

    median = 1. * data_filtered.median(axis=axis)

    if isinstance(data, pd.Series):
        madv = 1. * mad(data_filtered.dropna())
        c_data = data.sub(median).div(madv)
    else:
        inv_axis = (axis + 1) % 2
        madv = 1. * data_filtered.apply(lambda x: mad(x.dropna()), axis=axis)
        c_data = data.sub(median, axis=inv_axis).div(madv, axis=inv_axis)

    if clip is not None:
        return c_data.clip(-clip, clip)
    return c_data


def find_open_chromatin_probes(df_probes, df_atac):
    """
    Find probes in open chromatin regions.
    
    :param df_probes: DataFrame containing probe information with columns
                      'Chromosome', 'Start', 'End', 'Strand'.
    :param df_atac: DataFrame containing ATAC-seq open chromatin regions
                    with columns 'chrom', 'start', 'end'.
    :return: DataFrame of probes that overlap with regions of open chromatin,
             with columns renamed to ['chrom', 'start', 'end', 'strand']
             and index named 'id'.
    """
    matching_rows = []

    for _, probe in tqdm(df_probes.iterrows(),
                         total=df_probes.shape[0],
                         desc="Finding probes in open chromatin"):
        
        mask_exact_match = (
            (df_atac['chrom'] == probe['Chromosome']) &
            (df_atac['start'] == probe['Start']) &
            (df_atac['end'] == probe['End'])
        )
        
        mask_inside_region = (
            (df_atac['chrom'] == probe['Chromosome']) &
            (df_atac['start'] <= probe['Start']) &
            (df_atac['end']   >= probe['End'])
        )
        
        mask = mask_exact_match | mask_inside_region
        
        if df_atac[mask].shape[0] > 0:
            matching_rows.append(probe)

    df_matching_probes = pd.DataFrame(matching_rows)

    df_matching_probes = df_matching_probes.rename(
        columns={
            'Chromosome': 'chrom',
            'Start': 'start',
            'End': 'end',
            'Strand': 'strand'
        }
    )
    df_matching_probes.index.name = 'id'

    return df_matching_probes
    

def probes2genes(meth_df, gene_ref_path):
    """
    Annotate methylation probes with closest genes using bedtools closest.

    :param meth_df: pd.DataFrame
        Contains methylation probes with columns:
        ['chrom', 'start', 'end', 'id', 'score', 'strand'] + sample columns.
        Index is typically 'id' or something like that.
    :param gene_ref_path: str
        Path to gene reference BED file (e.g.
        'gencode.v19.protein_coding_gene.minimal_chr.bed.gz')
        with columns:
        #chr, gene_start, gene_end, gene_name, strand, ...
    :return: pd.DataFrame
        Probes annotated with closest gene and distance, with columns:
        ['chrom', 'start', 'end', 'id', 'score', 'strand',
         'gene_name', 'gene_distance', ...samples...]
        Index is 'un_id' = 'gene_name_id'.
    """

    # sort and select required columns for bedtools
    meth_df_copy = (
        meth_df
        .sort_values(['chrom', 'start'])
        .reset_index()[['chrom', 'start', 'end', 'id', 'score', 'strand']]
    )

    # read gene reference from the path passed in
    reference = read_dataset(gene_ref_path, index_col=False)

    # keep only chr1–22,X,Y style chromosomes and add a dummy score
    reference['score'] = 0
    reference = reference[reference['#chr'].str.contains(r'^(chr[0-9XY]+)$')]

    reference = reference.rename(
        columns={
            '#chr': 'chrom',
            'gene_start': 'start',
            'gene_end': 'end',
            'gene_name': 'gene'
        }
    )[['chrom', 'start', 'end', 'gene', 'score', 'strand']]

    reference = reference.sort_values(['chrom', 'start'])

    # temp BED files for bedtools
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_meth, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_ref:

        meth_bed_path = tmp_meth.name
        ref_bed_path = tmp_ref.name

        meth_df_copy.to_csv(tmp_meth, sep='\t', header=False, index=False)
        reference.to_csv(tmp_ref, sep='\t', header=False, index=False)

    bedtools_path = os.path.join(os.path.dirname(sys.executable), "bedtools")

    try:
        if not os.path.exists(bedtools_path):
            raise FileNotFoundError(
                f"bedtools not found at {bedtools_path}. Is it installed?"
            )

        cmd = [
            bedtools_path, "closest",
            "-a", meth_bed_path,
            "-b", ref_bed_path,
            "-t", "all",
            "-s",
            "-d",
        ]
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        out_buf = io.StringIO(result.stdout)
        annotated_metyl = pd.read_csv(out_buf, sep='\t', header=None)

        annotated_metyl = annotated_metyl.rename(
            columns={
                0:  'chrom',
                1:  'start',
                2:  'end',
                3:  'id',
                4:  'score',
                5:  'strand',
                6:  'gene_chrom',
                7:  'gene_start',
                8:  'gene_end',
                9:  'gene_name',
                10: 'gene_score',
                11: 'gene_strand',
                12: 'gene_distance',
            }
        )

        annotated_metyl = annotated_metyl[
            ['chrom', 'start', 'end', 'id', 'score',
             'strand', 'gene_name', 'gene_distance']
        ]

        annotated_metyl['un_id'] = (
            annotated_metyl['gene_name'] + '_' +
            annotated_metyl['id'].astype(str)
        )

        # merge back extra columns (e.g. samples) from original meth_df
        annotated_metyl = annotated_metyl.merge(
            meth_df.loc[:, ~meth_df.columns.isin(annotated_metyl.columns)].reset_index(),
            on='id',
            how='right'
        ).set_index('un_id')

        return annotated_metyl

    finally:
        for p in (meth_bed_path, ref_bed_path):
            try:
                os.remove(p)
            except OSError:
                pass


def plot_filters(df,
                 cols=['Stromal', 'Immune', 'epiCMIT.hypo'],
                 cut_offs=[[0.4, -0.4], [0.4, -0.4], [0.3, -0.3]],
                 xlabel='Correlation'):
    """
    Plot distribution of filtering columns with cutoffs.

    :param df: pd.DataFrame, methylation probes data with rows as probes and columns including patient and filter data
    :param cols: list, default=['Stromal', 'Immune', 'epiCMIT.hypo'], columns to plot for filtering
    :param cut_offs: list of lists, default=[[0.4, -0.4], [0.4, -0.4], [0.3, -0.3]], cutoffs for each filter; [upper, lower]
    :param xlabel: str, default='Correlation', label for the x-axis
    :return: None, displays the distribution plot with cutoffs
    """
    if len(cols) != len(cut_offs):
        raise ValueError("cols and cut_offs lists must be the same length")

    colors = {
        'Stromal': 'darkblue',
        'Immune': 'darkblue',
        'epiCMIT.hypo': 'brown'
    }

    for col, (cutoff1, cutoff2) in zip(cols, cut_offs):
        ax = sns.distplot(df[col], label=f'{col}: cutoff1={cutoff1}, cutoff2={cutoff2}')
        ymax = max([h.get_height() for h in ax.patches])
        
        ax.vlines(x=cutoff1, ymin=0, ymax=ymax, color=colors[col], linestyles='--', label=f'{col} Cutoff: {cutoff1}')
        ax.vlines(x=cutoff2, ymin=0, ymax=ymax, color=colors[col], linestyles='--', label=f'{col} Cutoff: {cutoff2}')
        
        ax.set_xlabel(xlabel)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.show()
    
def prelim_filters(df,
                   cols=['Stromal', 'Immune', 'epiCMIT.hypo'],
                   cut_offs=[[0.4, -0.4], [0.4, -0.4], [0.3, -0.3]]):
    
    """
    Apply upper and lower cutoffs for specified columns.

    :param df: pd.DataFrame, rows - methylation probes, columns include patient data and filtering columns
    :param cols: list, list of columns to filter by
    :param cut_offs: list of lists, each sublist contains [upper, lower] cutoffs for specific filters
    :return: pd.DataFrame with values falling within cutoffs for each column
    """
    if len(cols) != len(cut_offs):
        raise ValueError("cols and cut_offs lists must be the same length")

    conditions = []
    for c, (upper, lower) in zip(cols, cut_offs):
        if lower is not None:
            conditions.append((df[c] < upper) & (df[c] > lower))
        else:
            conditions.append(df[c] < upper)
    
    return df[np.all(conditions, axis=0)]


def remove_black_listed_samples(df,
                                cols=['chrom', 'start', 'end', 'id', 'score', 'strand', 'gene_name',
                                      'gene_distance', 'Stromal', 'Immune', 'epiCMIT.hyper', 'epiCMIT.hypo'],
                                black=False,
                                path_to_list='~/PPCG_meth_sample_white_list.tsv',
                                sep=',',
                                index_col=None):
    """
    Remove black-listed samples and drop NaN values.

    :param df: pd.DataFrame, rows - methylation probes, columns - [samples + {cols}]
    :param cols: list, list of columns different from sample columns
    :param black: bool, whether to exclude (black) or include (white) samples from the list
    :param path_to_list: str, path to black/white list
    :param sep: str, default=',', delimiter for list file
    :param index_col: int, default=None, column to use as index
    :return: pd.DataFrame, columns include Chromosome + Start + End + filtered sample names
    """
    
    merge_back = df.loc[:,df.columns.isin(cols)]
    add_cols = len(merge_back.columns)
    
    df = df.loc[:,~df.columns.isin(cols)]
    
    f=pd.read_csv(path_to_list).iloc[:,0].unique()
    
    if black:
        df1 = df.loc[:,~df.columns.isin(f)]
        print(len(df1.columns), 'samples left after excluding black-listed samples + {} additional columns'.format(add_cols))
    else:
        df1 = df.loc[:,df.columns.isin(f)]
        print(len(df1.columns), 'samples left after filtering with white list + {} additional columns'.format(add_cols))
        
    return pd.concat([merge_back, df1], axis=1)


def select_type_of_samples(df,
                           cols=['chrom', 'start', 'end', 'id', 'score', 'strand', 'gene_name',
                                      'gene_distance', 'Stromal', 'Immune', 'epiCMIT.hyper', 'epiCMIT.hypo'],
                           path_to_ann = '~/PPCG_annotation.w.20.ct.txt',
                           sep='\t',
                           index_col=0,
                           column='Meth_Sample_Type',
                           types=['TUM', 'MET']):
    
    """
    Select specific types of samples from methylation data.

    :param df: pd.DataFrame, rows - methylation probes, columns include samples and additional specified columns
    :param cols: list, list of additional columns apart from sample data
    :param path_to_ann: str, path to the annotation file
    :param sep: str, default='\t', separator for the annotation file
    :param index_col: int, default=0, column to use as index
    :param column: str, default='Meth_Sample_Type', column containing sample types
    :param types: list, default=['TUM', 'MET'], sample types to select
    :return: pd.DataFrame with filtered sample types
    """
    merge_back = df.loc[:,df.columns.isin(cols)]
    df = df.loc[:,~df.columns.isin(cols)]
    
    # index-sample names
    ann = pd.read_csv(path_to_ann, sep=sep, index_col=index_col)[column]
    samples = ann[ann.isin(types)].index
    
    df1 = df.loc[:,df.columns.isin(samples)]
    print(len(df1.columns), 'samples selected')
    
    return pd.concat([merge_back, df1], axis=1)


def pomerantz_state(
    me_filtered,
    reference_bed_path,
    output_png_path=None,
):
    """
    Annotate methylation data with Pomerantz chromatin states.

    :param me_filtered: pd.DataFrame
        Filtered methylation data with columns ['chrom', 'start', 'end']
        and index 'un_id' (or convertible to that).
    :param reference_bed_path: str
        Path to the Pomerantz reference BED file
        (columns: chrom, start, end, state).
    :param output_png_path: str or None
        Path to save barplot of chromatin state counts.
        If None, the plot is not saved (and not required).
    :return: pd.DataFrame
        me_filtered with added 'pomerantz' and 'pom_region' columns.
    """

    # prepare input BED from me_filtered
    me_pom = me_filtered.reset_index()[['chrom', 'start', 'end', 'un_id']]

    # read reference chromHMM / Pomerantz file from the *explicit* path
    po = pd.read_csv(
        reference_bed_path,
        sep='\t',
        header=None,
        usecols=[0, 1, 2, 3],
        names=['chrom', 'start', 'end', 'state'],
    )

    # temp BED files for bedtools
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_me, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_ref:

        me_bed_path = tmp_me.name
        ref_bed_path = tmp_ref.name

        me_pom.to_csv(tmp_me, sep='\t', header=False, index=False)
        po.to_csv(tmp_ref, sep='\t', header=False, index=False)

    bedtools_path = os.path.join(os.path.dirname(sys.executable), "bedtools")

    try:
        if not os.path.exists(bedtools_path):
            raise FileNotFoundError(
                f"bedtools not found at {bedtools_path}. "
                "Is it installed?"
            )

        cmd = [
            bedtools_path, "intersect",
            "-a", me_bed_path,
            "-b", ref_bed_path,
            "-loj",
        ]
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        out_buf = io.StringIO(result.stdout)
        intersection = pd.read_csv(out_buf, sep='\t', header=None)

        intersection.columns = [
            'chrom', 'start', 'end', 'un_id',
            'pom_chrom', 'pom_start', 'pom_end', 'pomerantz'
        ]

        intersection['pom_region'] = (
            intersection['pom_chrom']
            + ':' + intersection['pom_start'].astype(str)
            + '-' + intersection['pom_end'].astype(str)
        )

        intersection = intersection.set_index('un_id')[['pomerantz', 'pom_region']]

        # replace missing with a label
        intersection['pomerantz'] = intersection['pomerantz'].replace('.', 'Unknown_state')

        state_replacements = {
            'Active_prostate_lineage-specific_promoter': 'Act_PSP',
            'Bivalent_poised_promoter': 'Bival_P',
            'Active_non-prostate_lineage_promoter': 'Act_nPSP',
            'Active_prostate_lineage-specific_enhancer': 'Act_PSE',
            'Active_non-prostate_lineage_enhancer': 'Act_nPSE',
            'Primed_non-prostate_lineage-specific_enhancer': 'Primed_nPSE',
            'Primed_prostate_lineage_enhancer': 'Primed_PSE',
            'Bivalent_poised_enhancer': 'Bival_E',
            'Heterochromatin': 'Het',
            'Repressed_chromatin': 'Rep_chrom',
        }
        intersection['pomerantz'].replace(state_replacements, inplace=True)

        me_with_pom = pd.concat([me_filtered, intersection], axis=1)

        # optional plotting: only if user provided a path
        if output_png_path is not None:
            data_counts = me_with_pom['pomerantz'].value_counts()

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.bar(data_counts.index, data_counts.values, color='#5B5B5B')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_xticklabels(data_counts.index, rotation=90)
            ax.set_xlabel('Chromatin State')
            ax.set_ylabel('Count')

            fig.tight_layout()

            os.makedirs(os.path.dirname(output_png_path), exist_ok=True)
            fig.savefig(output_png_path, dpi=300)
            plt.close(fig)

        return me_with_pom

    finally:
        for p in (me_bed_path, ref_bed_path):
            try:
                os.remove(p)
            except OSError:
                pass



index_count = defaultdict(int)

def process_region(
    reg,
    merge_probes,
    me_filtered2,
    correlation_threshold=0.6,
    corr_plots_dir=None,
):
    """
    Process a single region and return merged or individual probe data.

    For regions that are merged (mean corr >= threshold), also save a
    correlation heatmap in corr_plots_dir / GENE_REGION.png
    and return its path in 'corr_plot_paths'.
    """
    global index_count

    results = {
        'indices': [],
        'chromosomes': [],
        'starts': [],
        'ends': [],
        'ids': [],
        'scores': [],
        'strands': [],
        'gene_names': [],
        'gene_distances': [],
        'mes': [],
        'pomerantz_list': [],
        'pom_region_list': [],
        'regions2_list': [],
        'corr_plot_paths': [],   # NEW: store path (or NaN) per resulting row
    }
    
    cols = ['chrom', 'start', 'end', 'id', 'score', 'strand',
            'gene_name', 'gene_distance', 'pomerantz', 'pom_region', 'regions2']
    
    try:
        region_mask = merge_probes['regions2'] == reg
        probes_to_merge = merge_probes.loc[region_mask].index

        correlation_matrix = (
            me_filtered2.loc[:, ~me_filtered2.columns.isin(cols)]
                        .loc[probes_to_merge]
                        .T
                        .corr()
        )
        
        if check_correlation_threshold(correlation_matrix, correlation_threshold):
            # ---- merged region ----
            index_base = (
                pd.Series(probes_to_merge[0])
                  .str.replace('_cg[0-9]+', '', regex=True)[0]
                + '_' +
                merge_probes.loc[region_mask, 'pomerantz'].values[0]
            )
            index_full = index_base if index_count[index_base] == 0 else index_base + str(index_count[index_base])
            index_count[index_base] += 1

            results['indices'].append(index_full)
            
            ch = me_filtered2.loc[probes_to_merge, 'chrom'].values[0]
            st = me_filtered2.loc[probes_to_merge, 'start'].min()
            en = me_filtered2.loc[probes_to_merge, 'end'].max()
            id_list = list(me_filtered2.loc[probes_to_merge, 'id'])
            score = me_filtered2.loc[probes_to_merge, 'score'].values[0]
            strand = me_filtered2.loc[probes_to_merge, 'strand'].values[0]
            gene_name = me_filtered2.loc[probes_to_merge, 'gene_name'].values[0]
            gene_distance = np.nan

            results['chromosomes'].append(ch)
            results['starts'].append(st)
            results['ends'].append(en)
            results['ids'].append(id_list)
            results['scores'].append(score)
            results['strands'].append(strand)
            results['gene_names'].append(gene_name)
            results['gene_distances'].append(gene_distance)
            
            me_i = me_filtered2.loc[probes_to_merge].iloc[:, 8:-3].mean(axis=0)
            results['mes'].append(me_i)
            
            pom_val = merge_probes.loc[region_mask, 'pomerantz'].values[0]
            pom_region = merge_probes.loc[region_mask, 'pom_region'].values[0]
            regions2_val = merge_probes.loc[region_mask, 'regions2'].values[0]

            results['pomerantz_list'].append(pom_val)
            results['pom_region_list'].append(pom_region)
            results['regions2_list'].append(regions2_val)

            # --------- save correlation heatmap for this merged region ----------
            if corr_plots_dir is not None:
                upper = correlation_matrix.where(
                    np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
                ).stack()
                upper = upper[upper.notna()]

                fig, ax = plt.subplots(figsize=(6, 5))
                sns.heatmap(
                    correlation_matrix,
                    vmin=-1, vmax=1,
                    cmap="vlag",
                    square=True,
                    cbar_kws={"shrink": 0.7},
                    ax=ax
                )
                ax.set_title(reg)
                fig.tight_layout()

                os.makedirs(corr_plots_dir, exist_ok=True)
                safe_name = reg.replace(':', '_').replace('/', '_')
                out_path = os.path.join(corr_plots_dir, f"{safe_name}.png")
                fig.savefig(out_path, dpi=300)
                plt.close(fig)

                results['corr_plot_paths'].append(out_path)
            else:
                # No directory provided → no file, but keep alignment
                results['corr_plot_paths'].append(np.nan)

        else:
            # ---- no merge: keep each probe separately; no individual plots ----
            for probe in probes_to_merge:
                cpg_id = probe.split('_')[1]
                index_base = (
                    pd.Series(probe)
                      .str.replace('_cg[0-9]+', '', regex=True)[0]
                    + '_' + cpg_id + '_' +
                    merge_probes.loc[probe, 'pomerantz']
                )
                index_full = index_base if index_count[index_base] == 0 else index_base + str(index_count[index_base])
                index_count[index_base] += 1

                results['indices'].append(index_full)

                ch = me_filtered2.loc[probe, 'chrom']
                st = me_filtered2.loc[probe, 'start']
                en = me_filtered2.loc[probe, 'end']
                id_single = me_filtered2.loc[probe, 'id']
                score = me_filtered2.loc[probe, 'score']
                strand = me_filtered2.loc[probe, 'strand']
                gene_name = me_filtered2.loc[probe, 'gene_name']
                gene_distance = me_filtered2.loc[probe, 'gene_distance']

                results['chromosomes'].append(ch)
                results['starts'].append(st)
                results['ends'].append(en)
                results['ids'].append(id_single)
                results['scores'].append(score)
                results['strands'].append(strand)
                results['gene_names'].append(gene_name)
                results['gene_distances'].append(gene_distance)

                me_i = me_filtered2.loc[probe].iloc[8:]
                results['mes'].append(me_i)
                
                results['pomerantz_list'].append(merge_probes.loc[probe, 'pomerantz'])
                results['pom_region_list'].append(merge_probes.loc[probe, 'pom_region'])
                results['regions2_list'].append(merge_probes.loc[probe, 'regions2'])

                # no per-probe plot
                results['corr_plot_paths'].append(np.nan)
                
    except Exception as e:
        print(f"Error processing region {reg}: {e}")

    return results



def check_correlation_threshold(correlation_matrix, threshold=0.6):
    """
    Check if the *mean* pairwise correlation in a matrix exceeds a given threshold.
    """
    upper_triangle = correlation_matrix.where(
        np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
    ).stack()

    upper_triangle = upper_triangle[upper_triangle.notna()]

    if upper_triangle.empty:
        return False

    mean_corr = upper_triangle.mean()
    return mean_corr >= threshold


def process_probes(
    me_filtered2,
    correlation_threshold=0.6,
    num_cores=4,
    corr_plots_dir=None,
):
    """
    Process and merge probes based on correlation thresholds using multiple cores.

    If corr_plots_dir is not None, correlation heatmaps for merged regions
    are saved there and the resulting DataFrame contains a 'corr_plot_path'
    column with the file path (NaN for unmerged/unique probes).
    """
    global index_count

    regions2 = me_filtered2.index.str.replace('_cg[0-9]+', '_', regex=True) + me_filtered2['pom_region']
    me_filtered2 = me_filtered2.copy()
    me_filtered2['regions2'] = regions2
    
    regions = me_filtered2['regions2'].value_counts()
    regions = regions[regions > 1].index

    merge_probes = me_filtered2[me_filtered2['regions2'].isin(regions)]
    
    unique_probes = me_filtered2[~me_filtered2.index.isin(merge_probes.index)].copy()
    unique_probes['pomerantz'] = unique_probes['pomerantz']
    # no correlation plot for unique probes → add column with NaN now
    unique_probes['corr_plot_path'] = np.nan
    unique_probes.index = unique_probes.index + '_' + unique_probes['pomerantz'].astype(str)
    
    index_count = defaultdict(int)
    
    results = Parallel(n_jobs=num_cores, backend="threading")(
        delayed(process_region)(
            reg,
            merge_probes,
            me_filtered2,
            correlation_threshold,
            corr_plots_dir,   # pass directory down
        )
        for reg in tqdm(merge_probes['regions2'].unique(),
                        desc="Merging probes within regions")
    )

    combined_results = defaultdict(list)
    for result in results:
        for key, value in result.items():
            combined_results[key].extend(value)
    
    combined_probes = pd.concat(combined_results['mes'], axis=1).T
    combined_probes.index = combined_results['indices']
    combined_probes['chrom'] = combined_results['chromosomes']
    combined_probes['start'] = combined_results['starts']
    combined_probes['end'] = combined_results['ends']
    combined_probes['id'] = combined_results['ids']
    combined_probes['score'] = combined_results['scores']
    combined_probes['strand'] = combined_results['strands']
    combined_probes['gene_name'] = combined_results['gene_names']
    combined_probes['gene_distance'] = combined_results['gene_distances']
    combined_probes['pomerantz'] = combined_results['pomerantz_list']
    combined_probes['pom_region'] = combined_results['pom_region_list']
    combined_probes['regions2'] = combined_results['regions2_list']
    combined_probes['corr_plot_path'] = combined_results['corr_plot_paths']

    # Ensure column order is consistent with me_filtered2 (+ corr_plot_path at the end if new)
    for col in me_filtered2.columns:
        if col not in combined_probes.columns:
            combined_probes[col] = np.nan

    # corr_plot_path may not exist in me_filtered2; keep it at the end
    cols_final = list(me_filtered2.columns)
    if 'corr_plot_path' not in cols_final:
        cols_final.append('corr_plot_path')

    combined_probes = combined_probes[cols_final]

    final_combined = pd.concat([unique_probes[cols_final], combined_probes])

    return final_combined


def resolve_duplicates(df):
    """
    Renames duplicate indices to ensure uniqueness.
    """
    df = df.copy()
    
    seen = defaultdict(int)
    new_index = []

    for idx in df.index:
        if seen[idx] == 0:
            new_index.append(idx)
        else:
            new_index.append(f"{idx}{seen[idx]}")
        seen[idx] += 1
    
    df.index = new_index
    return df


def annotate_probes_with_loops(
    probe_data,
    loop_bedpe_path,
    output_path=None,
):
    """
    Annotate probe data with loop information from a BEDPE loop file.

    :param probe_data: pd.DataFrame
        Contains probe data with 'chrom', 'start', 'end' and index = 'un_id'
        (or will be set to that name).
    :param loop_bedpe_path: str
        Path to loops BEDPE file (gzipped), readable via
        read_dataset(loop_bedpe_path, index_col=False), and containing columns:
        chrom1, start1, end1, chrom2, start2, end2, loop_id.
    :param output_path: str or None
        If given, save **all** probes (with loops as an annotation) to this path (TSV).
        If None, nothing is written to disk.
    :return: pd.DataFrame
        Same as input `probe_data` but with a new 'loops' column and index
        prefixed by 'Ls_' for probes that intersect at least one loop anchor.
    """

    probe_data = probe_data.copy()
    probe_data.index.name = 'un_id'

    loop_data = read_dataset(loop_bedpe_path, index_col=False)

    required_cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'loop_id']
    missing = [c for c in required_cols if c not in loop_data.columns]
    if missing:
        raise ValueError(f"Missing required columns in loop file: {missing}")

    loops1 = loop_data[['chrom1', 'start1', 'end1', 'loop_id']]
    loops2 = loop_data[['chrom2', 'start2', 'end2', 'loop_id']]

    probes_bed_df = probe_data.reset_index()[['chrom', 'start', 'end', 'un_id']]

    bedtools_path = os.path.join(os.path.dirname(sys.executable), "bedtools")
    if not os.path.exists(bedtools_path):
        raise FileNotFoundError(
            f"bedtools not found at {bedtools_path}. Is it installed?"
        )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_probes, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_loops1, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_loops2:

        probes_bed_path = tmp_probes.name
        loops1_bed_path = tmp_loops1.name
        loops2_bed_path = tmp_loops2.name

        probes_bed_df.to_csv(tmp_probes, sep='\t', header=False, index=False)
        loops1.to_csv(tmp_loops1, sep='\t', header=False, index=False)
        loops2.to_csv(tmp_loops2, sep='\t', header=False, index=False)

    try:
        cmd1 = [
            bedtools_path, "intersect",
            "-a", probes_bed_path,
            "-b", loops1_bed_path,
            "-wa", "-wb",
        ]
        result1 = subprocess.run(
            cmd1,
            check=True,
            capture_output=True,
            text=True,
        )

        if result1.stdout.strip():
            out_buf1 = io.StringIO(result1.stdout)
            intersection1 = pd.read_csv(out_buf1, sep='\t', header=None)
            intersection1.columns = [
                'chrom', 'start', 'end', 'un_id',
                'chrom1', 'start1', 'end1', 'loop_id1'
            ]
        else:
            intersection1 = pd.DataFrame(
                columns=['chrom', 'start', 'end', 'un_id',
                         'chrom1', 'start1', 'end1', 'loop_id1']
            )

        cmd2 = [
            bedtools_path, "intersect",
            "-a", probes_bed_path,
            "-b", loops2_bed_path,
            "-wa", "-wb",
        ]
        result2 = subprocess.run(
            cmd2,
            check=True,
            capture_output=True,
            text=True,
        )

        if result2.stdout.strip():
            out_buf2 = io.StringIO(result2.stdout)
            intersection2 = pd.read_csv(out_buf2, sep='\t', header=None)
            intersection2.columns = [
                'chrom', 'start', 'end', 'un_id',
                'chrom2', 'start2', 'end2', 'loop_id2'
            ]
        else:
            intersection2 = pd.DataFrame(
                columns=['chrom', 'start', 'end', 'un_id',
                         'chrom2', 'start2', 'end2', 'loop_id2']
            )

        combined_intersection = pd.concat(
            [
                intersection1[['un_id', 'loop_id1']],
                intersection2[['un_id', 'loop_id2']],
            ],
            axis=0,
            ignore_index=True,
        )

        if combined_intersection.empty:
            probe_data['loops'] = np.nan
        else:
            aggregated_loops = combined_intersection.groupby('un_id').agg(
                lambda x: list(set(x.dropna()))
            ).reset_index()

            loop_dict = aggregated_loops.set_index('un_id').to_dict('index')

            def get_loops(un_id):
                entry = loop_dict.get(un_id, {})
                l1 = entry.get('loop_id1', [])
                l2 = entry.get('loop_id2', [])
                if not isinstance(l1, list):
                    l1 = [] if pd.isna(l1) else [l1]
                if not isinstance(l2, list):
                    l2 = [] if pd.isna(l2) else [l2]
                return l1 + l2

            probe_data['loops'] = probe_data.index.map(get_loops)

            def normalize_loops(x):
                if isinstance(x, (list, np.ndarray, tuple)):
                    if len(x) == 0:
                        return np.nan
                    if len(x) == 1:
                        return x[0]
                    return list(x)
                return x if pd.notna(x) else np.nan

            probe_data['loops'] = probe_data['loops'].apply(normalize_loops)

        def has_loop(val):
            if isinstance(val, (list, np.ndarray, tuple)):
                return len(val) > 0
            return pd.notna(val)

        probe_data.index = [
            ('Ls_' + str(idx)) if has_loop(val) else str(idx)
            for idx, val in zip(probe_data.index, probe_data['loops'])
        ]

        if output_path is not None:
            probe_data.to_csv(output_path, sep='\t')

        return probe_data

    finally:
        for p in (probes_bed_path, loops1_bed_path, loops2_bed_path):
            try:
                os.remove(p)
            except OSError:
                pass



def slice_bigwig(input_path, output_path, chromosomes):
    """
    Slice a bigWig file to include only specific chromosomes.
    
    :param input_path: str, path to the original bigWig file
    :param output_path: str, path where the sliced bigWig file will be saved
    :param chromosomes: list of str, list of chromosome names to extract and include in the sliced bigWig file
    :return: None, the sliced bigWig file is saved at the specified output path
    """
    with pyBigWig.open(input_path) as bw:
        with pyBigWig.open(output_path, 'w') as new_bw:
            header = [(chrom, bw.chroms(chrom)) for chrom in chromosomes if chrom in bw.chroms()]
            new_bw.addHeader(header)
            for chrom in chromosomes:
                if chrom in bw.chroms():
                    intervals = bw.intervals(chrom)
                    for start, end, value in intervals:
                        new_bw.addEntries([chrom], [start], ends=[end], values=[value])
                else:
                    print(f"Chromosome {chrom} not found in the file {input_path}")
                    

def plot_genes_with_loops_bigwig(gene_data,
                                 loop_data,
                                 bigwig_files,
                                 gene1_name=None, gene2_name=None,
                                 bed_data=None):
    """
    Plot chromatin loops, gene structures, epigenetic marks, and BED data.
    
    :param gene_data: pandas.DataFrame, DataFrame containing gene annotations with columns ['gene', 'chrom', 'start', 'end', 'strand', 'type']
    :param loop_data: pandas.DataFrame, DataFrame containing chromatin loop data with columns ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    :param bigwig_files: list of str, list of three bigWig file paths in the order [H3K4me3, H3K27ac, CTCF]
    :param gene1_name: str, optional, name of the first gene to plot
    :param gene2_name: str, optional, name of the second gene to plot
    :param bed_data: pandas.DataFrame, optional, DataFrame containing BED file data with columns ['chrom', 'start', 'end'] to be plotted
    :return: None, displays the plot with chromatin loops, genes, epigenetic marks, and BED data
    """
    
    if bigwig_files is None or len(bigwig_files) != 3:
        print("Please provide a list of three bigWig file paths for H3k4me3, H3K27ac, and CTCF.")
        return

    bw_H3k4me3 = pyBigWig.open(bigwig_files[0])
    bw_H3K27ac = pyBigWig.open(bigwig_files[1])
    bw_CTCF = pyBigWig.open(bigwig_files[2])

    genes = []
    if gene1_name:
        gene1_filtered = gene_data[gene_data['gene'] == gene1_name]
        if not gene1_filtered.empty:
            genes.append((gene1_name, gene1_filtered))
        else:
            print(f"No data found for gene: {gene1_name}")

    if gene2_name:
        gene2_filtered = gene_data[gene_data['gene'] == gene2_name]
        if not gene2_filtered.empty:
            genes.append((gene2_name, gene2_filtered))
        else:
            print(f"No data found for gene: {gene2_name}")

    if not genes:
        print("No valid gene data provided.")
        return

    gene_chrom = genes[0][1]['chrom'].iloc[0]
    plot_start = min(g[1]['start'].min() for g in genes) - 10000
    plot_end = max(g[1]['end'].max() for g in genes) + 10000

    loops_filtered = loop_data[
        (loop_data['chrom1'] == gene_chrom) &
        (loop_data['chrom2'] == gene_chrom) &
        (loop_data['start1'] >= plot_start) & (loop_data['end1'] <= plot_end) &
        (loop_data['start2'] >= plot_start) & (loop_data['end2'] <= plot_end)
    ]

    fig, axs = plt.subplots(6, 1, figsize=(6, 3.5), sharex=True,
                            gridspec_kw={'height_ratios': [0.3, 0.1, 0.3, 0.3, 0.3, 0.3]})
    fig.subplots_adjust(hspace=0.005)

    for _, loop in loops_filtered.iterrows():
        start1, end1 = loop['start1'] / 1e6, loop['end1'] / 1e6
        start2, end2 = loop['start2'] / 1e6, loop['end2'] / 1e6
        mid1, mid2 = (start1 + end1) / 2, (start2 + end2) / 2
        width = abs(mid2 - mid1)

        axs[0].hlines(y=0.1, xmin=start1, xmax=end1, color="black", linewidth=1)
        axs[0].hlines(y=0.1, xmin=start2, xmax=end2, color="black", linewidth=1)
        axs[0].set_ylabel('chromatin\nloops', rotation=90)

        axs[0].add_patch(Arc(((mid1 + mid2) / 2, 0.1), width, 1.5, theta1=0, theta2=180, edgecolor='black', lw=1, alpha=0.7))

    axs[0].set_ylim(0, 1.5)
    axs[0].set_yticks([])

    if bed_data is not None:
        for _, row in bed_data.iterrows():
            chrom, start, end = row['chrom'], row['start'], row['end']
            if chrom == gene_chrom and plot_start <= start <= plot_end:
                axs[1].axvline(x=start / 1e6, color='red', linewidth=2)
        axs[1].set_ylim(0, 1)
        axs[1].set_yticks([])
        axs[1].set_ylabel('BED', rotation=90)

    y_positions = [0.7, 0.3]
    for idx, (gene_name, gene_filtered) in enumerate(genes):
        y_position = y_positions[idx]
        arrow_direction = 1 if gene_filtered['strand'].iloc[0] == '+' else -1
        arrow_interval = 150000

        axs[2].hlines(y=y_position, xmin=gene_filtered['start'].min() / 1e6, xmax=gene_filtered['end'].max() / 1e6, color='black', linewidth=1)
        for _, exon in gene_filtered[gene_filtered['type'] == 'exon'].iterrows():
            axs[2].add_patch(Rectangle((exon['start'] / 1e6, y_position - 0.05), (exon['end'] - exon['start']) / 1e6, 0.1, color='black'))

        for position in range(gene_filtered['start'].min(), gene_filtered['end'].max(), arrow_interval):
            axs[2].add_patch(FancyArrow(position / 1e6, y_position, arrow_direction * 1000 / 1e6, 0, width=0.02, head_width=0.08, head_length=1000 / 1e6, color='black'))
        axs[2].text((gene_filtered['start'].min() + gene_filtered['end'].max()) / 2 / 1e6, y_position - 0.15, gene_name, ha='center', va='top', fontsize=7, style='italic')

    axs[2].set_ylim(0, 1)
    axs[2].set_yticks([])
    axs[2].set_xlim(plot_start / 1e6, plot_end / 1e6)

    bigwig_bin_centers = np.linspace(plot_start, plot_end, plot_end - plot_start) / 1e6

    def plot_bigwig_data(bw, ax, color, label):
        values = bw.values(gene_chrom, plot_start, plot_end, numpy=True)
        values = np.nan_to_num(values)
        ax.fill_between(bigwig_bin_centers, values, color=color, alpha=1)
        ax.set_ylabel(label, rotation=90)
        ax.set_yticks([])
        ax.spines['right'].set_visible(True)

    plot_bigwig_data(bw_H3k4me3, axs[3], 'darkgreen', 'H3K4me3')
    
    plot_bigwig_data(bw_H3K27ac, axs[4], 'darkorange', 'H3K27ac')
    
    plot_bigwig_data(bw_CTCF, axs[5], 'brown', 'CTCF')

    axs[5].set_xlabel('Position (Mb)')

    plt.tight_layout()
    plt.show()

    bw_H3k4me3.close()
    bw_H3K27ac.close()
    bw_CTCF.close()
    
    

def read_exp_data(path_to_exp_ann = '~/ann_RNAseq_v1.0.tsv',
                  path_to_exp = '~/expressions.tsv',
                  sep='\t',
                  sep1='\t',
                  index_col=0,
                  index_col1=0,
                  s_type_column='RNA_Sample_Type',
                  types=['TUM', 'MET'],
                  exp_sample_names_col='PPCG_RNA_Assay_ID',
                  meth_sample_names_col='PPCG_Meth_Assay_ID'):
    
    """
    Read and match expression data to methylation sample names.

    :param path_to_exp_ann: str, path to RNAseq annotation file
    :param path_to_exp: str, path to RNAseq expression data file
    :param sep, sep1: str, delimiters for annotation and expression data files, respectively
    :param index_col, index_col1: int, columns to use as index for each file
    :param s_type_column: str, column with sample type in annotation
    :param types: list, sample types to include
    :param exp_sample_names_col: str, column name with RNAseq sample IDs
    :param meth_sample_names_col: str, column name with methylation sample IDs
    :return: pd.DataFrame with expression data matched to methylation samples
    """
    
    exp = pd.read_csv(path_to_exp, sep=sep, index_col=index_col)
    ann = pd.read_csv(path_to_exp_ann, sep=sep1, index_col=index_col1)
    ann = ann[ann[s_type_column].isin(types)]
    ann = ann[ann[exp_sample_names_col].isin(exp.columns)][[meth_sample_names_col,
                                                     exp_sample_names_col]].set_index(exp_sample_names_col)
    
    exp_meth = pd.concat([ann,
           exp.T
          ], axis=1).set_index(meth_sample_names_col).dropna()[
    pd.concat([ann,
           exp.T
          ], axis=1).set_index(meth_sample_names_col).dropna().index.notna()
    ].T
    exp_meth = exp_meth.loc[:,~exp_meth.columns.duplicated()]
    
    return exp_meth


def prepare_dict_for_linreg(exp_meth, me, median_scale, output_path):
    
    """
    Prepare a dictionary of gene annotations with methylation and expression data for linear regression.

    :param exp_meth: pd.DataFrame, combined expression and methylation data
    :param me: pd.DataFrame, methylation data
    :param median_scale: function, scaling function (e.g., median scaling)
    :param output_path: str, file path to save the dictionary as a pickle file
    :return: dict with methylation feature indices as keys and DataFrames with 'methylation' and 'expression' columns for each gene
    """
    warnings.filterwarnings('ignore')
    
    cols = ['chrom', 'start', 'end', 'id', 'score', 'strand', 'gene_name',
            'gene_distance', 'Stromal', 'Immune', 'epiCMIT.hyper', 'epiCMIT.hypo',
            'pomerantz', 'pom_region', 'regions2', 'loops']

    so = exp_meth.columns.intersection(me.loc[:, ~me.columns.isin(cols)].columns)
    exp = exp_meth.loc[:, so].loc[:, ~exp_meth.loc[:, so].columns.duplicated()]
    met = me.loc[:, ~me.columns.isin(cols)].loc[:, so]

    pattern_corrected = r'^(?:Ls_)?([A-Za-z0-9-]+)(?:_.*)?$'

    gene_ann_me = dict()

    # Loop through each methylation index, apply transformations, and store in dictionary
    for ind in tqdm(met.index, desc="Preparing methylation-expression pairs for linreg"):
        try:
            g_name = pd.Series(ind).str.replace(pattern_corrected, r'\1', regex=True)[0]
            e = exp.loc[g_name].astype(float)
            m = median_scale(met.loc[ind].astype(float))
            mer = pd.concat([m, e], axis=1).rename(columns={ind: 'methylation', g_name: 'expression'})
            mer = mer.dropna()
            gene_ann_me[ind] = mer
        except:
            pass
    
    # Save the dictionary as a pickle file
    with open(output_path, 'wb') as fp:
        pickle.dump(gene_ann_me, fp, protocol=pickle.HIGHEST_PROTOCOL)

    return gene_ann_me


def prepare_TCGA_dict(
    me_with_loops,
    exp_path,
    output_path,
):
    """
    Prepare a dictionary of gene annotations with methylation and TCGA expression data
    for linear regression.

    :param me_with_loops: pd.DataFrame
        Methylation data AFTER all processing (pomerantz, regions2, loops, etc.).
    :param exp_path: str
        Path to TCGA expression matrix.
    :param output_path: str
        Path to save the resulting dictionary as a pickle file.
    :return: dict
        Keys: methylation feature indices
        Values: DataFrame with columns ['methylation', 'expression'].
    """
    warnings.filterwarnings('ignore')

    cols = [
        'chrom', 'start', 'end', 'id', 'score', 'strand', 'gene_name',
        'gene_distance', 'Stromal', 'Immune', 'epiCMIT.hyper', 'epiCMIT.hypo',
        'pomerantz', 'pom_region', 'regions2', 'loops'
    ]

    exp_all = read_dataset(exp_path, index_col=0)

    me = me_with_loops.copy()
    sample_cols_me = me.loc[:, ~me.columns.isin(cols)].columns

    so = exp_all.columns.intersection(sample_cols_me)

    exp = exp_all.loc[:, so].loc[:, ~exp_all.loc[:, so].columns.duplicated()]
    met = me.loc[:, ~me.columns.isin(cols)].loc[:, so]

    pattern_corrected = r'^(?:Ls_)?([A-Za-z0-9-]+)(?:_.*)?$'

    gene_ann_me = {}

    for ind in tqdm(met.index, desc="Preparing methylation-expression pairs for linreg"):
        try:
            g_name = pd.Series(ind).str.replace(pattern_corrected, r'\1', regex=True)[0]

            e = exp.loc[g_name].astype(float)
            m = met.loc[ind].astype(float)

            mer = pd.concat([m, e], axis=1)
            mer.columns = ['methylation', 'expression']
            mer = mer.dropna()

            if not mer.empty:
                gene_ann_me[ind] = mer
        except Exception:
            continue

    with open(output_path, 'wb') as fp:
        pickle.dump(gene_ann_me, fp, protocol=pickle.HIGHEST_PROTOCOL)

    return gene_ann_me

    
def lin_reg(
    linreg_dict,
    scale_methylation=True,
    scale_expression=False,
    plot=True,
    top_n=20,
    forest_kwargs=None,
    show_plot=False,
    save_plot=True,
    plot_dir="forest_linreg",
):
    """
    Perform linear regression for methylation and expression data, with optional median scaling,
    then filter by adjusted p-value and effect size, and optionally make a forest plot.

    :param linreg_dict: dict
        Keys = methylation feature indices
        Values = DataFrames with columns ['methylation', 'expression'].
    :param scale_methylation: bool
        If True (default), apply median_scale to methylation before regression.
    :param scale_expression: bool
        If True, apply median_scale to expression before regression.
    :param plot: bool
        If True, create a forest plot of the top_n most negative regression coefficients.
    :param top_n: int
        Number of most negative associations to show in the forest plot (default 20).
    :param forest_kwargs: dict or None
        Extra keyword arguments passed to forest_plot (e.g. {'title': 'My title'}).
    :param show_plot: bool
        If True, show the forest plot (plt.show()).
    :param save_plot: bool
        If True (default), save the forest plot into plot_dir.
    :param plot_dir: str
        Directory where forest plot will be saved if save_plot is True.
        Default: "forest_linreg".
    :return: pandas.DataFrame
        Columns: coef, p_value, q95l, q95h, rsquared, eCpGs, adjusted_p_value, gene,
        filtered to adjusted_p_value < 0.05 and coef < -0.4.
    """

    pattern_corrected = r'^(?:Ls_)?([A-Za-z0-9-]+)(?:_.*)?$'

    full_genes_meth = pd.DataFrame()

    for g in tqdm(list(linreg_dict.keys()), desc="Running linear regression"):
        try:
            df = linreg_dict[g]

            X = df[['methylation']].astype(float)
            y = df[['expression']].astype(float)

            if scale_methylation:
                X = median_scale(X)
            if scale_expression:
                y = median_scale(y)

            model = sm.OLS(endog=y, exog=X, data=df)
            results = model.fit()

            r_squared = pd.DataFrame({'rsquared': [results.rsquared]}, index=['methylation'])

            gene_sum = pd.concat(
                [
                    pd.DataFrame(results.params[[0]]).rename(columns={0: 'coef'}),
                    pd.DataFrame(results.pvalues[[0]]).rename(columns={0: 'p_value'}),
                    results.conf_int().rename(columns={0: 'q95l', 1: 'q95h'}).iloc[[0]],
                    r_squared,
                ],
                axis=1,
            )

            full_genes_meth = pd.concat(
                [full_genes_meth, gene_sum.rename(index={'methylation': g})],
                axis=0,
            )
        except Exception:
            print(g)

    full_genes_meth = full_genes_meth.sort_values('coef')

    full_genes_meth['eCpGs'] = full_genes_meth.index

    full_genes_meth['adjusted_p_value'] = multipletests(
        full_genes_meth['p_value'],
        method='fdr_bh'
    )[1].astype(float)

    full_genes_meth['gene'] = full_genes_meth['eCpGs'].str.replace(
        pattern_corrected, r'\1', regex=True
    )

    full_genes_meth = full_genes_meth[
        (full_genes_meth['adjusted_p_value'] <= 0.05) &
        (full_genes_meth['coef'] < -0.4)
    ]

    if plot and not full_genes_meth.empty:
        forest_kwargs = forest_kwargs or {}

        plot_df = full_genes_meth.nsmallest(
            min(top_n, full_genes_meth.shape[0]), 'coef'
        ).copy()

        if 'title' not in forest_kwargs:
            forest_kwargs['title'] = f"Top {min(top_n, plot_df.shape[0])} negative associations"

        save_path = None
        if save_plot:
            os.makedirs(plot_dir, exist_ok=True)
            save_path = os.path.join(
                plot_dir,
                f"forest_linreg_top{min(top_n, plot_df.shape[0])}.png"
            )

        forest_plot(
            conclusion_metrics=plot_df,
            metrics='Coefficient',
            save_path=save_path,
            show=show_plot,
            **forest_kwargs
        )

    return full_genes_meth
    

def forest_plot(conclusion_metrics,
                title,
                colour_interval=False,
                figsize=(4, 5),
                ylimits=(-8, 1),
                feature_order=None,
                metrics='Coefficient',
                display_n=False,
                font_size=6,
                save_path=None,
                show=True):
    """
    param: conclusion_metrics, pd.DataFrame with the validation metrics
    param: title, str
    param: figsize, tuple
    param: ylimits, tuple, limits for the y axis, makes the plot more readable
    param: metrics, str, whether it's hazard or odds ratio or coefficient
    param: font_size, int
    param: save_path, str or None - if not None, save figure to this path
    param: show, bool - if True, display figure with plt.show()
    """
    if metrics == 'Coefficient':
        markers_type = 'prognosis'
    elif metrics in ['odds ratio', 'Response coefficient']:
        markers_type = 'prediction'
    else:
        raise Exception('Unknown metric type')
        
    q = conclusion_metrics.copy()
    
    if feature_order is not None:
        q = q.loc[feature_order]

    if feature_order is None:
        feature_order = q.nsmallest(q.shape[0], 'coef').gene

    q['p_value'] = 'P = ' + round(q.p_value, 3).astype(str)
    if display_n:
        q['p_value'] = q['p_value'].str.cat(q.n_per_marker, sep='  ')
    elif not display_n:
        pass
    else:
        raise ValueError('display_n argument must be boolean')
        
    if colour_interval:
        good_colour = '#000096'
        bad_colour = '#E61E1E'
        neutral_colour = '#000000'
    else:
        good_colour = '#000000'
        bad_colour = '#000000'
        neutral_colour = '#000000'
           
    p = (p9.ggplot(q)
         + p9.geom_pointrange(p9.aes(x='gene', y='coef', ymin='q95l', ymax='q95h'))
         + p9.coord_flip()
         + p9.scale_x_discrete(limits=feature_order)
         + p9.geom_hline(yintercept=0, color='#E61E1E', size=1, linetype="dashed")
         + p9.labs(x='gene', y=metrics)
         + p9.geom_text(p9.aes(x='gene', y='coef', label='p_value'),
                        size=font_size, nudge_x=0.4, nudge_y=0.1)
         + p9.ggtitle(title)
         + p9.ylim(ylimits)
         + p9.theme_bw()
         + p9.theme(
             text=p9.element_text(size=font_size),
             axis_title=p9.element_text(color='black', size=font_size),
             axis_text=p9.element_text(color='black', size=font_size),
             plot_title=p9.element_text(color='black', size=font_size)
         )
        )

    fig = p.draw()
    fig.set_size_inches(figsize[0], figsize[1])

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def slice_probes_with_loops_by_linreg(
    linreg_results,
    probes_path,
    output_path=None,
    join_linreg=True,
):
    """
    Read 'probes_with_loops' file and slice it based on linreg results.

    :param linreg_results: pd.DataFrame
        Output of lin_reg(), index = eCpG IDs.
    :param probes_path: str
        Path to the probes-with-loops file.
    :param output_path: str or None
        If not None, save sliced probes table (TSV) to this path.
    :param join_linreg: bool
        If True, join regression metrics onto the probes DataFrame.
    :return: pd.DataFrame
    """

    probes = pd.read_csv(probes_path, sep='\t', index_col=0)

    probes.index = probes.index.astype(str)
    linreg_results = linreg_results.copy()
    linreg_results.index = linreg_results.index.astype(str)

    common_ids = probes.index.intersection(linreg_results.index)

    sliced = probes.loc[common_ids].copy()

    if join_linreg:
        reg_cols = [c for c in ['coef', 'p_value', 'adjusted_p_value', 'q95l', 'q95h', 'rsquared', 'gene', 'eCpGs']
                    if c in linreg_results.columns]
        sliced = sliced.join(linreg_results[reg_cols], how='left')

    if output_path is not None:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        sliced.to_csv(output_path, sep='\t')

    return sliced



def categorize_meth(
    meth_series,
    down_threshold=.2,
    up_threshold=.7,
    plot=False,
    quant_low=0.3,
    quant_high=0.7,
    kde_grid_points=1000,
    distance_frac=0.06,
    height_quantile=0.17,
    prominence_frac=0.01,
    peaks_dir="peaks",
):
    import re

    """
    Categorize methylation levels for a single probe using thresholds, KDE peak detection, or GMM clustering.

    :param meth_series: pandas.Series, beta-values of one probe across samples (non-null, numeric)
    :param down_threshold: float, threshold below which all values are labeled 'Not/Low_methyl'
    :param up_threshold: float, threshold above which all values are labeled 'Highly_methyl'
    :param plot: bool, if True, saves KDE+hist plots with peak lines into peaks_dir
    :param quant_low: float, fallback lower quantile for defining 'Not/Low_methyl' if peak/GMM fails
    :param quant_high: float, fallback upper quantile for defining 'Highly_methyl' if peak/GMM fails
    :param kde_grid_points: int, number of grid points for KDE calculation
    :param distance_frac: float, minimum distance between peaks as a fraction of kde_grid_points
    :param height_quantile: float, quantile threshold to filter low peaks based on KDE height
    :param prominence_frac: float, fraction of KDE std used to define minimum peak prominence
    :param peaks_dir: str or None, directory to save peak plots. If None, no files are written.
    :return: tuple of two elements:
             - pandas.Series with methylation categories: 'Not/Low_methyl', 'Medium_methyl', or 'Highly_methyl'
             - pandas.Series with up to three peak values: 'peak_1', 'peak_2', 'peak_3' (NaN if not found)
    """

    feature = meth_series.name
    meth_series = meth_series[meth_series.notna()].astype(float)

    # Work in a small DataFrame for convenience
    feature_values_i = pd.DataFrame(meth_series)
    feature_values_i['value'] = feature_values_i[feature]

    # Init peak result container
    peak_result = pd.Series(
        [np.nan, np.nan, np.nan],
        index=['peak_1', 'peak_2', 'peak_3'],
        name=feature
    )

    if feature_values_i[feature].max() < down_threshold:
        feature_values_i[feature] = 'Not/Low_methyl'
        return feature_values_i[feature], peak_result

    if feature_values_i[feature].min() >= up_threshold:
        feature_values_i[feature] = 'Highly_methyl'
        return feature_values_i[feature], peak_result

    # KDE
    kde = gaussian_kde(feature_values_i['value'])
    x = np.linspace(feature_values_i['value'].min(),
                    feature_values_i['value'].max(),
                    kde_grid_points)
    pdf = kde(x)

    height_threshold = np.quantile(pdf, height_quantile)
    peak_distance = int(kde_grid_points * distance_frac)
    peak_prominence = np.std(pdf) * prominence_frac

    maxima, properties = signal.find_peaks(
        pdf,
        height=height_threshold,
        distance=peak_distance,
        prominence=peak_prominence
    )

    peak_vals = x[maxima]
    sorted_peaks = np.sort(peak_vals)[::-1][:3]  # top 3, descending
    for i, val in enumerate(sorted_peaks):
        peak_result.iloc[i] = val

    num_peaks = len(sorted_peaks)

    # Decide thresholds based on number of peaks
    if num_peaks == 1:
        # Try 3-component GMM
        gmm = GaussianMixture(n_components=3, covariance_type='full')
        gmm.fit(meth_series.to_numpy().reshape(-1, 1))
        labels = gmm.predict(meth_series.to_numpy().reshape(-1, 1))
        feature_values_i['labels'] = labels

        if len(np.unique(labels)) == 3:
            means = feature_values_i.groupby('labels')[feature].mean().sort_values()
            down_cl_threshold = means.iloc[0]
            up_cl_threshold = means.iloc[2]
        else:
            # fallback to quantiles
            q_low = feature_values_i[feature].quantile(quant_low)
            q_high = feature_values_i[feature].quantile(quant_high)
            down_cl_threshold, up_cl_threshold = q_low, q_high

    elif num_peaks == 2:
        warnings.warn(f'Feature {feature} has bimodal distribution. Review manually.')
        down_cl_threshold, up_cl_threshold = sorted_peaks[1], sorted_peaks[0]

    elif num_peaks == 3:
        warnings.warn(f'Feature {feature} has multimodal distribution. Review manually.')
        down_cl_threshold = sorted_peaks[2]
        up_cl_threshold = np.median(sorted_peaks)

    else:
        
        q_low = feature_values_i[feature].quantile(quant_low)
        q_high = feature_values_i[feature].quantile(quant_high)
        down_cl_threshold, up_cl_threshold = q_low, q_high

    feature_values_i[feature] = np.where(
        feature_values_i['value'] < down_cl_threshold,
        'Not/Low_methyl',
        np.where(
            feature_values_i['value'] > up_cl_threshold,
            'Highly_methyl',
            'Medium_methyl'
        )
    )

    if plot and num_peaks > 0 and peaks_dir is not None:
        os.makedirs(peaks_dir, exist_ok=True)

        safe_name = re.sub(r'[^A-Za-z0-9_.-]', '_', str(feature))
        out_path = os.path.join(peaks_dir, f"{safe_name}.png")

        fig, ax = plt.subplots()
        sns.histplot(
            feature_values_i['value'],
            bins=40,
            stat='density',
            color='gray',
            kde=True,
            ax=ax
        )
        for val in sorted_peaks:
            if not pd.isna(val):
                ax.axvline(x=val, color='red', linestyle='--')
        ax.set_title(feature)
        ax.set_xlabel("Beta-value")
        ax.set_xlim(0, 1)
        fig.tight_layout()
        fig.savefig(out_path, dpi=300)
        plt.close(fig)

    return feature_values_i[feature], peak_result



def methylation(df,
                path=False,
                sep='\t',
                header=0,
                index_col=0,
                comment=None,
                save=True,
                output_dir="data",
                base_name="methylation",
                **kwargs):
    """
    Load and process methylation data, categorizing methylation levels and detecting peak values.

    Uses only:
      - 'chrom', 'start', 'end' as region metadata
      - all other columns (not in a fixed metadata list) as sample beta-values.

    :return: (full_categorized_df, peak_df)
             full_categorized_df: chrom/start/end + per-sample methylation categories
             peak_df: peak_1/2/3 per feature
    """

    if path:
        beta_values = pd.read_csv(df, sep=sep, header=header, index_col=index_col, comment=comment)
    else:
        beta_values = df.copy()

    non_sample_meta_cols = [
        'chrom', 'start', 'end',
        'id', 'score', 'strand', 'gene_name', 'gene_distance',
        'Stromal', 'Immune', 'epiCMIT.hyper', 'epiCMIT.hypo',
        'pomerantz', 'pom_region', 'regions2', 'loops',
        'coef', 'p_value', 'adjusted_p_value', 'q95l', 'q95h',
        'rsquared', 'gene', 'eCpGs', 'corr_plot_path'
    ]

    non_sample_meta_cols = [c for c in non_sample_meta_cols if c in beta_values.columns]

    sample_cols = [c for c in beta_values.columns if c not in non_sample_meta_cols]

    if len(sample_cols) == 0:
        raise ValueError("No sample columns detected: all columns are in the metadata list.")

    beta_only = beta_values[sample_cols].astype(float)

    region_cols = [c for c in ['chrom', 'start', 'end'] if c in beta_values.columns]
    metadata = beta_values[region_cols].copy()

    categorized_df = pd.DataFrame(index=beta_only.index, columns=beta_only.columns)
    peak_df = pd.DataFrame(index=beta_only.index, columns=['peak_1', 'peak_2', 'peak_3'])

    for row in tqdm(beta_only.index, desc="Categorizing methylation, detecting peaks"):
        categ, peaks = categorize_meth(beta_only.loc[row], **kwargs)
        categorized_df.loc[row] = categ
        peak_df.loc[row] = peaks

    full_categorized_df = pd.concat([metadata, categorized_df], axis=1)

    if save:
        os.makedirs(output_dir, exist_ok=True)
        labeled_path = os.path.join(output_dir, f"{base_name}_labeled.tsv")
        peaks_path = os.path.join(output_dir, f"{base_name}_peaks.tsv")

        full_categorized_df.to_csv(labeled_path, sep='\t')
        peak_df.to_csv(peaks_path, sep='\t')

    return full_categorized_df, peak_df


def filter_probes_by_methylation(methylation,
                                 threshold=0.8,
                                 save=True,
                                 output_path=None):
    """
    Filter probes based on the proportion of cases with 'Highly_methyl' or 'Medium_methyl' methylation levels.
    """
    filtered_probes = []

    for i in methylation.index:
        probe_data = methylation.loc[i, methylation.columns[3:]]
        counts = probe_data.value_counts()
        
        highly_methyl_proportion = counts.get('Highly_methyl', 0) / len(probe_data)
        medium_methyl_proportion = counts.get('Medium_methyl', 0) / len(probe_data)
        
        if highly_methyl_proportion < threshold and medium_methyl_proportion < threshold:
            filtered_probes.append(i)
    
    filtered_methylation = methylation.loc[filtered_probes]

    if save and output_path is not None:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        filtered_methylation.to_csv(output_path, sep='\t')

    return filtered_methylation
