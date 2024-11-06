import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
import pickle
import pyBigWig
import statistics
import statistics as st
import seaborn as sns
sns.set_style("whitegrid")
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow, Arc
from collections import defaultdict
from sklearn.mixture import GaussianMixture
from statannot import add_stat_annotation
from statsmodels.stats.multitest import multipletests
from scipy.stats import gaussian_kde
from joblib import Parallel, delayed
import plotnine as p9
from scipy import signal
import pybedtools

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


gene_ref_paths = {'hg38': '~/38gencode.33.hugo.bed', 
                  'hg19': '~/gencode.v19.protein_coding_gene.minimal_chr.bed.gz'}
ctb_ref_paths = {'hg38': '~/38cytoBand.txt', 
                 'hg19': '~/19cytoBand.txt'}


def find_open_chromatin_probes(df_probes, df_atac):
    """
    Find probes in open chromatin regions.
    
    :param df_probes: DataFrame containing probe information with columns 'Chromosome', 'Start', 'End'.
    :param df_atac: DataFrame containing ATAC-seq open chromatin regions with columns 'chrom', 'start', 'end'.
    :return: DataFrame of probes that overlap with regions of open chromatin.
    """
    matching_rows = []

    for _, probe in tqdm(df_probes.iterrows(), total=df_probes.shape[0], desc="Finding matching probes"):
        
        mask_exact_match = (df_atac['chrom'] == probe['Chromosome']) & \
                           (df_atac['start'] == probe['Start']) & \
                           (df_atac['end'] == probe['End'])
        
        mask_inside_region = (df_atac['chrom'] == probe['Chromosome']) & \
                             (df_atac['start'] <= probe['Start']) & \
                             (df_atac['end'] >= probe['End'])
        
        mask = mask_exact_match | mask_inside_region
        
        if df_atac[mask].shape[0] > 0:
            matching_rows.append(probe)

    df_matching_probes = pd.DataFrame(matching_rows)
    return df_matching_probes


def probes2genes(meth_df, ref='hg19'):
    
    """
    Annotate methylation probes with closest genes.

    :param meth_df: pd.DataFrame, contains methylation probes with columns 'chrom', 'start', 'end', 'id', 'score', and 'strand'
    :param ref: str, default='hg19', genome reference for gene annotation (e.g., 'hg19' or 'hg38')
    :return: pd.DataFrame, probes annotated with closest gene and distance, indexed by unique gene-probe identifier
    """
                
    meth_df_copy = meth_df.sort_values(['chrom', 'start']).reset_index()[['chrom', 'start', 'end', 'id', 'score', 'strand']]  
    
    methyl = pybedtools.BedTool.from_dataframe(meth_df_copy)
    
    reference = read_dataset(gene_ref_paths[ref], index_col=False)
    
    reference['score']=0
    
    reference = reference[reference['#chr'].str.contains(r'^(chr[0-9XY]+)$')]
    
    reference = reference.rename(
        columns={'#chr':'chrom', 
                 'gene_start':'start', 
                 'gene_end':'end',
                 'gene_name':'gene'})[['chrom', 'start', 'end', 'gene', 'score', 'strand']]
    
    reference = reference.sort_values(['chrom', 'start'])
    
    reference = pybedtools.BedTool.from_dataframe(reference)
    
    annotated_metyl = methyl.closest(reference, t='all', s=True, d=True).to_dataframe(header=None)
    
    annotated_metyl = annotated_metyl.rename(
        columns={0:'chrom', 1:'start', 2:'end', 3:'id', 4:'score',
                 5:'strand', 6:'gene_chrom', 7:'gene_start', 8:'gene_end',
                 9:'gene_name', 10:'gene_score', 11:'gene_strand', 12:'gene_distance'})
    
    annotated_metyl=annotated_metyl[['chrom', 'start', 'end', 'id', 'score', 
                                     'strand', 'gene_name', 'gene_distance']]

    annotated_metyl['un_id'] = annotated_metyl['gene_name']+'_'+annotated_metyl['id'].astype(str)
    
    annotated_metyl = annotated_metyl.merge(meth_df.loc[:,~meth_df.columns.isin(annotated_metyl.columns)].reset_index(), 
                                      on='id', 
                                      how='right').set_index('un_id')


    return annotated_metyl

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


def pomerantz_state(me_filtered, reference_bed_path):
    
    """
    Annotate methylation data with Pomerantz chromatin states.

    :param me_filtered: pd.DataFrame, filtered methylation data with ['chrom', 'start', 'end', 'un_id']
    :param reference_bed_path: str, path to the Pomerantz reference BED file
    :return: pd.DataFrame with added 'pomerantz' and 'pom_region' columns for chromatin states
    """
    
    me_pom = me_filtered.reset_index()[['chrom', 'start', 'end', 'un_id']]
    me_forpomer = pybedtools.BedTool.from_dataframe(me_pom)
    
    # Load reference Pomerantz bed file
    po = pd.read_csv(reference_bed_path, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'state'])
    reference = pybedtools.BedTool.from_dataframe(po)

    intersection = me_forpomer.intersect(reference, loj=True).to_dataframe(header=None)
    intersection.columns = ['chrom', 'start', 'end', 'un_id', 'pom_chrom', 'pom_start', 'pom_end', 'pomerantz']
    
    intersection['pom_region'] = intersection['pom_chrom'] + ':' + intersection['pom_start'].astype(str) + '-' + intersection['pom_end'].astype(str)
    intersection = intersection.set_index('un_id')[['pomerantz', 'pom_region']]
    intersection['pomerantz'] = intersection['pomerantz'].replace('.', 'Unknown_state')
    
    # Replace long state names with shorter codes
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
        'Repressed_chromatin': 'Rep_chrom'
    }
    intersection['pomerantz'].replace(state_replacements, inplace=True)

    me_filtered = pd.concat([me_filtered, intersection], axis=1)

    return me_filtered

def process_region(reg, merge_probes, me_filtered2, correlation_threshold=0.6):
    
    """
    Process a single region and return merged or individual probe data.

    :param reg: str, region to process
    :param merge_probes: pd.DataFrame, contains probes to be merged
    :param me_filtered2: pd.DataFrame, methylation data
    :param correlation_threshold: float, default=0.6, threshold for merging probes within the same region
    :return: dict with keys 'indices', 'chromosomes', 'starts', 'ends', 'ids', 'scores', 'strands', 'gene_names', 'gene_distances', 'mes', 'pomerantz_list', 'pom_region_list', 'regions2_list'
    """
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
        'regions2_list': []
    }
    
    cols = ['chrom', 'start', 'end', 'id', 'score', 'strand', 
            'gene_name', 'gene_distance', 'pomerantz', 'pom_region', 'regions2']
    
    try:
        region_mask = merge_probes['regions2'] == reg
        probes_to_merge = merge_probes.loc[region_mask].index

        correlation_matrix = me_filtered2.loc[:, ~me_filtered2.columns.isin(cols)].loc[probes_to_merge].T.corr()
        
        # Check correlation threshold
        if check_correlation_threshold(correlation_matrix, correlation_threshold):
            index_base = pd.Series(probes_to_merge[0]).str.replace('_cg[0-9]+', '')[0] + '_' + merge_probes.loc[region_mask, 'pomerantz'].values[0]
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
            
            results['pomerantz_list'].append(merge_probes.loc[region_mask, 'pomerantz'].values[0])
            results['pom_region_list'].append(merge_probes.loc[region_mask, 'pom_region'].values[0])
            results['regions2_list'].append(merge_probes.loc[region_mask, 'regions2'].values[0])
            
        else:
            # cases where probes do not meet correlation threshold
            for probe in probes_to_merge:
                cpg_id = probe.split('_')[1]
                index_base = pd.Series(probe).str.replace('_cg[0-9]+', '')[0] + '_' + cpg_id + '_' + merge_probes.loc[probe, 'pomerantz']
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
                results['ids'].append(id_single)  # Use single id directly
                results['scores'].append(score)
                results['strands'].append(strand)
                results['gene_names'].append(gene_name)
                results['gene_distances'].append(gene_distance)

                me_i = me_filtered2.loc[probe].iloc[8:]
                results['mes'].append(me_i)
                
                results['pomerantz_list'].append(merge_probes.loc[probe, 'pomerantz'])
                results['pom_region_list'].append(merge_probes.loc[probe, 'pom_region'])
                results['regions2_list'].append(merge_probes.loc[probe, 'regions2'])
                
    except Exception as e:
        print(f"Error processing region {reg}: {e}")

    return results

def check_correlation_threshold(correlation_matrix, threshold=0.6):
    
    """
    Check if all correlations in a matrix exceed a given threshold.

    :param correlation_matrix: pd.DataFrame, correlation matrix for probes
    :param threshold: float, default=0.6, minimum correlation threshold
    :return: bool, True if all correlations exceed threshold, otherwise False
    """
    upper_triangle = correlation_matrix.where(
        np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
    ).stack()
    return (upper_triangle > threshold).all()

def process_probes(me_filtered2, correlation_threshold=0.6, num_cores=16):
    
    """
    Process and merge probes based on correlation thresholds using multiple cores.

    :param me_filtered2: pd.DataFrame, methylation data with necessary columns for probe processing
    :param correlation_threshold: float, default=0.6, threshold for merging probes within the same region
    :param num_cores: int, default=16, number of cores for parallel processing
    :return: pd.DataFrame with merged and unique probes
    """
    regions2 = me_filtered2.index.str.replace('_cg[0-9]+', '_') + me_filtered2['pom_region']
    me_filtered2['regions2'] = regions2
    
    regions = me_filtered2['regions2'].value_counts()[me_filtered2['regions2'].value_counts() > 1].index
    merge_probes = me_filtered2[me_filtered2['regions2'].isin(regions)]
    
    unique_probes = me_filtered2[~me_filtered2.index.isin(merge_probes.index)]
    unique_probes['pomerantz'] = unique_probes['pomerantz']
    unique_probes.index = unique_probes.index + '_' + unique_probes['pomerantz']
    
    global index_count
    index_count = defaultdict(int)
    
    results = Parallel(n_jobs=num_cores)(delayed(process_region)(reg, 
                                                                 merge_probes, 
                                                                 me_filtered2, 
                                                                 correlation_threshold
                                                                ) for reg in tqdm(merge_probes['regions2'].unique()))

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

    for col in me_filtered2.columns:
        if col not in combined_probes.columns:
            combined_probes[col] = np.nan

    combined_probes = combined_probes[me_filtered2.columns]

    final_combined = pd.concat([unique_probes, combined_probes])

    return final_combined


def resolve_duplicates(df):
    
    """
    Renames duplicate indices to ensure uniqueness.

    :param df: pd.DataFrame, DataFrame with potentially duplicate indices
    :return: pd.DataFrame with unique indices
    """
    df = df.copy()
    
    # Track seen indices and their counts
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


def annotate_probes_with_loops(probe_data, loop_data):
    
    """
    Annotate probe data with loop information from genomic loop data.

    :param probe_data: pd.DataFrame, contains probe data with 'chrom', 'start', 'end', 'un_id'
    :param loop_data: pd.DataFrame, loop information with 'chrom1', 'start1', 'end1', etc.
    :return: pd.DataFrame with probe data and new 'loops' column for associated loop IDs
    """
    probe_data.index.name='un_id'
    loops1 = loop_data[['chrom1', 'start1', 'end1', 'loop_id']]
    loops2 = loop_data[['chrom2', 'start2', 'end2', 'loop_id']]

    loops_bed1 = pybedtools.BedTool.from_dataframe(loops1)
    loops_bed2 = pybedtools.BedTool.from_dataframe(loops2)

    prob_gene_bed = pybedtools.BedTool.from_dataframe(probe_data.reset_index()[['chrom', 'start', 'end', 'un_id']])

    intersection1 = prob_gene_bed.intersect(loops_bed1, wa=True, wb=True).to_dataframe(
        names=['chrom', 'start', 'end', 'un_id', 'chrom1', 'start1', 'end1', 'loop_id1']
    )
    intersection2 = prob_gene_bed.intersect(loops_bed2, wa=True, wb=True).to_dataframe(
        names=['chrom', 'start', 'end', 'un_id', 'chrom2', 'start2', 'end2', 'loop_id2']
    )

    combined_intersection = pd.concat([
        intersection1[['un_id', 'loop_id1']], 
        intersection2[['un_id', 'loop_id2']]
    ])

    aggregated_loops = combined_intersection.groupby('un_id').agg(lambda x: list(set(x.dropna()))).reset_index()
    loop_dict = aggregated_loops.set_index('un_id').to_dict('index')

    probe_data['loops'] = probe_data.index.map(
        lambda x: loop_dict.get(x, {}).get('loop_id1', []) + loop_dict.get(x, {}).get('loop_id2', [])
    )

    probe_data['loops'] = probe_data['loops'].apply(lambda x: x[0] if len(x) == 1 else (x if x else np.nan))

    probe_data.index = probe_data.apply(
        lambda row: 'Ls_' + row.name if row['loops'] is not np.nan else row.name, axis=1
    )

    return probe_data


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
    for ind in tqdm(met.index):
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

def lin_reg(linreg_dict):
    
    """
    Perform linear regression for methylation and expression data, adjusting p-values.

    :param linreg_dict: dict, dictionary where keys are methylation feature indices and values are DataFrames with 'methylation' and 'expression' columns
    :return: pandas.DataFrame with regression results including coefficients, p-values, confidence intervals, and adjusted p-values
    """
    
    pattern_corrected = r'^(?:Ls_)?([A-Za-z0-9-]+)(?:_.*)?$'

    full_genes_meth = pd.DataFrame()

    for g in tqdm(list(linreg_dict.keys())):
        try:
            model = sm.OLS(endog=linreg_dict[g][['expression']], 
                           exog=median_scale(linreg_dict[g][['methylation']]),
                          data=linreg_dict[g])
            results = model.fit()

            r_squared = pd.DataFrame({'rsquared': [results.rsquared]}, index=['methylation'])

            gene_sum = pd.concat([pd.DataFrame(results.params[[0]]).rename(columns={0:'coef'}),
                                  pd.DataFrame(results.pvalues[[0]]).rename(columns={0:'p_value'}),
                                  results.conf_int().rename(columns={0:'q95l', 1:'q95h'}).iloc[[0]],
                                  r_squared], axis=1)

            full_genes_meth = pd.concat([full_genes_meth, gene_sum.rename(index={'methylation':g})], axis=0)
        except:
            print(g)

    full_genes_meth = full_genes_meth.sort_values('coef')
    
    full_genes_meth['eCpGs'] = full_genes_meth.index
    
    full_genes_meth['adjusted_p_value'] = multipletests(full_genes_meth['p_value'], method='fdr_bh')[1].astype(float)
    
    full_genes_meth['gene'] = full_genes_meth['eCpGs'].str.replace(pattern_corrected, r'\1', regex=True)
    
    return full_genes_meth

def forest_plot(conclusion_metrics,
                       title,
                       colour_interval=False,
                       figsize=(2.140, 2.194),
                       ylimits=(-8, 1),
                       feature_order=None,
                       metrics='Coefficient',
                       display_n=False,
                       text_family='Arial',
                       font_size=6):
    """
    param: conclusion_metrics, pd.DataFrame with the validation metrics
    param: title, str
    param: figsize, tuple
    param: ylimits, tuple, limits for the y axis, makes the plot more readable
    param: matrics, str, whether it's hazard or odds ratio
    param: text_family, str
    param: font_size, int, size of the font
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
        good_colour='#000096'
        bad_colour='#E61E1E'
        neutral_colour='#000000'
    else:
        good_colour='#000000'
        bad_colour='#000000'
        neutral_colour='#000000'
           
    p = (p9.ggplot(q)
         + p9.geom_pointrange(p9.aes(x='gene', y='coef', ymin='q95l', ymax='q95h'))
         + p9.coord_flip()
         + p9.scale_x_discrete(limits=feature_order)
         + p9.geom_hline(yintercept=0, color='#E61E1E', size=1, linetype="dashed")
         + p9.labs(x='gene', y=metrics)
         + p9.geom_text(p9.aes(x='gene', y='coef', label='p_value'),
                        size=font_size, nudge_x=0.4, nudge_y=0.1, family=text_family)
         + p9.ggtitle(title)
         + p9.ylim(ylimits)
         + p9.theme_bw()
         + p9.theme(
             text=p9.element_text(family=text_family, size=font_size),
             axis_title=p9.element_text(color='black', size=font_size),
             axis_text=p9.element_text(color='black', size=font_size),
             plot_title=p9.element_text(color='black', size=font_size)
         )
        )

    fig = p.draw()
    fig.set_size_inches(figsize[0], figsize[1])



def categorize_meth(meth_series, down_threshold=.2, up_threshold=.7, plot=False,
                    quant_low=0.3, quant_high=0.7):
    """
    Categorize methylation levels as 'Not/Low_methyl', 'Medium_methyl', or 'Highly_methyl'.

    :param meth_series: pandas.Series, one methylation feature (use only TUM/MET samples)
    :param down_threshold: float, threshold for low methylation (values below this are 'Not/Low_methyl')
    :param up_threshold: float, threshold for high methylation (values above this are 'Highly_methyl')
    :param plot: bool, if True, plots distribution with categories
    :param quant_low: float, lower quantile for medium methylation category
    :param quant_high: float, upper quantile for medium methylation category
    :return: pandas.Series with categorized methylation levels for each feature
    """
    feature = meth_series.name
    meth_series = meth_series[meth_series.notna()]
    meth_series = meth_series.astype(float)
    feature_values_i = pd.DataFrame(meth_series)
    feature_values_i['value'] = feature_values_i[feature]

    # annotate features with low methylation range [0:0.2]
    if feature_values_i[feature].max() < down_threshold:

        feature_values_i[feature] = 'Not/Low_methyl'

        if plot:
            sns.distplot(feature_values_i['value'],
                         hist=False, rug=True, color='#FFD3D3')
            plt.legend(labels=['Not/Low_methyl'])
            plt.xlabel(f'{feature} methylation')
            plt.xlim(0, 1)

        categorised = feature_values_i[feature]

    # annotate features with low methylation range [0.7:]
    elif feature_values_i[feature].min() >= up_threshold:

        feature_values_i[feature] = 'Highly_methyl'

        if plot:
            sns.distplot(feature_values_i['value'],
                         hist=False, rug=True, color='#D50000')
            plt.legend(labels=['Highly_methyl'])
            plt.xlabel(f'{feature} methylation')
            plt.xlim(0, 1)

        categorised = feature_values_i[feature]

    else:
        kde = gaussian_kde(feature_values_i[feature])
        x = np.linspace(feature_values_i[feature].min(),
                        feature_values_i[feature].max(), 1000)
        pdf = kde(x)
        
        height_threshold = np.quantile(pdf, 0.4)  #was 1.1
    
        prominence = np.std(pdf) * 0.03
        
        maxima, _ = signal.find_peaks(pdf, height=height_threshold, 
                                      prominence=prominence
                                     )

        # Identify the number and location of the peaks
        if len(maxima) == 1:
            peak_vals = pd.Series([x[maxima[0]]], index=['peak_1'])
        elif len(maxima) == 2:
            if pdf[maxima[0]] >= pdf[maxima[1]]:
                peak_vals = pd.Series([x[maxima[0]], x[maxima[1]]], index=['peak_1', 'peak_2'])
            else:
                peak_vals = pd.Series([x[maxima[1]], x[maxima[0]]], index=['peak_1', 'peak_2'])
        else:
            peak_vals = pd.Series([x[maxima[0]], x[maxima[1]], x[maxima[2]]], index=['peak_1', 'peak_2', 'peak_3'])

        # Save the peak values to a separate dataframe if the peak height is above a certain threshold
        peak_df = pd.DataFrame()
        for peak in peak_vals.index:
            if pdf[int(np.round((peak_vals[peak] - feature_values_i[feature].min()) / (
                    feature_values_i[feature].max() - feature_values_i[feature].min()) * 999))] > 1.1:
                peak_df[peak] = pd.Series([peak_vals[peak]])
        peak_df.index = [feature]

        # 1 peak
        if len(peak_df.columns) == 1:

            gmm = GaussianMixture(n_components=3, covariance_type='full')
            gmm.fit(meth_series.to_numpy().reshape(-1, 1))
            labels = gmm.predict(meth_series.to_numpy().reshape(-1, 1))

            feature_values_i['labels'] = labels

            if len(feature_values_i['labels'].unique()) == 3:

                cl0 = (feature_values_i[feature][feature_values_i.labels == 0]).mean()
                cl1 = (feature_values_i[feature][feature_values_i.labels == 1]).mean()
                cl2 = (feature_values_i[feature][feature_values_i.labels == 2]).mean()

                up_cl_threshold = np.median([cl0, cl1, cl2])
                down_cl_threshold = np.min([cl0, cl1, cl2])

                feature_values_i[feature] = np.where(feature_values_i[feature] < down_cl_threshold,
                                                     'Not/Low_methyl',
                                                     np.where(feature_values_i[feature] > up_cl_threshold,
                                                              'Highly_methyl', 'Medium_methyl'))
                categorised = feature_values_i[feature]

                if plot:
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Not/Low_methyl']['value'],
                                 hist=False, rug=True, color='#FFD3D3', kde_kws={'linestyle': '-'})
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Medium_methyl']['value'],
                                 hist=False, rug=True, color='#FF8787', kde_kws={'linestyle': '-'})
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Highly_methyl']['value'],
                                 hist=False, rug=True, color='#D50000', kde_kws={'linestyle': '-'})
                    plt.legend(labels=['Not/Low_methyl', 'Medium_methyl', 'Highly_methyl'])
                    plt.axvline(x=up_cl_threshold, color='black', linestyle=':')
                    plt.axvline(x=down_cl_threshold, color='black', linestyle=':')
                    plt.xlabel(f'{feature} methylation GMM')
                    plt.xlim(0, 1)
                else:
                    pass

            else:
                # for cases when GMM didn't work
                q_low = feature_values_i[feature].quantile(quant_low)
                q_high = feature_values_i[feature].quantile(quant_high)

                feature_values_i[feature] = np.where(feature_values_i[feature] < q_low,
                                                     'Not/Low_methyl',
                                                     np.where(feature_values_i[feature] > q_high,
                                                              'Highly_methyl', 'Medium_methyl'))
                categorised = feature_values_i[feature]

                if plot:
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Not/Low_methyl']['value'],
                                 hist=False, rug=True, color='#FFD3D3', kde_kws={'linestyle': '-'})
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Medium_methyl']['value'],
                                 hist=False, rug=True, color='#FF8787', kde_kws={'linestyle': '-'})
                    sns.distplot(feature_values_i[feature_values_i[feature] == 'Highly_methyl']['value'],
                                 hist=False, rug=True, color='#D50000', kde_kws={'linestyle': '-'})
                    plt.legend(labels=['Not/Low_methyl', 'Medium_methyl', 'Highly_methyl'])
                    plt.axvline(x=q_high, color='black', linestyle=':')
                    plt.axvline(x=q_low, color='black', linestyle=':')
                    plt.xlabel(f'{feature} methylation quantiles')
                    plt.xlim(0, 1)
                else:
                    pass


        # 2 peaks
        elif len(peak_df.columns) == 2:

            warnings.warn(f'This feature {feature} has bimodal distribution you may want to review it manually')

            up_cl_threshold = peak_df.loc[feature].max()
            down_cl_threshold = peak_df.loc[feature].min()

            feature_values_i[feature] = np.where(feature_values_i[feature] < down_cl_threshold,
                                                 'Not/Low_methyl',
                                                 np.where(feature_values_i[feature] >= up_cl_threshold,
                                                          'Highly_methyl', 'Medium_methyl'))
            if plot:
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Not/Low_methyl']['value'],
                             hist=False, rug=True, color='#FFD3D3', kde_kws={'linestyle': '-'})
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Medium_methyl']['value'],
                             hist=False, rug=True, color='#FF8787', kde_kws={'linestyle': '-'})
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Highly_methyl']['value'],
                             hist=False, rug=True, color='#D50000', kde_kws={'linestyle': '-'})
                plt.legend(labels=['Not/Low_methyl', 'Medium_methyl', 'Highly_methyl'])
                plt.axvline(x=up_cl_threshold, color='black', linestyle=':')
                plt.axvline(x=down_cl_threshold, color='black', linestyle=':')
                plt.xlabel(f'{feature} methylation KDE+finder')
                plt.xlim(0, 1)
            else:
                pass

            categorised = feature_values_i[feature]

        # 3 peaks
        elif len(peak_df.columns) == 3:

            warnings.warn(f'This feature {feature} has multimodal distribution you may want to review it manually')

            up_cl_threshold = peak_df.loc[feature].median()
            down_cl_threshold = peak_df.loc[feature].min()

            feature_values_i[feature] = np.where(feature_values_i[feature] < down_cl_threshold,
                                                 'Not/Low_methyl',
                                                 np.where(feature_values_i[feature] > up_cl_threshold,
                                                          'Highly_methyl', 'Medium_methyl'))
            if plot:
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Not/Low_methyl']['value'],
                             hist=False, rug=True, color='#FFD3D3', kde_kws={'linestyle': '-'})
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Medium_methyl']['value'],
                             hist=False, rug=True, color='#FF8787', kde_kws={'linestyle': '-'})
                sns.distplot(feature_values_i[feature_values_i[feature] == 'Highly_methyl']['value'],
                             hist=False, rug=True, color='#D50000', kde_kws={'linestyle': '-'})
                plt.legend(labels=['Not/Low_methyl', 'Medium_methyl', 'Highly_methyl'])
                plt.axvline(x=up_cl_threshold, color='black', linestyle=':')
                plt.axvline(x=down_cl_threshold, color='black', linestyle=':')
                plt.xlabel(f'{feature} methylation KDE+finder')
                plt.xlim(0, 1)
            else:
                pass

            categorised = feature_values_i[feature]

    return categorised

def methylation(df, path=False, sep='\t', header=0, index_col=0, comment=None, **kwargs):
    
    """
    Load and process methylation data, categorizing methylation levels.

    :param df: str or pandas.DataFrame, path to methylation data (if path=True) or DataFrame with methylation data
               (if path=False); rows are features, columns include 'chrom', 'start', 'end', and sample names
    :param path: bool, if True, reads data from the file path specified in df; if False, uses the DataFrame df
    :param sep: str, delimiter to use (passed to pandas.read_csv)
    :param header: int, row number to use as column names (passed to pandas.read_csv)
    :param index_col: int, column to use as row names (passed to pandas.read_csv)
    :param comment: str, character indicating lines to skip (passed to pandas.read_csv)
    :return: pandas.DataFrame with methylation data categorized into 'Not/Low_methyl', 'Medium methyl', 'Highly methyl'
    """
    if path:
        beta_values = pd.read_csv(df, sep=sep, index_col=index_col)
    else:
        beta_values = df

    regions = beta_values[['chrom', 'start', 'end']]

    beta_values = beta_values.loc[:, ~beta_values.columns.isin(['chrom', 'start', 'end'])].astype(float)

    processed_beta_val = beta_values.apply(lambda x: categorize_meth(x, **kwargs), axis=1)

    processed_beta_val = pd.concat([regions, processed_beta_val], axis=1)

    return processed_beta_val

def filter_probes_by_methylation(methylation, threshold=0.7):
    """
    Filter probes based on the proportion of cases with 'Highly_methyl' or 'Medium_methyl' methylation levels.
    
    This function filters out probes (rows) where 70% or more of the samples belong to either 
    'Highly_methyl' or 'Medium_methyl' categories. Probes with more diverse methylation patterns 
    are retained for further analysis.

    :param methylation: pandas.DataFrame, methylation data where rows represent probes with columns for 'chrom', 
                        'start', 'end', and each sample's methylation level.
    :param threshold: float, the proportion threshold for removing probes. If 70% or more of the samples for a 
                      probe are classified as either 'Highly_methyl' or 'Medium_methyl', the probe is removed.
                      Default is 0.7.
    :return: pandas.DataFrame with 'chrom', 'start', 'end', and all sample values for the probes that meet 
             the filtering criteria.
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
    
    return filtered_methylation

