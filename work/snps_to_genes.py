import time
import pandas as pd

# Load SNP p-value association table from file
def load_SNP_pvals(snp_pval_file, delimiter='\t', header=False, cols='0,1,2,3'):
    # Check for valid 'cols' parameter
    try:
        cols_idx = [int(c) for c in cols.split(',')]
    except:
        raise ValueError('Invalid column index string')
    # Load gene_pos_file
    if header:
        SNP_summary = pd.read_csv(snp_pval_file, delimiter=delimiter)
    else:
        SNP_summary = pd.read_csv(snp_pval_file, delimiter=delimiter, header=-1)
    # Check gene positions table format
    if (SNP_summary.shape[1] < 4) | (max(cols_idx) >  SNP_summary.shape[1]-1):
        raise ValueError('Not enough columns in SNP Summary File')
    # Construct gene position table
    SNP_summary = SNP_summary[cols_idx]
    SNP_summary.columns = ['Marker', 'Chr', 'Pos', 'P-Value']
    return SNP_summary

# Load gene positions from file
def load_gene_pos(gene_pos_file, delimiter='\t', header=False, cols='0,1,2,3'):
    # Check for valid 'cols' parameter
    try:
        cols_idx = [int(c) for c in cols.split(',')]
    except:
        raise ValueError('Invalid column index string')
    # Load gene_pos_file
    if header:
        gene_positions = pd.read_csv(gene_pos_file, delimiter=delimiter)
    else:
        gene_positions = pd.read_csv(gene_pos_file, delimiter=delimiter, header=-1)
    # Check gene positions table format
    if (gene_positions.shape[1] < 4) | (max(cols_idx) >  gene_positions.shape[1]-1):
        raise ValueError('Not enough columns in Gene Positions File')
    # Construct gene position table
    gene_positions = gene_positions[cols_idx]
    gene_positions.columns = ['Gene', 'Chr', 'Start', 'End']
    return gene_positions.set_index('Gene')

# Assigning GWAS p-values to genes - Minimum P Method
#1. For each gene in the genome (or as defined by the Gene Positions file), we will collect all SNPs within a specified genomic distance from the gene body (transcription start site to transcription end site). The SNP must fall within the specified genomic distance (up or downstream of the gene body). This distance is given as kilobases, (e.g. if 'window' is set to 5, this will collect all SNPs within 5kb of the gene body.
#2. Each gene is then assigned the minimum of all the p-values across all SNPs falling within the specified window.

def min_p(SNP_summary, gene_positions, window):
    starttime = time.time()
    dist = window*1000
    genelist = list(gene_positions.index)
    min_p_list = []
    SNP_summary['Chr']=SNP_summary['Chr'].astype(str)
    for gene in genelist:
        gene_info = gene_positions.ix[gene]
        chrom = str(gene_info['Chr'])
        start = gene_info['Start']
        stop = gene_info['End']
        # Get all SNPs on same chromosome
        SNP_summary_filt1 = SNP_summary[SNP_summary['Chr']==chrom]
        # Get all SNPs after window start position
        SNP_summary_filt2 = SNP_summary_filt1[SNP_summary_filt1['Pos'] >= (start-dist)]
        # Get all SNPs before window end position
        SNP_summary_filt3 = SNP_summary_filt2[SNP_summary_filt2['Pos'] <= (stop+dist)]
        # Get min_p statistics for this gene
        if len(SNP_summary_filt3) >= 1:
            min_p_data = SNP_summary_filt3.ix[SNP_summary_filt3['P-Value'].argmin()]
            min_p_list.append([gene, chrom, start, stop, SNP_summary_filt3.shape[0], min_p_data['Marker'], int(min_p_data['Pos']), min_p_data['P-Value']])
        else:
            min_p_list.append([gene, chrom, start, stop, 0, None, None, None])
    min_p_table = pd.DataFrame(min_p_list, columns = ['Gene', 'Chr', 'Gene Start', 'Gene End', 'nSNPs', 'TopSNP', 'TopSNP Pos', 'TopSNP P-Value'])
    min_p_table['SNP Distance'] = abs(min_p_table['TopSNP Pos'].subtract(min_p_table['Gene Start']))
    min_p_table = min_p_table.dropna().sort_values(by=['TopSNP P-Value', 'Chr', 'Gene Start'])
    print "P-Values assigned to genes:", time.time()-starttime, 'seconds'
    return min_p_table

snp_summary_file = 'pgc.scz.full.2012-04.txt'
snp_summary = load_SNP_pvals(snp_summary_file, delimiter='\t', header=False, cols='0,1,2,7')

multiple_marks = snp_summary.index.value_counts()[snp_summary.index.value_counts() > 1].index

gene_pos_file = 'glist-hg18_proteinCoding.txt'
hg18_gene_pos = load_gene_pos(gene_pos_file, delimiter='\t', header=False)

min_p_table = min_p(snp_summary, hg18_gene_pos, 10)

min_p_table.to_csv('scz_gene_10k.txt',sep='\t')
