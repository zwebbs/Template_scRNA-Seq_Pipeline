# File Name: doublet_detection.py
# Created By: HK and ZW
# Created on: 2022-08-16
# Purpose: implements doublet and multiplet detection for 10x data

# Library + module imports
# ----------------------------------------------------------------------------
import argparse
import gzip
import numpy as np
from pathlib import Path
import scrublet as scr
from scipy.io import mmread


# Function Definitions
# ----------------------------------------------------------------------------
def load_genes(filename, column=1):
    genes_out = []
    with gzip.open(filename, 'rb') as fobj:
            for line in fobj:
                gene_id = line.decode('utf-8').strip().split('\t')[column]
                genes_out.append(gene_id)
    return genes_out

def load_barcodes(filename):
    barcodes_out = []
    with gzip.open(filename, 'rb') as fobj:
        for line in fobj:
            barcodes_out.append(line.decode('utf-8').strip())
    return barcodes_out

def write_doublets_out(filename, doublet_scores, doublet_predictions, barcodes):
    with open(filename, 'w') as fobj:
        data_out = zip(barcodes, doublet_scores, doublet_predictions)
        for b,s,p in data_out:
            fobj.write(f"{b},{s},{p}\n")

# MAIN execution block
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    # setup commandline parser
    descr = "Detects multiplets for removal in 10X scRNASeq data"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-d","--datadir", type=str,
        help="input directory for filtered counts matrix from cellranger")
    parser.add_argument("-n", "--sample-name", type=str,
        help="name of the sample for writing plots and output files.")
    parser.add_argument("-p", "--plotdir", type=str,
        help="output directory for any generated plots")
    parser.add_argument("-o", "--outdir", type=str,
        help="output directory for any generated files"),
    parser.add_argument("-r", "--prior-doublet-rate", type=float,
        help="prior for the model on how frequently doublets occur")

    # parser user arguments
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    plotdir = Path(args.plotdir).resolve()


    # Load counts matrix and gene list
    input_dir = Path(args.datadir).resolve()
    if not (input_dir.exists() and input_dir.is_dir()):
        raise ValueError("Error! data directory does not exist")

    counts_matrix = mmread(input_dir / 'matrix.mtx.gz').T.tocsc()
    genes = load_genes(input_dir / 'features.tsv.gz', column=0)
    barcodes = load_barcodes(input_dir / 'barcodes.tsv.gz')

    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))

    # Initialize Scrublet object
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.prior_doublet_rate)

    # Run the default pipeline, which includes: 1) Doublet simulation, 2) Normalization, gene filtering, rescaling, PCA,
    # 3) Doublet score calculation, 4) Doublet score threshold detection and doublet calling
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                        min_gene_variability_pctl=85, n_prin_comps=30)

    # Plot doublet score histograms for observed transcriptomes and simulated doublets
    ## if automatic threshold detection doesn't work well, you can adjust the threshold with the call_doublets() function.
    #scrub.call_doublets(threshold=0.25)
    hist_fig, hist_ax = scrub.plot_histogram()
    hist_fig.savefig(plotdir / (args.sample_name + "_scrublet_hist.png"))

    # Get 2-D embedding to visualize the results
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

    # Plot doublet predictions on 2-D embedding
    umap_fig, umap_ax = scrub.plot_embedding('UMAP', order_points=True)
    umap_fig.savefig(plotdir / (args.sample_name + "_scrublet_umap.png"))

    # write out the doublet scores and predictions for each cell
    write_doublets_out((outdir / (args.sample_name + "_scrublet_doublets.csv")),
        doublet_scores=doublet_scores, doublet_predictions=predicted_doublets,
        barcodes=barcodes)


    #args.sample_name
    print(doublet_scores)
    print(predicted_doublets)

############ load_genes function from Scrublet package
# def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
#     gene_list = []
#     gene_dict = {}
#
#     with open(filename) as f:
#         for iL in range(skip_rows):
#             f.readline()
#         for l in f:
#             gene = l.strip('\n').split(delimiter)[column]
#             if gene in gene_dict:
#                 gene_dict[gene] += 1
#                 gene_list.append(gene + '__' + str(gene_dict[gene]))
#                 if gene_dict[gene] == 2:
#                     i = gene_list.index(gene)
#                     gene_list[i] = gene + '__1'
#             else:
#                 gene_dict[gene] = 1
#                 gene_list.append(gene)
#     return gene_list
