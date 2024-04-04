# Import modules
import os
import argparse
import matplotlib.pyplot as plt
import scanpy as sc
import skimage

# Make sure the working dir is src
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Set arguments parsing
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--matrix", type=str, required=True)
arg_parser.add_argument("--run_sample", type=str, required=True)
arg_parser.add_argument(
    "--file_ext",
    type=str,
    required=False,
    choices=["png", "jpg", "jpeg", "svg", "pdf", "eps", "tif", "tiff"],
    default="svg",
)
args = arg_parser.parse_args()

# Parse the arguments
matrix_path = f"../{args.matrix}"
run_sample = args.run_sample
file_ext = args.file_ext

output_path = f"../output/{run_sample}/a_preprocessing/"

# Turn matplotlib interactive mode off
with plt.rc_context({"interactive": False, "savefig.format": file_ext}):

    # Read the matrix file
    adata = sc.read_10x_h5(matrix_path)

    # Create the AnnData object
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.obs["n_genes"] = (adata.X > 0).sum(1)

    # QC plots, pre-filtering
    sc.pl.violin(adata, "n_counts", log=True, show=False)
    plt.savefig(f"{output_path}/{run_sample}_pre_filter_vln", bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, "n_counts", "n_genes", show=False)
    plt.savefig(f"{output_path}/{run_sample}_pre_filter_scatter", bbox_inches="tight")
    plt.close()

    # Filter the data
    # Remove cells with less than 200 detected genes
    # Remove genes present in less than 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Run scrublet to detected putative doublets and remove them
    sc.pp.scrublet(adata, expected_doublet_rate=0.03)

    sc.pl.scrublet_score_distribution(adata, show=False)
    plt.savefig(f"{output_path}/{run_sample}_scrublet", bbox_inches="tight")
    plt.close()

    scrublet_scores = adata.obs["doublet_score"]
    scrublet_predictions = adata.obs["predicted_doublet"]
    predicted_doublets = scrublet_scores[scrublet_predictions == True].index
    adata = adata[~adata.obs.index.isin(predicted_doublets)]

    # QC plots, post-filtering
    sc.pl.violin(adata, "n_counts", log=True, show=False)
    plt.savefig(f"{output_path}/{run_sample}_post_filter_vln", bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, "n_counts", "n_genes", show=False)
    plt.savefig(f"{output_path}/{run_sample}_post_filter_scatter", bbox_inches="tight")
    plt.close()

# Save preprocessed adata as a .h5ad file
adata.write(f"{output_path}/{run_sample}_preprocessed_adata.h5ad", compression="gzip")
