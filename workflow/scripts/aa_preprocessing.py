# Import modules
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import skimage

# Set arguments parsing
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--h5_matrix", type=str, required=True)
arg_parser.add_argument("--sample_name", type=str, required=True)
arg_parser.add_argument("--results_dir", type=str, required=True)
arg_parser.add_argument(
    "--plot_ext",
    type=str,
    required=False,
    choices=["png", "jpg", "jpeg", "svg", "pdf", "eps", "tif", "tiff"],
    default="svg",
)

# Parse the arguments
args = arg_parser.parse_args()
h5_matrix = args.h5_matrix
sample_name = args.sample_name
results_dir = args.results_dir
plot_ext = args.plot_ext

# Turn matplotlib interactive mode off
with plt.rc_context({"interactive": False, "savefig.format": plot_ext}):
    # Read the matrix file
    adata = sc.read_10x_h5(f"{h5_matrix}")

    # Create the AnnData object
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.obs["n_genes"] = (adata.X > 0).sum(1)

    # QC plots, pre-filtering
    sc.pl.violin(adata, "n_counts", log=True, show=False)
    plt.savefig(f"{results_dir}/{sample_name}_pre_filter_vln", bbox_inches="tight")
    plt.close()
    sc.pl.scatter(adata, "n_counts", "n_genes", show=False)
    plt.savefig(f"{results_dir}/{sample_name}_pre_filter_scatter", bbox_inches="tight")
    plt.close()

    # Filter the data
    # Remove cells with less than 200 detected genes
    # Remove genes present in less than 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Run scrublet to estimate doublets rate and filter doublets out
    adata_sim = sc.pp.scrublet_simulate_doublets(adata)
    sc.pp.scrublet(adata, adata_sim)

    sc.pl.scrublet_score_distribution(adata, show=False)
    plt.savefig(f"{results_dir}/{sample_name}_scrublet", bbox_inches="tight")
    plt.close()

    scrublet_scores = adata.obs["doublet_score"]
    scrublet_predictions = adata.obs["predicted_doublet"]
    predicted_doublets = scrublet_scores[scrublet_predictions == True].index
    adata = adata[~adata.obs.index.isin(predicted_doublets)]

    # QC plots, post-filtering
    sc.pl.violin(adata, "n_counts", log=True, show=False)
    plt.savefig(f"{results_dir}/{sample_name}_post_filter_vln", bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, "n_counts", "n_genes", show=False)
    plt.savefig(f"{results_dir}/{sample_name}_post_filter_scatter", bbox_inches="tight")
    plt.close()

# Save preprocessed adata as a .h5ad file
adata.write(f"{results_dir}/{sample_name}_preprocessed_adata.h5ad", compression="gzip")
