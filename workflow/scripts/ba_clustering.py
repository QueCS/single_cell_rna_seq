# Import modules
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import skimage

# Set arguments parsing
arg_parser = argparse.ArgumentParser()
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
sample_name = args.sample_name
results_dir = args.results_dir
plot_ext = args.plot_ext

# Read the preprocessed adata .h5ad file
adata = sc.read_h5ad(
    f"../results/{sample_name}/aa_preprocessing/{sample_name}_preprocessed_adata.h5ad"
)

# Turn matplotlib interactive mode off
with plt.rc_context({"interactive": False, "savefig.format": plot_ext}):
    # Normalization
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Find the 2k HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(f"{results_dir}/{sample_name}_hvgs", bbox_inches="tight")
    plt.close()

    # Compute PCA
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig(f"{results_dir}/{sample_name}_pca_variance", bbox_inches="tight")
    plt.close()

    # Compute kNN Graph and UMAP
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50, metric="cosine")
    sc.tl.umap(adata, min_dist=0.3, maxiter=200)

    # Compute and Plot the Leiden clustering of the UMAP
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # Overlap QC metrics on UMAP
    metrics = ["leiden", "doublet_score", "n_counts", "n_genes"]
    for metric in metrics:
        sc.pl.umap(adata, color=metric, show=False)
        plt.savefig(f"{results_dir}/{sample_name}_umap_{metric}", bbox_inches="tight")
        plt.close()

# Save clustered adata as a .h5ad file
adata.write(f"{results_dir}/{sample_name}_clustered_adata.h5ad", compression="gzip")
