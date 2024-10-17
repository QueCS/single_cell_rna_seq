# Import modules
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib_venn as vlt
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
    f"../results/{sample_name}/ba_clustering/{sample_name}_clustered_adata.h5ad"
)

# Turn matplotlib interactive mode off
with plt.rc_context({"interactive": False, "savefig.format": plot_ext}):
    # Determine how many clusters are present
    n_clusters = len(set(adata.obs["leiden"]))
    clusters = [[str(i)] for i in range(n_clusters)]

    # Compute DGE and export associated plots
    test_methods = ["t-test", "t-test_overestim_var", "wilcoxon"]
    for test_method in test_methods:
        sc.tl.rank_genes_groups(
            adata, "leiden", method=test_method, key_added=test_method
        )
        for cluster in clusters:
            sc.pl.rank_genes_groups(adata, key=test_method, groups=cluster, show=False)
            plt.savefig(
                f"{results_dir}/{sample_name}_{test_method}_cluster{cluster[0]}_top20_hvg_scores",
                bbox_inches="tight",
            )
            plt.close()

    for cluster in clusters:
        wilcoxon_names = adata.uns["wilcoxon"]["names"][cluster[0]][:100]
        tt_names = adata.uns["t-test"]["names"][cluster[0]][:100]
        ttov_names = adata.uns["t-test_overestim_var"]["names"][cluster[0]][:100]
        vlt.venn3(
            [set(wilcoxon_names), set(tt_names), set(ttov_names)],
            (
                "rank-sum test",
                "t-test",
                "t-test\n(overestimated variance)",
            ),
            layout_algorithm=vlt.layout.venn3.DefaultLayoutAlgorithm(
                fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1)
            ),
        )
        plt.savefig(
            f"{results_dir}/{sample_name}_cluster{cluster[0]}_top100_hvg_intersection"
        )
        plt.close()

# Save analysed adata as a .h5ad file
adata.write(f"{results_dir}/{sample_name}_dge_adata.h5ad", compression="gzip")
