# Import modules
import argparse
import pandas as pd
import glob
import matplotlib.pyplot as plt

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

# Loop to find all .csv files in +1 directories and merge them into a single {sample_name}.csv using the same columns
merged_data = pd.DataFrame()
for csv in glob.glob(f"{results_dir}/*/*.csv"):
    df = pd.read_csv(csv)
    if merged_data.empty:
        merged_data = df
    else:
        common_cols = list(set(merged_data.columns) & set(df.columns))
        if len(common_cols) > 0:
            merged_data = pd.merge(
                merged_data, df[common_cols], on=common_cols, how="outer"
            )

# Turn matplotlib interactive mode off
with plt.rc_context({"interactive": False}):
    # Plot a histogram showing n_clusters frequencies
    cluster_counts = merged_data["n_clusters"].value_counts().sort_index()
    plt.bar(
        cluster_counts.index, cluster_counts.values, color="grey", edgecolor="black"
    )
    plt.grid(axis="y")
    plt.title("Frequency distribution of cluster count")
    plt.xlabel("Cluster count")
    plt.ylabel("Frequency")
    plt.xticks(cluster_counts.index)
    plt.savefig(
        f"{results_dir}/{sample_name}_multik_cluster_count_frequency.{plot_ext}",
        bbox_inches="tight",
    )
    plt.close()

# Save merged_data as a .csv
merged_data.sort_values(by=["subsampling_iteration", "leiden_resolution"])
merged_data.to_csv(f"{results_dir}/{sample_name}_multik_results.csv", index=False)
