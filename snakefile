import os

# Extract dataset IDs from .h5 matrix files within data/
runs_samples = sorted(list([f.split("_filtered_feature_bc_matrix.h5")[0] for f in os.listdir("data/") if f.endswith("_filtered_feature_bc_matrix.h5")]))

# List script outputs
rule all:
    input:
        expand("output/{run_sample}/a_preprocessing/{run_sample}_pre_filter_vln.svg", run_sample=runs_samples),
        expand("output/{run_sample}/a_preprocessing/{run_sample}_pre_filter_scatter.svg", run_sample=runs_samples),
        expand("output/{run_sample}/a_preprocessing/{run_sample}_scrublet.svg", run_sample=runs_samples),
        expand("output/{run_sample}/a_preprocessing/{run_sample}_post_filter_vln.svg", run_sample=runs_samples),
        expand("output/{run_sample}/a_preprocessing/{run_sample}_post_filter_scatter.svg", run_sample=runs_samples),
        expand("output/{run_sample}/a_preprocessing/{run_sample}_preprocessed_adata.h5ad", run_sample=runs_samples),

rule preprocessing:
    input:
        matrix = "data/{run_sample}_filtered_feature_bc_matrix.h5"
    output:
        pre_filter_vln = "output/{run_sample}/a_preprocessing/{run_sample}_pre_filter_vln.svg",
        pre_filter_scatter = "output/{run_sample}/a_preprocessing/{run_sample}_pre_filter_scatter.svg",
        scrublet_plot = "output/{run_sample}/a_preprocessing/{run_sample}_scrublet.svg",
        post_filter_vln = "output/{run_sample}/a_preprocessing/{run_sample}_post_filter_vln.svg",
        post_filter_scatter = "output/{run_sample}/a_preprocessing/{run_sample}_post_filter_scatter.svg",
        preprocessed_adata = "output/{run_sample}/a_preprocessing/{run_sample}_preprocessed_adata.h5ad"
    params:
        file_ext = "svg"
    shell:
        """
        mkdir -p output/{wildcards.run_sample}/a_preprocessing
        python3 src/a_preprocessing.py --matrix {input.matrix} --run_sample {wildcards.run_sample} --file_ext {params.file_ext}
        """
