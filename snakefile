import pandas as pd

# Extract sample data from data/data.csv
data_df = pd.read_csv("data/data.csv")
sample_names = data_df["sample_name"]

# List script outputs
rule all:
    input:
        expand("output/{sample_name}/a_preprocessing/{sample_name}_pre_filter_vln.svg", sample_name=sample_names),
        expand("output/{sample_name}/a_preprocessing/{sample_name}_pre_filter_scatter.svg", sample_name=sample_names),
        expand("output/{sample_name}/a_preprocessing/{sample_name}_scrublet.svg", sample_name=sample_names),
        expand("output/{sample_name}/a_preprocessing/{sample_name}_post_filter_vln.svg", sample_name=sample_names),
        expand("output/{sample_name}/a_preprocessing/{sample_name}_post_filter_scatter.svg", sample_name=sample_names),
        expand("output/{sample_name}/a_preprocessing/{sample_name}_preprocessed_adata.h5ad", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_hvgs.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_pca_variance.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_umap_louvain.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_umap_doublet_score.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_umap_n_counts.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_umap_n_genes.svg", sample_name=sample_names),
        expand("output/{sample_name}/b_clustering/{sample_name}_clustered_adata.h5ad", sample_name=sample_names)

rule preprocessing:
    output:
        pre_filter_vln = "output/{sample_name}/a_preprocessing/{sample_name}_pre_filter_vln.svg",
        pre_filter_scatter = "output/{sample_name}/a_preprocessing/{sample_name}_pre_filter_scatter.svg",
        scrublet_plot = "output/{sample_name}/a_preprocessing/{sample_name}_scrublet.svg",
        post_filter_vln = "output/{sample_name}/a_preprocessing/{sample_name}_post_filter_vln.svg",
        post_filter_scatter = "output/{sample_name}/a_preprocessing/{sample_name}_post_filter_scatter.svg",
        preprocessed_adata = "output/{sample_name}/a_preprocessing/{sample_name}_preprocessed_adata.h5ad"
    params:
        file_ext = "svg"
    shell:
        """
        mkdir -p output/{wildcards.sample_name}/a_preprocessing
        python3 src/a_preprocessing.py --sample_name {wildcards.sample_name} --file_ext {params.file_ext}
        """

rule clustering:
    input:
        preprocessed_adata = rules.preprocessing.output.preprocessed_adata
    output:
        hvgs_plot = "output/{sample_name}/b_clustering/{sample_name}_hvgs.svg",
        pca_variance_plot = "output/{sample_name}/b_clustering/{sample_name}_pca_variance.svg",
        umap_louvain = "output/{sample_name}/b_clustering/{sample_name}_umap_louvain.svg",
        umap_doublet_score = "output/{sample_name}/b_clustering/{sample_name}_umap_doublet_score.svg",
        umap_n_counts = "output/{sample_name}/b_clustering/{sample_name}_umap_n_counts.svg",
        umap_n_genes = "output/{sample_name}/b_clustering/{sample_name}_umap_n_genes.svg",
        clustered_adata = "output/{sample_name}/b_clustering/{sample_name}_clustered_adata.h5ad"
    params:
        file_ext = "svg"
    shell:
        """
        mkdir -p output/{wildcards.sample_name}/b_clustering
        python3 src/b_clustering.py --sample_name {wildcards.sample_name} --file_ext {params.file_ext}
        """