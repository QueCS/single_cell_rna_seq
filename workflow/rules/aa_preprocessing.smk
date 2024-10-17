rule preprocessing:
    input:
        h5_matrix = "../resources/{sample_name}_filtered_feature_bc_matrix.h5"
    output:
        preprocessed_adata = "../results/{sample_name}/aa_preprocessing/{sample_name}_preprocessed_adata.h5ad"
    params:
        results_dir = "../results/{sample_name}/aa_preprocessing",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"]
    shell:
        """
        mkdir -p {params.results_dir}
        .venv/bin/python3 scripts/aa_preprocessing.py \
        --h5_matrix {input.h5_matrix} \
        --results_dir {params.results_dir} \
        --sample_name {params.sample_name} \
        --plot_ext {params.plot_ext}
        """