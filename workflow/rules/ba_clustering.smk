rule clustering:
    input:
        rules.preprocessing.output
    output:
        clustered_adata = "../results/{sample_name}/ba_clustering/{sample_name}_clustered_adata.h5ad"
    params:
        results_dir = "../results/{sample_name}/ba_clustering",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"]
    shell:
        """
        mkdir -p {params.results_dir}
        .venv/bin/python3 scripts/ba_clustering.py \
        --results_dir {params.results_dir} \
        --sample_name {params.sample_name} \
        --plot_ext {params.plot_ext}
        """