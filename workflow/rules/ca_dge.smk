rule dge:
    input:
        rules.clustering.output
    output:
        dge_adata = "../results/{sample_name}/ca_dge/{sample_name}_dge_adata.h5ad"
    params:
        results_dir = "../results/{sample_name}/ca_dge",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"]
    shell:
        """
        mkdir -p {params.results_dir}
        .venv/bin/python3 scripts/ca_dge.py \
        --results_dir {params.results_dir} \
        --sample_name {params.sample_name} \
        --plot_ext {params.plot_ext}
        """