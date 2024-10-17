rule multik:
    input:
        rules.preprocessing.output,
        rules.clustering.output
    output:
        multik_adata = "../results/{sample_name}/da_multik/{sample_name}_{subsampling_iteration}/{sample_name}_iter_{subsampling_iteration}_res_{leiden_resolution}multik_adata.h5ad"
    params:
        results_dir = "../results/{sample_name}/da_multik/{sample_name}_{subsampling_iteration}",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"],
        subsampling_fraction = config["multik"]["subsampling_fraction"]
    shell:
        """
        mkdir -p {params.results_dir}
        .venv/bin/python3 scripts/da_multik.py \
        --sample_name {params.sample_name} \
        --results_dir {params.results_dir} \
        --plot_ext {params.plot_ext} \
        --subsampling_iteration {wildcards.subsampling_iteration} \
        --subsampling_fraction {params.subsampling_fraction} \
        --leiden_resolution {wildcards.leiden_resolution}
        """