rule multik_compute:
    input:
        rules.preprocessing.output
    output:
        multik_compute_results = "../results/{sample_name}/da_multik/{sample_name}_iter_{subsampling_iteration}/{sample_name}_iter_{subsampling_iteration}_res_{leiden_resolution}_multik_results.csv"
    params:
        results_dir = "../results/{sample_name}/da_multik/{sample_name}_iter_{subsampling_iteration}",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"],
        subsampling_fraction = config["multik"]["subsampling_fraction"]
    shell:
        """
        mkdir -p {params.results_dir}
        .venv/bin/python3 scripts/da_multik_compute.py \
        --sample_name {params.sample_name} \
        --results_dir {params.results_dir} \
        --plot_ext {params.plot_ext} \
        --subsampling_iteration {wildcards.subsampling_iteration} \
        --subsampling_fraction {params.subsampling_fraction} \
        --leiden_resolution {wildcards.leiden_resolution}
        """

rule multik_collect:
    input:
        expand(rules.multik_compute.output, sample_name=sample_names, subsampling_iteration=subsampling_iterations, leiden_resolution=leiden_resolutions)
    output:
        multik_collect_results = "../results/{sample_name}/da_multik/{sample_name}_multik_results.csv",
    params:
        results_dir = "../results/{sample_name}/da_multik",
        sample_name = "{sample_name}",
        plot_ext = config["misc"]["plot_ext"]
    shell:
        """
        .venv/bin/python3 scripts/db_multik_collect.py --sample_name {params.sample_name} --results_dir {params.results_dir} --plot_ext {params.plot_ext}
        """