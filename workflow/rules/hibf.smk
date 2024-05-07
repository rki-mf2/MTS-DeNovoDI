'''
Build vector of database (first IDX)
'''
rule build_vect_hibf:
    input:
       input_database = config['database']
    output:
       vect_idx = "results/megaX/HIBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect"
    log:
        stdout="results/logs/megaX/HIBF/vect.stdout", stderr="results/logs/megaX/HIBF/vect.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        output_dir = 'results/megaX/HIBF/',
        mutation_score = config['GENERAL']['score'],
        minimiser = config['GENERAL']['minimiser'],
        window_size = config['GENERAL']['window_size'],
    message: "MegaX build vect idx mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX build_vect -i {input.input_database} -m {params.matrix}\
         -t {threads} -s {params.mutation_score} -q {params.k_mer} -o {params.output_dir}\
          -w {params.window_size} -Z {params.minimiser} > {log.stdout} 2> {log.stderr}"

'''
Build HIBF of target vect (first IDX)
'''
rule build_hibf:
    input:
       vect_idx = "results/megaX/HIBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect"
    output:
       hibf = "results/megaX/HIBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + "_HIBF_files.log"
    log:
        stdout="results/logs/megaX/HIBF/HIBF.stdout", stderr="results/logs/megaX/HIBF/HIBF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        hibf_name = "results/megaX/HIBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".hibf",
        disable_estimation = config['HIBF']['disable_estimation'],
        disable_rearrangement = config['HIBF']['disable_rearrangement'],
        max_rearrangement = config['HIBF']['max_rearrangement'],
        maximum_false_positive_rate = config['HIBF']['maximum_false_positive_rate'],
        alpha = config['HIBF']['alpha'],
        sketch_bits = config['HIBF']['sketch_bits'],
        split_size = config['HIBF']['split_size'],
        hash_functions = config['HIBF']['hashes']
    message: "MegaX build HIBF idx mode!"
    threads: config['HIBF']['threads']
    shell:
        "./workflow/tools/megaX build_hibf -v {input.vect_idx} -F {params.hibf_name}\
        -E {params.disable_estimation} -R {params.disable_rearrangement} -g {params.max_rearrangement}\
        -r {params.maximum_false_positive_rate} -p {params.alpha} -S {params.sketch_bits}\
        -t {threads} -a {params.hash_functions} -M {params.split_size} > {log.stdout} 2> {log.stderr} &&\
         sed -i \"s|^|$(pwd)/|\" {output.hibf}"

'''
Count query sequences into HIBF(s)
'''
rule count_hibf:
    input:
       hibf_set = "results/megaX/HIBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + "_HIBF_files.log",
       query = "data/peptides/{sample}.fasta"
    output:
       counting_results = "results/megaX/HIBF/{sample}_results_HIBF.log"
    log:
        stdout="results/logs/megaX/HIBF/{sample}_results_HIBF.stdout", stderr="results/logs/megaX/HIBF/{sample}_results_HIBF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        threshold = config['HIBF']['threshold'],
    message: "MegaX count HIBF mode!"
    threads: config['HIBF']['threads']
    shell:
        "./workflow/tools/megaX count_hibf -F {input.hibf_set} -f {input.query}\
        -q {params.k_mer} -d {params.threshold} > {log.stdout} 2> {log.stderr}"

'''
Build HIBF of directly from input sequences database(no mutation)
'''
rule build_hibf_ref:
    input:
       input_database = config['database']
    output:
       hibf = "results/megaX/HIBF_REF/hibf_ref_" + str(config['GENERAL']['k_mer']) + "_HIBF_files.log"
    log:
        stdout="results/logs/megaX/HIBF_REF/HIBF_REF.stdout", stderr="results/logs/megaX/HIBF_REF/HIBF_REF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        hibf_name = "results/megaX/HIBF_REF/hibf_ref_" + str(config['GENERAL']['k_mer']) + ".hibf",
        disable_estimation = config['HIBF']['disable_estimation'],
        disable_rearrangement = config['HIBF']['disable_rearrangement'],
        max_rearrangement = config['HIBF']['max_rearrangement'],
        maximum_false_positive_rate = config['HIBF']['maximum_false_positive_rate'],
        alpha = config['HIBF']['alpha'],
        sketch_bits = config['HIBF']['sketch_bits'],
        split_size = config['HIBF']['split_size'],
        hash_functions = config['HIBF']['hashes']
    message: "MegaX build HIBF idx mode from input database!"
    threads: config['HIBF']['threads']
    shell:
        "./workflow/tools/megaX hibf_ref -i {input.input_database} -F {params.hibf_name}\
        -E {params.disable_estimation} -R {params.disable_rearrangement} -g {params.max_rearrangement}\
        -r {params.maximum_false_positive_rate} -p {params.alpha} -S {params.sketch_bits}\
        -t {threads} -a {params.hash_functions} -M {params.split_size} > {log.stdout} 2> {log.stderr} &&\
         sed -i \"s|^|$(pwd)/|\" {output.hibf}"

'''
Count query sequences into HIBF(s)
'''
rule count_hibf_ref:
    input:
       hibf_set = "results/megaX/HIBF_REF/hibf_ref_" + str(config['GENERAL']['k_mer']) + "_HIBF_files.log",
       query = "data/peptides/{sample}.fasta"
    output:
       counting_results = "results/megaX/HIBF_REF/{sample}_results_HIBF.log"
    log:
        stdout="results/logs/megaX/HIBF_REF/{sample}_results_HIBF.stdout", stderr="results/logs/megaX/HIBF_REF/{sample}_results_HIBF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        threshold = config['HIBF']['threshold'],
    message: "MegaX count HIBF_REF mode!"
    threads: config['HIBF']['threads']
    shell:
        "./workflow/tools/megaX count_hibf -F {input.hibf_set} -f {input.query}\
        -q {params.k_mer} -d {params.threshold} > {log.stdout} 2> {log.stderr}"