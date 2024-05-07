'''
Build vector of database (first IDX)
'''
rule build_vect_ibf:
    input:
       input_database = config['database']
    output:
       vect_idx = "results/megaX/IBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect"
    log:
        stdout="results/logs/megaX/IBF/vect.stdout", stderr="results/logs/megaX/IBF/vect.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        output_dir = 'results/megaX/IBF/',
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
Build IBF of target vect (first IDX)
'''
rule build_ibf:
    input:
       vect_idx = "results/megaX/IBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect"
    output:
       ibf = "results/megaX/IBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + "_" + str(config['IBF']['hashes']) + ".ibf"
    log:
        stdout="results/logs/megaX/IBF/IBF.stdout", stderr="results/logs/megaX/IBF/IBF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        output_dir = 'results/megaX/IBF/',
        hash_functions = config['IBF']['hashes']
    message: "MegaX build IBF idx mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX build_ibf -v {input.vect_idx} -o {params.output_dir}\
          -a {params.hash_functions} > {log.stdout} 2> {log.stderr}"

'''
Count IBF Rule
'''
rule count_ibf:
    input:
       ibf = "results/megaX/IBF/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + "_" + str(config['IBF']['hashes']) + ".ibf",
       query = "data/peptides/{sample}.fasta"
    output:
       counting_results = "results/megaX/IBF/{sample}_results_IBF.log"
    log:
        stdout="results/logs/megaX/IBF/{sample}_results_IBF.stdout", stderr="results/logs/megaX/IBF/{sample}_results_IBF.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
    message: "MegaX count IBF mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX count_ibf -F {input.ibf} -f {input.query} -q {params.k_mer} > {log.stdout} 2> {log.stderr}"

