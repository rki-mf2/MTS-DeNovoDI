'''
Build vector of database (first IDX)
'''
rule build_vect_bs:
    input:
       input_database = config['database']
    output:
       vect_idx = "results/megaX/BS/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect"
    log:
        stdout="results/logs/megaX/BS/vect.stdout", stderr="results/logs/megaX/BS/vect.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        output_dir = 'results/megaX/BS/',
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
Count query into vector of database
'''
rule count_vect_bs:
    input:
       vect_idx = "results/megaX/BS/hashes_vect_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".vect",
       query = "data/peptides/{sample}.fasta"
    output:
       counting_results = "results/megaX/BS/{sample}" +  "_results_BS.log"
    log:
        stdout="results/logs/megaX/BS/{sample}_results_BS.stdout", stderr="results/logs/megaX/BS/{sample}_results_BS.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
    message: "MegaX binary search!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX binary_search -v {input.vect_idx} -f {input.query}\
         -t {threads} -q {params.k_mer} > {log.stdout} 2> {log.stderr}"