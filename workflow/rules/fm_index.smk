'''
Build FM Index of database (first IDX)
'''
rule build_fm_idx:
    input:
       input_database = config['database']
    output:
       fm_idx = "results/megaX/FM_IDX/fm_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".idx"
    log:
        stdout="results/logs/megaX/FM_IDX/fm.stdout", stderr="results/logs/megaX/FM_IDX/fm.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        mutation_score = config['GENERAL']['score'],
        number_of_references = config['FM_INDEX']['number_of_references'],
    message: "MegaX build fm idx mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX build_fm -i {input.input_database} -m {params.matrix}\
         -t {threads} -s {params.mutation_score} -q {params.k_mer} -n {params.number_of_references}\
          -X {output.fm_idx}> {log.stdout} 2> {log.stderr}"

'''
Count query into fm index of the targeted database
'''
rule count_fm_idx:
    input:
       fm_idx = "results/megaX/FM_IDX/fm_" + str(config['GENERAL']['k_mer']) + "_" + str(config['GENERAL']['score']) + ".idx",
       query = "data/peptides/{sample}.fasta"
    output:
       counting_results = "results/megaX/FM_IDX/{sample}" +  "_results_FM_index.log"
    log:
        stdout="results/logs/megaX/FM_IDX/{sample}_results_FM_index.stdout", stderr="results/logs/megaX/FM_IDX/{sample}_results_FM_index.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        number_of_references = config['FM_INDEX']['number_of_references'],
    message: "MegaX fm idx search mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX count_fm -X {input.fm_idx} -f {input.query}\
         -t {threads} -q {params.k_mer} -n {params.number_of_references} > {log.stdout} 2> {log.stderr}"