def GET_PEPTIDES_FILE():
    """
    Load peptides mode
            :param
                None
            :return
                peptides_fasta (str): path to peptides fasta file
    """

    if(config['sample_type'] == 'mgf'):
        print("[INFO] Start processing MS input data...")

        if(config['DE_NOVO_SEQ']['algorithm'] == 'casanovo'):
            casanovo_model = config['DE_NOVO_SEQ']['casanovo_model']
            try:
                with open(casanovo_model, 'r') as file:
                    print("[INFO] Casanovo Model found:", casanovo_model)
            except FileNotFoundError:
                print("[ERROR] Casanovo model is not accessible:", casanovo_model)
            peptides_file = "results/de_novo/casanovo/{sample}.fasta"
        elif(config['DE_NOVO_SEQ']['algorithm'] == 'smsnet'):
            peptides_file = "results/de_novo/smsnet/{sample}.fasta"
        elif(config['DE_NOVO_SEQ']['algorithm'] == 'both'):
            peptides_file = "results/de_novo/combined/{sample}.fasta"

    elif(config['sample_type'] == 'novor'):
        peptides_file = "results/de_novo/novor/{sample}.fasta"

    elif(config['sample_type'] == 'peptide_shaker'):
        peptides_file = "results/groundtruth/peptide_shaker/{sample}.fasta"
        
    elif(config['sample_type'] == 'fa'):
        query = lambda wildcards:SAMPLES.at[wildcards.sample,'file_path']
        peptides_file = query
    
    return peptides_file


def get_file_name(file):
    """
    Get file name without path
            :param
                file
            :return
                file_name (str): file name string extracted from complete path
    """
    import os

    return os.path.basename(file)

'''
Copy peptides file(s) to the corresponding pipeline folder
'''
rule cp_peptides:
    input:
       input_file = GET_PEPTIDES_FILE()
    output:
       fasta_output = "data/peptides/{sample}.fasta"
    log:
        stdout="results/logs/peptides/{sample}.stdout", stderr="results/logs/peptides/{sample}.stderr"
    conda:
        "../envs/megax.yaml"
    message: "Copy peptides to target folder."
    shell:
        "cp {input.input_file} {output.fasta_output} > {log.stdout} 2> {log.stderr}"

'''
Prepare input statistics for megaX usages
'''
rule stat:
    input:
       input_database = config['database']
    output:
       id_map = "results/megaX/stat/reference_id_map.log",
       hist = "results/megaX/stat/reference_hist.log",
       lengths = "results/megaX/stat/lengths.log",
    log:
        stdout="results/logs/megaX/stat/database.stdout", stderr="results/logs/megaX/stat/database.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        output_dir = 'results/megaX/stat/'
    message: "megaX stat mode!"
    shell:
        "./workflow/tools/megaX stat -q {params.k_mer} -i {input.input_database} -m {params.matrix}\
        -o {params.output_dir} > {log.stdout} 2> {log.stderr}"

'''
Print mutation statistics
'''
rule mutate_stat:
    input:
       input_database = config['database']
    output:
       mutate_output = "results/megaX/mutate_stat/mutate_stat.log",
    log:
        stderr="results/logs/megaX/mutate_stat/database.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        mutation_score = config['GENERAL']['score'],
        minimiser = config['GENERAL']['minimiser'],
        matrix = config['matrix'],
    message: "megaX mutate_stat mode!"
    shell:
        "./workflow/tools/megaX mutate_stat -q {params.k_mer} -s {params.mutation_score}\
         -i {input.input_database} -m {params.matrix} -Z {params.minimiser} > {output.mutate_output} 2> {log.stderr}"


'''
Write mutated database to output file
'''
rule write_db:
    input:
       input_database = config['database']
    output:
       mutated_db_file = "results/megaX/mutated_db/mutated_" + get_file_name(config['database']),
    log:
        stdout="results/logs/megaX/mutated_db/database.stdout", stderr="results/logs/megaX/mutated_db/database.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        matrix = config['matrix'],
        output_dir = 'results/megaX/mutated_db/',
        mutation_score = config['GENERAL']['score'],
    message: "megaX write_db mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX write_db -q {params.k_mer} -i {input.input_database} -m {params.matrix}\
        -t {threads} -o {params.output_dir} -s {params.mutation_score} > {log.stdout} 2> {log.stderr}"

'''
Write mutated database to output file
'''
rule simulate_peptides:
    input:
       input_database = config['database']
    output:
       simulated_peptides = "results/megaX/simulate_peptides/simulated_" + get_file_name(config['database']),
    log:
        stdout="results/logs/megaX/simulate_peptides/simulate_peptides.stdout", stderr="results/logs/megaX/simulate_peptides/simulate_peptides.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        errors = config['GENERAL']['simulation_errors'],
        target_peptides = config['GENERAL']['number_of_peptides'],
        output_dir = 'results/megaX/simulate_peptides/',
    message: "megaX simulate_peptides mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX simulate_peptides -i {input.input_database} -n {params.target_peptides} -o {params.output_dir}\
         -d {params.errors} > {log.stdout} 2> {log.stderr}"

'''
Write mutated database of simulated sequences to sequence length (instead of k_mers)
'''
rule mutate_seq_len:
    input:
       input_database = config['database']
    output:
       mutated_db_file = "results/megaX/mutated_db_len/mutated_" + get_file_name(config['database']),
    log:
        stdout="results/logs/megaX/mutated_db_len/simulate_peptides.stdout", stderr="results/logs/megaX/mutated_db_len/simulate_peptides.stderr"
    conda:
        "../envs/megax.yaml"
    params:
        number_of_muatetd_seq = config['GENERAL']['number_of_muatetd_seq'],
        matrix = config['matrix'],
        mutation_score = config['GENERAL']['score'],
        output_dir = 'results/megaX/mutated_db_len/',
    message: "megaX mutate_seq_len mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX mutate_seq_len -i {input.input_database} -N {params.number_of_muatetd_seq} -o {params.output_dir}\
         -m {params.matrix} -s {params.mutation_score} > {log.stdout} 2> {log.stderr}"