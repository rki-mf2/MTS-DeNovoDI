'''
Run Casanovo Sequencing Algorithm
'''
rule casanovo_seq:
    input:
       mgf_input =lambda wildcards:SAMPLES.at[wildcards.sample,'file_path'],
    output:
       mztab_output = "results/de_novo/casanovo/{sample}.mztab",
    log:
        stdout="results/logs/de_novo/casanovo/{sample}.stdout", stderr="results/logs/de_novo/casanovo/{sample}.stderr",
    conda:
        "../envs/casanovo.yaml"
    params:
        gpu_id = config["DE_NOVO_SEQ"]["casanovo_gpu"],
        model = config["DE_NOVO_SEQ"]["casanovo_model"],
        mode = config["DE_NOVO_SEQ"]["casanovo_mode"],
        tool_config = "configs/casanovo_tool_config.yaml"
    threads: 50 #maximum; if cores < threads: threads = min(threads, cores)
    message: "Executing casanovo on the sample file {input}"
    shell:
        "{params.gpu_id} casanovo {params.mode} -o {output.mztab_output} --model {params.model} --config {params.tool_config} {input.mgf_input}\
        > {log.stdout} 2> {log.stderr}"


'''
Filter Casanovo Peptides
'''
rule casanovo_fasta:
    input:
       mztab_input = "results/de_novo/casanovo/{sample}.mztab",
       #mztab_input = "/home/lutfia/scratch/casanovo/filtered_host.mztab"
    output:
       fasta_output = "results/de_novo/casanovo/{sample}.fasta",
    conda:
        "../envs/env.yaml"
    params:
        cutoff = config["DE_NOVO_SEQ"]["peptide_cutoff"],
    message: "Generating fasta file from casanovo peptides."
    script:
        "../scripts/processing/process_casanovo.py"