'''
Process Novor output file
'''
rule novor_sequences:
    input:
       novor_path  = lambda wildcards:SAMPLES.at[wildcards.sample,'file_path'],
    output:
       output_file =  "results/de_novo/novor/{sample}.fasta",
    log:
        stdout="results/logs/de_novo/novor/{sample}.stdout", stderr="results/logs/de_novo/novor/{sample}.stderr",
    conda:
        "../envs/env.yaml"
    params:
        cutoff = config["DE_NOVO_SEQ"]["peptide_cutoff"],
    threads: 50 #maximum; if cores < threads: threads = min(threads, cores)
    message: "Executing sequences extraction from Novor file: {input}"
    shell:
        "python ./workflow/scripts/processing/process_novor.py --novor_path {input.novor_path} --output_file {output.output_file} --cutoff {params.cutoff} > {log.stdout} 2> {log.stderr}"
