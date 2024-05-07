'''
Process PeptideShaker Report output file
'''
rule groundtruth_sequences:
    input:
       peptide_shaker_path  = lambda wildcards:SAMPLES.at[wildcards.sample,'file_path'],
    output:
       output_file =  "results/groundtruth/peptide_shaker/{sample}.fasta",
    log:
        stdout="results/logs/groundtruth/peptide_shaker/{sample}.stdout", stderr="results/logs/groundtruth/peptide_shaker/{sample}.stderr",
    conda:
        "../envs/env.yaml"
    threads: 50 #maximum; if cores < threads: threads = min(threads, cores)
    message: "Executing sequences extraction from groundtruth peptide_shaker file: {input}"
    shell:
        "python ./workflow/scripts/processing/process_peptideshaker.py --peptide_shaker_path {input.peptide_shaker_path} --output_file {output.output_file} > {log.stdout} 2> {log.stderr}"
