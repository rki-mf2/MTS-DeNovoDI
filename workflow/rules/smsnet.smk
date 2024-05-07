'''
Run SMSNet Sequencing Algorithm
'''
rule smsnet_seq:
    input:
       mgf_input =lambda wildcards:SAMPLES.at[wildcards.sample,'file_path'],
    output:
       output_dir =  directory("results/de_novo/smsnet/{sample}"),
    log:
        stdout="results/logs/de_novo/smsnet/{sample}.stdout", stderr="results/logs/de_novo/smsnet/{sample}.stderr",
    conda:
        "../envs/smsnet.yaml"
    params:
        gpu_id = config["DE_NOVO_SEQ"]["smsnet_gpu"],
        model = config["DE_NOVO_SEQ"]["smsnet_model"]
    threads: 50 #maximum; if cores < threads: threads = min(threads, cores)
    message: "Executing smsnet on the sample file {input}"
    shell:
        "{params.gpu_id} python ./workflow/de_novo_seq/SMSNet/run.py --model_dir {params.model} --inference_input_file {input.mgf_input} --inference_output_file {output.output_dir} > {log.stdout} 2> {log.stderr}"

'''
Filter Smsnet Peptides
'''
rule smsnet_fasta:
    input:
       input_file = rules.smsnet_seq.output,
    output:
       fasta_output = "results/de_novo/smsnet/{sample}.fasta",
    conda:
        "../envs/env.yaml"
    params:
        cutoff = config["DE_NOVO_SEQ"]["peptide_cutoff"],
    message: "Generating fasta file from smsnet peptides."
    script:
        "../scripts/processing/process_smsnet.py"