'''
Run both sequencing algorithms with combining output peptides
'''
rule combine_denovo:
    input:
       #smsnet_path = directory("results/de_novo/smsnet/{sample}"),
       smsnet_path = rules.smsnet_seq.output,
       casanovo_path  = "results/de_novo/casanovo/{sample}.mztab",
    output:
       output_file =  "results/de_novo/combined/{sample}.fasta",
    log:
        stdout="results/logs/de_novo/combined/{sample}.stdout", stderr="results/logs/de_novo/combined/{sample}.stderr",
    conda:
        "../envs/env.yaml"
    params:
        cutoff = config["DE_NOVO_SEQ"]["peptide_cutoff"],
        sample_name = "{sample}",
    threads: 50 #maximum; if cores < threads: threads = min(threads, cores)
    message: "Executing both sequencing tools on the sample file {input}"
    shell:
        "python ./workflow/scripts/processing/combine_denovo.py --sample_name {params.sample_name} --casanovo_path {input.casanovo_path} --smsnet_path {input.smsnet_path} --output_path {output.output_file} --cutoff {params.cutoff} > {log.stdout} 2> {log.stderr}"
