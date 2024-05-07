'''
Histogram of raw database
'''
rule hist:
    input:
        input_hist = "results/megaX/stat/reference_hist.log"
    output:
        output_hist = "results/reports/reference_hist.pdf"
    conda:
        "../envs/env.yaml"
    log:
        log_file = "results/logs/reports/reference_hist.log",
    params: 
        depth = config['VISUALIZATION']['depth']
    script:
        "../scripts/visualization/hist.R"

'''
Visualizing peptides counts and lengths
'''
rule peptide_visualization:
    input:
        input_peptides = "data/peptides/{sample}.fasta"
    output:
        output_vis = "results/reports/peptides_vis_{sample}.pdf"
    conda:
        "../envs/env.yaml"
    log:
        log_file = "results/logs/reports/peptides_vis_{sample}.log",
    params: 
        k_mer = config['GENERAL']['k_mer'],
        sample_name = "{sample}"
    script:
        "../scripts/visualization/peptides_vis.R"

'''
Taxonomy visualization
'''
rule tax_visualization:
    input:
        input_taxonomy = "results/megaX/classifiation/{sample}_classification_report.txt"
    output:
        output_vis = "results/reports/{sample}_classification_report.pdf"
    conda:
        "../envs/env.yaml"
    log:
        log_file = "results/logs/reports/{sample}_classification_report.log",
    params: 
        k_mer = config['GENERAL']['k_mer'],
        depth = config['VISUALIZATION']['depth'],
        sample_name = "{sample}"
    script:
        "../scripts/visualization/taxonomy.R"

'''
Assignment visualization
'''
rule assignment_visualization:
    input:
        input_assignment = "results/megaX/classifiation/{sample}_counts.log",
    output:
        output_vis = "results/reports/{sample}_assignment_report.pdf"
    conda:
        "../envs/env.yaml"
    log:
        log_file = "results/logs/reports/{sample}_assignment_report.log",
    params: 
        k_mer = config['GENERAL']['k_mer'],
        depth = config['VISUALIZATION']['depth'],
        sample_name = "{sample}"
    script:
        "../scripts/visualization/assignment.R"