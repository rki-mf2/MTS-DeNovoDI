def GET_RESULTS():
    """
        Get results binary file from different IDX types
                :param
                    None
                :return
                    file_name (str): path to the results file
        """
    counting_results = ""
    if (config['INDEX_TYPE'] == 'IBF'):
        if not os.path.exists("results/megaX/IBF/"):
            os.makedirs("results/megaX/IBF/")
        counting_results = "results/megaX/IBF/{sample}" +  "_results_IBF.log"
        
    elif (config['INDEX_TYPE'] == 'HIBF'):
        if not os.path.exists("results/megaX/HIBF/"):
            os.makedirs("results/megaX/HIBF/")
        counting_results = "results/megaX/HIBF/{sample}" +  "_results_HIBF.log"

    elif (config['INDEX_TYPE'] == 'HIBF_REF'):
        if not os.path.exists("results/megaX/HIBF_REF/"):
            os.makedirs("results/megaX/HIBF_REF/")
        counting_results = "results/megaX/HIBF_REF/{sample}_results_HIBF.log"

    elif (config['INDEX_TYPE'] == 'BS'):
        if not os.path.exists("results/megaX/BS/"):
            os.makedirs("results/megaX/BS/")
        counting_results = "results/megaX/BS/{sample}" +  "_results_BS.log"

    elif (config['INDEX_TYPE'] == 'FM_IDX'):
        if not os.path.exists("results/megaX/FM_IDX/"):
            os.makedirs("results/megaX/FM_IDX/")
        counting_results = "results/megaX/FM_IDX/{sample}" +  "_results_FM_index.log"
    else:
        print("Error: please choose one of the provided IDX types: ['IBF', 'HIBF', 'HIBF_REF', 'BS', 'FM_IDX']")
        exit()
    return counting_results

'''
Evaluation Rule
'''
rule evaluate:
    input:
       counting_results = GET_RESULTS(),
       reference_id_map = "results/megaX/stat/reference_id_map.log",
       reference_hist = "results/megaX/stat/reference_hist.log",
       lengths = "results/megaX/stat/lengths.log",
       query = "data/peptides/{sample}.fasta",
    output:
       evaluation_results = "results/megaX/classifiation/{sample}_counts.log",
       assignment_results = "results/megaX/classifiation/{sample}_assigned_peptides.log",
    log:
        stdout="results/logs/megaX/classifiation/{sample}_count.stdout", stderr="results/logs/megaX/classifiation/{sample}_count.stderr",
    conda:
        "../envs/megax.yaml"
    params:
        k_mer = config['GENERAL']['k_mer'],
        classification_threshold = config['CLASSIFICATION']['classification_threshold'],
    message: "megaX evaluation mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX evaluate -C {input.counting_results} -V {output.evaluation_results} -f {input.query}\
         -q {params.k_mer} -D {params.classification_threshold} -I {input.reference_id_map} > {log.stdout} 2> {log.stderr} &&\
         cp results/megaX/stat/assigned_peptides_idx.log {output.assignment_results}"

'''
Classification Rule
'''
rule classification:
    input:
       evaluation_results = "results/megaX/classifiation/{sample}_counts.log",
       assignment_results = "results/megaX/classifiation/{sample}_assigned_peptides.log",
    output:
        classifiation_results = "results/megaX/classifiation/{sample}_classification_report.txt",
       
    log:
        stdout="results/logs/megaX/classifiation/{sample}_classification_report.stdout", stderr="results/logs/megaX/classifiation/{sample}_classification_report.stderr",
    conda:
        "../envs/megax.yaml"
    message: "megaX classification mode!"
    threads: config['GENERAL']['threads']
    shell:
        "./workflow/tools/megaX classification -V {input.evaluation_results} -T {output.classifiation_results} > {log.stdout} 2> {log.stderr}"
