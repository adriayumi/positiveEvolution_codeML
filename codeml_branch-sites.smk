# conda environment: codeml

workdir: '/opt/adri/PR1/pipeline_codeml'

SAMPLES, = glob_wildcards('/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_nt-renamed.fasta')

ruleorder: mafft > macse > iqtree > M0 > bsA > bsA1

rule all:
    input:
        expand('/opt/adri/PR1/pipeline_codeml/null_model/M0.0.2_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/null_model/M0.0.7_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/null_model/M0.1.2_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.0.2_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.0.7_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.1.2_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.0.2_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.0.7_{sample}/out', sample=SAMPLES),
        expand('/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.1.2_{sample}/out', sample=SAMPLES)

rule mafft:
    input:
        '/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_renamed.fasta'
    output:
        '/opt/adri/PR1/pipeline_codeml/mafft/{sample}.mafft.out'
    shell:
        'mafft --localpair --maxiterate 1000 --inputorder {input} > {output}'

rule macse:
    input:
        AA = rules.mafft.output,
        NT = '/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_nt-renamed.fasta',
    output:
        '/opt/adri/PR1/pipeline_codeml/macse/{sample}.macse.out'
    shell:
        'java -jar /opt/adri/programas/macse_v2.01.jar -prog reportGapsAA2NT -align_AA {input.AA} -seq {input.NT} -out_NT {output}'

rule iqtree:
    input:
        rules.mafft.output,
    output:
        '/opt/adri/PR1/pipeline_codeml/mafft/{sample}.mafft.out.treefile'
    shell:
        'iqtree -s {input} -nt AUTO -bb 1000 -m TEST '
 
rule M0:
    input:
        TREEFILE = '/opt/adri/PR1/pipeline_codeml/mafft/{sample}.mafft.out.treefile',
        ALIGNMENT = '/opt/adri/PR1/pipeline_codeml/macse/{sample}.macse.out',
    output:
        '/opt/adri/PR1/pipeline_codeml/null_model/M0.0.2_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/null_model/M0.0.7_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/null_model/M0.1.2_{sample}/out',
    script:
        'run_model0.py'
        
rule bsA:
    input:
        TREEFILE = '/opt/adri/PR1/pipeline_codeml/mafft/{sample}.mafft.out.treefile',
        ALIGNMENT = '/opt/adri/PR1/pipeline_codeml/macse/{sample}.macse.out'
    output:
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.0.2_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.0.7_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA.1.2_{sample}/out',
    script:
        'run_modelbsA.py'
        
rule bsA1:
    input:
        TREEFILE = '/opt/adri/PR1/pipeline_codeml/mafft/{sample}.mafft.out.treefile',
        ALIGNMENT = '/opt/adri/PR1/pipeline_codeml/macse/{sample}.macse.out'
    output:
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.0.2_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.0.7_{sample}/out',
        '/opt/adri/PR1/pipeline_codeml/bioC_branchsite/workdir/bsA1.1.2_{sample}/out',
    script:
        'run_modelbsA1.py'
        
