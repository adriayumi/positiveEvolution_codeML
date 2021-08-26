# conda environment: codeml

workdir: './'

SAMPLES, = glob_wildcards('/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_nt-renamed.fasta')

rule all:
    input:
        expand('mafft/{sample}.mafft.out', sample=SAMPLES),
        expand('macse/{sample}.macse.out', sample=SAMPLES),
        expand('mafft/{sample}.mafft.out.treefile', sample=SAMPLES),

rule mafft:
    input:
        '/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_renamed.fasta'
    output:
        'mafft/{sample}.mafft.out'
    threads: 10
    shell:
        'mafft --localpair --maxiterate 1000 --inputorder {input} > {output}'


rule macse:
    input:
        AA = 'mafft/{sample}.mafft.out',
        NT = '/opt/adri/PR1/PR1_orthogroups_fastas/fastas_headers_renomeados/{sample}_nt-renamed.fasta',
    output:
        'macse/{sample}.macse.out'
    threads: 10
    shell:
        'java -jar /opt/adri/programas/macse_v2.01.jar -prog reportGapsAA2NT -align_AA {input.AA} -seq {input.NT} -out_NT {output}'

rule iqtree:
    input:
        'mafft/{sample}.mafft.out', 
    output:
        'mafft/{sample}.mafft.out.treefile'
    threads: 10
    shell:
        'iqtree -s {input} -nt AUTO -bb 1000 -m TEST '
