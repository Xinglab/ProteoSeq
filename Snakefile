import snakemake.utils
import pandas as pd
import os

snakemake.utils.min_version('6.5.3')

configfile: 'snakemake_config.yaml'

onsuccess:
    print('workflow success')

onerror:
    print('workflow error')

DEFAULT_MEM_MB=16 * 1024  # 16 GB
DEFAULT_TIME_HOURS=12

if not os.path.exists(config['msgf_path']):
    exit('MS-GF+ is not set!')

if (not os.path.exists(config['RNA_gtf'])) and (not os.path.exists(config['protein_fasta'])):
    exit('Neither RNA_gtf or protein_fasta should be input!')

if (not config['reverse']) and (not config['shuffle']):
    exit('Please select at least one method (reverse or shuffle) to construct the decoy database!')

DATASETS = sorted(config['ms_datasets'])
MS_PATH_DICT = {}
DATASET_CONFIG = {}
DATASET_ALL_MS_RUNS = {}
if DATASETS:
    for dataset in DATASETS:
        ms_sample_inputs = config['ms_datasets'][dataset]
        mzml_dir = ms_sample_inputs.get('mzml_dir')
        DATASET_CONFIG[dataset] = ms_sample_inputs.get('msgf_config')
        for f in os.listdir(mzml_dir):
            if not f.endswith('.mzML'):
                continue
            ms_run = f[:-5]
            MS_PATH_DICT[(dataset, ms_run)] = os.path.join(mzml_dir, f)
            if ms_run not in DATASET_ALL_MS_RUNS.setdefault(dataset, []):
                DATASET_ALL_MS_RUNS[dataset].append(ms_run)


def all_input(wildcards):
    inputs = dict()
    inputs['protein'] = os.path.join(config['reference_dir'], 'protein_db.fa')
    inputs['index'] = os.path.join(config['reference_dir'], 'protein_db.index.txt')
    inputs['buildSA'] = os.path.join(config['reference_dir'], 'buildSA.done')
    if len(DATASET_ALL_MS_RUNS) == 0:
        return inputs
    inputs['novel_psms'] = os.path.join(config['results_dir'], 'all_sig_psms_count.novel.txt')
    inputs['uniprot_psms'] = os.path.join(config['results_dir'], 'all_sig_psms_count.uniprot.txt')
    inputs['gencode_psms'] = os.path.join(config['results_dir'], 'all_sig_psms_count.gencode.txt')
    if (not config['protein_fasta']) or (config['protein_fasta'] == ''):
        if config['RNA_gtf'] != '':
            inputs['novel_psms_mapping'] = os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.rna_sample.gtf')
            if (config['enable_visualization']):
                inputs['vis_dir']=os.path.join(config['visualization_dir'], 'plot_novel_peps.done')
    if config['predict_gene_by_PSM_count']:
        inputs['gene_abun']=os.path.join(config['results_dir'], 'all_gene_level_count.txt')
    if config['report_intensity']:
        inputs['intensity'] = os.path.join(config['results_dir'], 'all_sig_psms_intensity.log2.txt')
    if (config['predict_gene_by_intensity']) and (config['report_intensity']):
        inputs['gene_abun_intensity']=os.path.join(config['results_dir'], 'all_gene_level_intensity.log2.txt')
    return inputs


localrules: all
rule all:
    input:
        unpack(all_input),


def get_rna_gtf(wildcards):
    if not config['predict_ORF']:
        return config['RNA_gtf']
    return os.path.join(config['reference_dir'], 'protein_predicted.gtf')


def get_genome_fasta(wildcards):
    if (config['genome_fasta_path']) and (config['genome_fasta_path'] != ''):
        return config['genome_fasta_path']
    return os.path.join(config['reference_dir'], 'GRCh38.primary_assembly.genome.fa')


def get_genome_gtf(wildcards):
    if (config['GENCODE_gtf_path']) and (config['GENCODE_gtf_path'] != ''):
        return config['GENCODE_gtf_path']
    return os.path.join(config['reference_dir'], 'gencode.v39.annotation.gtf')


def get_gencode_protein(wildcards):
    if (config['GENCODE_protein_path']) and (config['GENCODE_protein_path'] != ''):
        return config['GENCODE_protein_path']
    return os.path.join(config['reference_dir'], 'gencode.v39.pc_transcripts.fa')


def get_uniprot_protein(wildcards):
    if (config['UniProt_protein_path']) and (config['UniProt_protein_path'] != ''):
        return config['UniProt_protein_path']
    return os.path.join(config['reference_dir'], 'uniprot_sprot_wVarsplic_human.PC.fa')


rule download_genome_fasta:
    output:
        genome_fasta=os.path.join(config['reference_dir'], 'GRCh38.primary_assembly.genome.fa'),
    log:
        out=os.path.join(config['log_dir'], 'download_genome_fasta.log'),
    params:
        genome_fasta_url=config['genome_fasta_url'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl \'{params.genome_fasta_url}\' -o {output.genome_fasta}.gz >> {log.out} 2>&1'
        ' && '
        'gunzip {output.genome_fasta}.gz >> {log.out} 2>&1'


rule download_gencode_gtf:
    output:
        gencode_gtf=os.path.join(config['reference_dir'], 'gencode.v39.annotation.gtf'),
    log:
        out=os.path.join(config['log_dir'], 'download_gencode_gtf.log'),
    params:
        gencode_gtf_url=config['GENCODE_gtf_url'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl \'{params.gencode_gtf_url}\' -o {output.gencode_gtf}.gz > {log.out} 2>&1'
        ' && '
        'gunzip {output.gencode_gtf}.gz >> {log.out} 2>&1'


rule download_gencode_protein:
    output:
        gencode_protein=os.path.join(config['reference_dir'], 'gencode.v39.pc_transcripts.fa'),
    log:
        out=os.path.join(config['log_dir'], 'download_gecode_protein.log'),
    params:
        gencode_protein_url=config['GENCODE_protein_url'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl \'{params.gencode_protein_url}\' -o {output.gencode_protein}.gz > {log.out} 2>&1'
        ' && '
        'gunzip {output.gencode_protein}.gz >> {log.out} 2>&1'


rule download_uniprot_reference:
    output:
        uniprot=os.path.join(config['reference_dir'], 'uniprot_sprot.fasta'),
        uniprotvar=os.path.join(config['reference_dir'], 'uniprot_sprot_varsplic.fasta'),
    log:
        out=os.path.join(config['log_dir'], 'download_uniprot_reference.log'),
    params:
        dir=config['reference_dir'],
        uniprot_url=config['UniProt_protein_url'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl \'{params.uniprot_url}\' | tar zxvf - -C {params.dir} \'uniprot_sprot.fasta.gz\' \'uniprot_sprot_varsplic.fasta.gz\' > {log.out} 2>&1'
        ' && '
        'gunzip {output.uniprot}.gz >> {log.out} 2>&1'
        ' && '
        'gunzip {output.uniprotvar}.gz >> {log.out} 2>&1'


rule select_species_in_uniprot:
    input:
        uniprot_fasta=os.path.join(config['reference_dir'], 'uniprot_sprot.fasta'),
        uniprotvar_fasta=os.path.join(config['reference_dir'], 'uniprot_sprot_varsplic.fasta'),
    output:
        combined_fasta=os.path.join(config['reference_dir'], 'uniprot_sprot_wVarsplic.fa'),
        fasta=os.path.join(config['reference_dir'], 'uniprot_sprot_wVarsplic_human.fa'),
    params:
        species='9606',
        scripts=config['scripts_dir'],
        conda_wrapper=config['conda_wrapper'],
    log:
        out=os.path.join(config['log_dir'], 'download_uniprot_reference.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cat {input.uniprot_fasta} > {output.combined_fasta}\" >> {log.out}'
        ' && '
        'echo \"Command: cat {input.uniprotvar_fasta} >> {output.combined_fasta}\" >> {log.out}'
        ' && '
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/select_species_in_uniprot.py'
        ' -i {output.combined_fasta} -o {output.fasta} --species {params.species}\" >> {log.out}'
        ' && '
        'cat {input.uniprot_fasta} > {output.combined_fasta}'
        ' && '
        'cat {input.uniprotvar_fasta} >> {output.combined_fasta}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/select_species_in_uniprot.py'
        ' -i {output.combined_fasta} -o {output.fasta} --species {params.species}'
        ' >> {log.out} 2>&1'


rule filter_ptc_in_uniprot:
    input:
        fasta=os.path.join(config['reference_dir'], 'uniprot_sprot_wVarsplic_human.fa'),
    output:
        fasta=os.path.join(config['reference_dir'], 'uniprot_sprot_wVarsplic_human.PC.fa'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'download_uniprot_reference.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/filter_ptc.py'
        ' -i {input.fasta} -o {output.fasta}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/filter_ptc.py'
        ' -i {input.fasta} -o {output.fasta}'
        ' >> {log.out} 2>&1'


rule annotate_gtf:
    input:
        sample_gtf=config['RNA_gtf'],
        genome_gtf=get_genome_gtf,
    output:
        gtf=temp(os.path.join(config['reference_dir'], '1_annotated_wCDS.gtf')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/annotate_gtf.py -v 1'
        ' -i {input.sample_gtf} -a {input.genome_gtf} -o {output.gtf} --selected_features CDS'
        ' > {log.out}\"'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_gtf.py -v 1'
        ' -i {input.sample_gtf} -a {input.genome_gtf} -o {output.gtf} --selected_features CDS'
        ' >> {log.out} 2>&1'


rule select_pc_transcripts:
    input:
        gtf=os.path.join(config['reference_dir'], '1_annotated_wCDS.gtf'),
    output:
        pc_gtf=temp(os.path.join(config['reference_dir'], '1_basic_pc.gtf')),
        other_gtf=temp(os.path.join(config['reference_dir'], '1_non_basic_pc.gtf')),
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    params:
        tx_type='protein_coding',
        tx_tag='basic',
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/select_transcripts_in_gtf.py'
        ' -i {input.gtf} -o {output.pc_gtf} --output_discarded {output.other_gtf}'
        ' --tx_type {params.tx_type} --tx_tag {params.tx_tag}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/select_transcripts_in_gtf.py'
        ' -i {input.gtf} -o {output.pc_gtf} --output_discarded {output.other_gtf}'
        ' --tx_type {params.tx_type} --tx_tag {params.tx_tag}'
        ' >> {log.out} 2>&1'


rule select_nmd_transcripts:
    input:
        gtf=os.path.join(config['reference_dir'], '1_non_basic_pc.gtf'),
    output:
        nmd_gtf=temp(os.path.join(config['reference_dir'], '1_nmd.gtf')),
        other_gtf=temp(os.path.join(config['reference_dir'], '1_non_basic_pc_nmd.gtf')),
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    params:
        tx_type='nonsense_mediated_decay',
        tx_wCDS='True',
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/select_transcripts_in_gtf.py'
        ' -i {input.gtf} -o {output.nmd_gtf}'
        ' --output_discarded {output.other_gtf} --tx_type {params.tx_type} --tx_wCDS {params.tx_wCDS}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/select_transcripts_in_gtf.py'
        ' -i {input.gtf} -o {output.nmd_gtf}'
        ' --output_discarded {output.other_gtf} --tx_type {params.tx_type} --tx_wCDS {params.tx_wCDS}'
        ' >> {log.out} 2>&1'


rule translate_basic_pc:
    input:
        genome_fasta=get_genome_fasta,
        gtf=os.path.join(config['reference_dir'], '1_basic_pc.gtf'),
    output:
        protein_fasta=os.path.join(config['reference_dir'], 'protein.basic_pc.fa'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/translate_by_gtf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.protein_fasta}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/translate_by_gtf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.protein_fasta}'
        ' >> {log.out} 2>&1'


rule predict_orf:
    input:
        genome_fasta=get_genome_fasta,
        gtf=os.path.join(config['reference_dir'], '1_non_basic_pc_nmd.gtf'),
    output:
        gtf=temp(os.path.join(config['reference_dir'], '1_predicted_orf.gtf')),
        protein_fasta=os.path.join(config['reference_dir'], 'protein.predicted_orf.fa'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/predict_orf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.gtf} --output_protein {output.protein_fasta}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/predict_orf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.gtf} --output_protein {output.protein_fasta}'
        ' >> {log.out} 2>&1'


rule predict_nmd:
    input:
        genome_fasta=get_genome_fasta,
        gtf=os.path.join(config['reference_dir'], '1_non_basic_pc_nmd.gtf'),
    output:
        gtf=temp(os.path.join(config['reference_dir'], '1_predicted_nmd.gtf')),
        protein_fasta=os.path.join(config['reference_dir'], 'protein.predicted_nmd.fa'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_nmd.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/predict_orf.py -v 1'
        ' --keep_nmd True'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.gtf} --output_protein {output.protein_fasta}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/predict_orf.py -v 1'
        ' --keep_nmd True'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.gtf} --output_protein {output.protein_fasta}'
        ' >> {log.out} 2>&1'


rule add_orf_tag_to_pc_gtf:
    input:
        pc_gtf=os.path.join(config['reference_dir'], '1_basic_pc.gtf'),
    output:
        pc_gtf=temp(os.path.join(config['reference_dir'], '1_basic_pc.anno.gtf')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.pc_gtf} -o {output.pc_gtf} --tag ORF:GENCODE\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.pc_gtf} -o {output.pc_gtf} --tag ORF:GENCODE'
        ' >> {log.out} 2>&1'


rule add_orf_tag_to_novel_gtf:
    input:
        novel_gtf=os.path.join(config['reference_dir'], '1_predicted_orf.gtf'),
    output:
        novel_gtf=temp(os.path.join(config['reference_dir'], '1_predicted_orf.anno.gtf')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:novel\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:novel'
        ' >> {log.out} 2>&1'


rule combine_pc_gtf:
    input:
        annotated_gtf=os.path.join(config['reference_dir'], '1_basic_pc.anno.gtf'),
        novel_gtf=os.path.join(config['reference_dir'], '1_predicted_orf.anno.gtf'),
    output:
        gtf=temp(os.path.join(config['reference_dir'], 'protein_predicted.gtf.temp')),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    shell:
        'echo \"Command: cat {input.annotated_gtf} > {output.gtf}\" >> {log.out}'
        ' && '
        'echo \"Command: cat {input.novel_gtf} >> {output.gtf}\" >> {log.out}'
        ' && '
        'cat {input.annotated_gtf} > {output.gtf}'
        ' && '
        'cat {input.novel_gtf} >> {output.gtf}'


rule add_orf_tag_to_annotated_nmd_gtf:
    input:
        novel_gtf=os.path.join(config['reference_dir'], '1_nmd.gtf'),
    output:
        novel_gtf=temp(os.path.join(config['reference_dir'], '1_nmd.anno.gtf')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:NMD\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:NMD'
        ' >> {log.out} 2>&1'


rule add_orf_tag_to_predicted_nmd_gtf:
    input:
        novel_gtf=os.path.join(config['reference_dir'], '1_predicted_nmd.gtf'),
    output:
        novel_gtf=temp(os.path.join(config['reference_dir'], '1_predicted_nmd.anno.gtf')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:NMD\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/add_tag_to_gtf.py'
        ' -i {input.novel_gtf} -o {output.novel_gtf} --tag ORF:NMD'
        ' >> {log.out} 2>&1'


rule combine_nmd_gtf:
    input:
        annotated_gtf=os.path.join(config['reference_dir'], '1_nmd.anno.gtf'),
        novel_gtf=os.path.join(config['reference_dir'], '1_predicted_nmd.anno.gtf'),
    output:
        gtf=temp(os.path.join(config['reference_dir'], 'nmd_predicted.gtf.temp')),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    shell:
        'echo \"Command: cat {input.annotated_gtf} > {output.gtf}\" >> {log.out}'
        ' && '
        'echo \"Command: cat {input.novel_gtf} >> {output.gtf}\" >> {log.out}'
        ' && '
        'cat {input.annotated_gtf} > {output.gtf}'
        ' && '
        'cat {input.novel_gtf} >> {output.gtf}'


rule combine_pc_nmd_gtf:
    input:
        pc_gtf=os.path.join(config['reference_dir'], 'protein_predicted.gtf.temp'),
        nmd_gtf=os.path.join(config['reference_dir'], 'nmd_predicted.gtf.temp'),
    output:
        gtf=temp(os.path.join(config['reference_dir'], 'protein_predicted.wnmd.gtf.temp')),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    shell:
        'echo \"Command: cat {input.pc_gtf} > {output.gtf}\" >> {log.out}'
        ' && '
        'echo \"Command: cat {input.nmd_gtf} >> {output.gtf}\" >> {log.out}'
        ' && '
        'cat {input.pc_gtf} > {output.gtf}'
        ' && '
        'cat {input.nmd_gtf} >> {output.gtf}'


def select_final_predicted_gtf(wildcards):
    if (not config['discard_NMD']): #include NMD in protein db
        return os.path.join(config['reference_dir'], 'protein_predicted.wnmd.gtf.temp')
    return os.path.join(config['reference_dir'], 'protein_predicted.gtf.temp')


rule get_final_predicted_gtf:
    input:
        temp_db=select_final_predicted_gtf,
    output:
        final_db=os.path.join(config['reference_dir'], 'protein_predicted.gtf'),
    log:
        out=os.path.join(config['log_dir'], 'predict_orf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cat {input.temp_db} > {output.final_db}\" >> {log.out}'
        ' && '
        'cat {input.temp_db} > {output.final_db}'


rule translate_pc:
    input:
        genome_fasta=get_genome_fasta,
        gtf=get_rna_gtf,
    output:
        protein_fasta=os.path.join(config['reference_dir'], 'protein_from_gtf.fa'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'translate_protein.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/translate_by_gtf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.protein_fasta}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/translate_by_gtf.py -v 1'
        ' -g {input.genome_fasta} -i {input.gtf} -o {output.protein_fasta}'
        ' >> {log.out} 2>&1'


def get_sample_protein(wildcards):
    if config['protein_fasta'] != '':
        return config['protein_fasta']
    return os.path.join(config['reference_dir'], 'protein_from_gtf.fa')


rule collapse_protein_db:
    input:
        pc_fasta=get_sample_protein,
    output:
        fasta=temp(os.path.join(config['reference_dir'], 'protein_db.woCanonical.temp')),
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp0')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'collapse_protein_seqs.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/collapse_protein_seqs.py'
        ' -i {input.pc_fasta} -o {output.fasta}'
        ' -ot {output.index_table} --add_tags ORF_type:Novel\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/collapse_protein_seqs.py'
        ' -i {input.pc_fasta} -o {output.fasta}'
        ' -ot {output.index_table} --add_tags ORF_type:Novel'
        ' >> {log.out} 2>&1'


def get_canonical_protein_fasta(wildcards):
    canonical_ref_list = []
    if config['merge_UniProt_protein']:
        canonical_ref_list.append(get_uniprot_protein(wildcards))
    if config['merge_GENCODE_protein']:
        canonical_ref_list.append(get_gencode_protein(wildcards))
    if (config['customized_canonical_protein_path']) and (config['customized_canonical_protein_path'] != ''):
        canonical_ref_list.append(config['customized_canonical_protein_path'])
    return canonical_ref_list


def get_canonical_protein_tag(wildcards):
    tag_list = []
    if config['merge_UniProt_protein']:
        tag_list.append('UniProt')
    if config['merge_GENCODE_protein']:
        tag_list.append('GENCODE')
    if (config['customized_canonical_protein_path']) and (config['customized_canonical_protein_path'] != ''):
        tag_list.append('Canonical')
    return ','.join(tag_list)


rule collapse_protein_db_with_canonical:
    input:
        pc_fasta=get_sample_protein,
        canonical_fasta=get_canonical_protein_fasta,
    output:
        fasta=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.temp')),
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp0')),
    params:
        input_tag_list=get_canonical_protein_tag,
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'collapse_protein_seqs.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/collapse_protein_seqs.py'
        ' -i {input.canonical_fasta} {input.pc_fasta}'
        ' -o {output.fasta} -ot {output.index_table} --add_tags ORF_type:{params.input_tag_list},Novel\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/collapse_protein_seqs.py'
        ' -i {input.canonical_fasta} {input.pc_fasta}'
        ' -o {output.fasta} -ot {output.index_table} --add_tags ORF_type:{params.input_tag_list},Novel'
        ' >> {log.out} 2>&1'


rule annotate_index_table_wocanonical_by_genome_gtf:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp0'),
        gtf=get_genome_gtf,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp1')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}'
        ' >> {log.out} 2>&1'


rule annotate_index_table_wocanonical_by_rna_gtf:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp1'),
        gtf=get_rna_gtf,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp2')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}'
        ' >> {log.out} 2>&1'


rule annotate_index_table_wcanonical_by_genome_gtf:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp0'),
        gtf=get_genome_gtf,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp1')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}'
        ' >> {log.out} 2>&1'
        

rule annotate_index_table_wcanonical_by_rna_gtf:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp1'),
        gtf=get_rna_gtf,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp2')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_gtf.py'
        ' -i {input.index_table} -o {output.index_table} -g {input.gtf}'
        ' >> {log.out} 2>&1'


rule annotate_index_table_wcanonical_wo_rna_gtf_by_uniprot:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp1'),
        uniprot=get_uniprot_protein,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp3')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_uniprot.py'
        ' -i {input.index_table} -o {output.index_table} -f {input.uniprot}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_uniprot.py'
        ' -i {input.index_table} -o {output.index_table} -f {input.uniprot}'
        ' >> {log.out} 2>&1'


rule annotate_index_table_wcanonical_w_rna_gtf_by_uniprot:
    input:
        index_table=os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp2'),
        uniprot=get_uniprot_protein,
    output:
        index_table=temp(os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp4')),
    log:
        out=os.path.join(config['log_dir'], 'annotate_protein_db.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/annotate_table_by_uniprot.py'
        ' -i {input.index_table} -o {output.index_table} -f {input.uniprot}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/annotate_table_by_uniprot.py'
        ' -i {input.index_table} -o {output.index_table} -f {input.uniprot}'
        ' >> {log.out} 2>&1'


def select_final_protein_db(wildcards):
    if (not config['merge_UniProt_protein']) and (not config['merge_GENCODE_protein']):
        if (not config['customized_canonical_protein_path']) or (config['customized_canonical_protein_path'] == ''):
            return os.path.join(config['reference_dir'], 'protein_db.woCanonical.temp')
    return os.path.join(config['reference_dir'], 'protein_db.wCanonical.temp')


rule get_final_db:
    input:
        temp_db=select_final_protein_db,
    output:
        final_db=os.path.join(config['reference_dir'], 'protein_db.fa'),
    log:
        out=os.path.join(config['log_dir'], 'get_final_protein_db.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cat {input.temp_db} > {output.final_db}\" >> {log.out}'
        ' && '
        'cat {input.temp_db} > {output.final_db}'


def select_final_index(wildcards):
    if (not config['merge_UniProt_protein']) and (not config['merge_GENCODE_protein']):
        if (config['protein_fasta']) and (config['protein_fasta'] != ''): #Ignore config['RNA_gtf'], only annotate with gencode_gtf
            return os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp1')
        else: #Annotate with gencode_gtf, then rna_gtf
            return os.path.join(config['reference_dir'], 'protein_db.woCanonical.index.temp2')
    elif (not config['merge_UniProt_protein']):
        if (config['protein_fasta']) and (config['protein_fasta'] != ''): #Ignore config['RNA_gtf'], only annotate with gencode_gtf
            return os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp1')
        else: #Annotate with gencode_gtf, then rna_gtf
            return os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp2')
    else:
        if (config['protein_fasta']) and (config['protein_fasta'] != ''): #Ignore config['RNA_gtf'], only annotate with gencode_gtf, then uniprot
            return os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp3')
        else: #Annotate with gencode_gtf, then rna_gtf, then uniprot
            return os.path.join(config['reference_dir'], 'protein_db.wCanonical.index.temp4')


rule get_final_index:
    input:
        temp_index=select_final_index,
    output:
        final_index=os.path.join(config['reference_dir'], 'protein_db.index.txt'),
    log:
        out=os.path.join(config['log_dir'], 'get_final_index_db.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cp {input.temp_index} {output.final_index}\" > {log.out}'
        ' && '
        'cp {input.temp_index} {output.final_index}'


def get_decoy_method(wildcards):
    decoy_method = ''
    if config['reverse']:
        decoy_method += ' --reverse'
    if config['shuffle']:
        decoy_method += ' --shuffle'
    return decoy_method


rule get_decoy:
    input:
        protein_fasta=os.path.join(config['reference_dir'], 'protein_db.fa'),
    output:
        decoy_fasta=os.path.join(config['reference_dir'], 'protein_db.mimic.fa'),
    log: out=os.path.join(config['log_dir'], 'buildSA.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
        decoy_method = get_decoy_method,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: python {params.scripts}/reverse_protein_seqs.py'
        ' -i {input.protein_fasta} -o {output.decoy_fasta} --prefix Random {params.decoy_method}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/reverse_protein_seqs.py'
        ' -i {input.protein_fasta} -o {output.decoy_fasta} --prefix Random {params.decoy_method}'
        ' >> {log.out} 2>&1'


rule build_sa:
    input:
        protein_fasta=os.path.join(config['reference_dir'], 'protein_db.fa'),
        decoy_fasta=os.path.join(config['reference_dir'], 'protein_db.mimic.fa'),
    output:
        touch(os.path.join(config['reference_dir'], 'buildSA.done')),
    log: out=os.path.join(config['log_dir'], 'buildSA.log'),
    params:
        cmd=config['msgf_path'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: java -Xmx3500M -cp {params.cmd} edu.ucsd.msjava.msdbsearch.BuildSA -d {input.protein_fasta} -tda 0\"'
        ' >> {log.out}'
        ' && '
        'java -Xmx3500M -cp {params.cmd} edu.ucsd.msjava.msdbsearch.BuildSA -d {input.protein_fasta} -tda 0'
        ' >> {log.out} 2>&1'
        ' && '
        'echo \"Command: java -Xmx3500M -cp {params.cmd} edu.ucsd.msjava.msdbsearch.BuildSA -d {input.decoy_fasta} -tda 0\"'
        ' >> {log.out}'
        'java -Xmx3500M -cp {params.cmd} edu.ucsd.msjava.msdbsearch.BuildSA -d {input.decoy_fasta} -tda 0'
        ' >> {log.out} 2>&1'


rule create_work_dir:
    output:
        touch(os.path.join(config['work_dir'], 'mk_work_dir.done')),


rule run_msgf_target:
    input:
        mzml=lambda wildcards: MS_PATH_DICT[(wildcards.dataset, wildcards.ms_run)],
        protein_fasta=os.path.join(config['reference_dir'], 'protein_db.fa'),
        decoy_fasta=os.path.join(config['reference_dir'], 'protein_db.mimic.fa'),
        buildSA=os.path.join(config['reference_dir'], 'buildSA.done'),
        msgf_config=lambda wildcards: DATASET_CONFIG[wildcards.dataset],
    output:
        mzid_target=os.path.join(config['work_dir'], '{dataset}/{ms_run}.mzid'),
        mzid_decoy=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.mimic.mzid'),
    group: 'msgfpercolator'
    log:
        out=os.path.join(config['log_dir'], '{dataset}/{ms_run}_msgf.log'),
    params:
        cmd=config['msgf_path'],
        mem=config['msgf_mem_gb'],
    threads: config['msgf_threads']
    resources:
        mem_mb=config['msgf_mem_gb'] * 1024,
        time_hours=config['msgf_time_hr'],
    shell:
        'echo \"Command: java -Xmx{params.mem}g -jar {params.cmd} -s {input.mzml} -o {output.mzid_target} -d {input.protein_fasta}'
        ' -addFeatures 1 -tda 0 -conf {input.msgf_config} -thread {threads}\" > {log.out}'
        ' && '
        'java -Xmx{params.mem}g -jar {params.cmd} -s {input.mzml} -o {output.mzid_target} -d {input.protein_fasta}'
        ' -addFeatures 1 -tda 0 -conf {input.msgf_config} -thread {threads} >> {log.out} 2>&1'
        ' && '
        'echo \"Command: java -Xmx{params.mem}g -jar {params.cmd} -s {input.mzml} -o {output.mzid_decoy} -d {input.decoy_fasta}'
        ' -addFeatures 1 -tda 0 -conf {input.msgf_config} -thread {threads}\" >> {log.out}'
        ' && '
        'java -Xmx{params.mem}g -jar {params.cmd} -s {input.mzml} -o {output.mzid_decoy} -d {input.decoy_fasta}'
        ' -addFeatures 1 -tda 0 -conf {input.msgf_config} -thread {threads} >> {log.out} 2>&1'


rule run_msgf2pin:
    input:
        target_mzid=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.mzid'),
        decoy_mzid=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.mimic.mzid'),
        target_fasta=os.path.join(config['reference_dir'], 'protein_db.fa'),
        decoy_fasta=os.path.join(config['reference_dir'], 'protein_db.mimic.fa'),
    output:
        pin=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.pin'),
    group: 'msgfpercolator'
    params:
        conda_wrapper=config['conda_wrapper'],
    log:
        out=os.path.join(config['log_dir'], '{dataset}/{ms_run}_msgf2pin.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} msgf2pin {input.target_mzid} {input.decoy_mzid} -v 5'
        ' -o {output.pin}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} msgf2pin {input.target_mzid} {input.decoy_mzid} -v 5'
        ' -o {output.pin}'
        ' >> {log.out} 2>&1'


rule run_percolator:
    input:
        pin=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.pin'),
    output:
        psms=temp(os.path.join(config['work_dir'], '{dataset}', '{ms_run}.temp')),
    group: 'msgfpercolator'
    log:
        out=os.path.join(config['log_dir'], '{dataset}/{ms_run}_percolator.log'),
    params:
        train_fdr=config['train_fdr'],
        test_fdr=config['test_fdr'],
        conda_wrapper=config['conda_wrapper'],
    threads: config['percolator_threads']
    resources:
        mem_mb=config['percolator_mem_gb'] * 1024,
        time_hours=config['percolator_time_hr'],
    shell:
        'echo \"Command: bash {params.conda_wrapper} percolator {input.pin}'
        ' -m {output.psms} --num-threads {threads}'
        ' --trainFDR {params.train_fdr} --testFDR {params.test_fdr}'
        ' --only-psms --search-input separate\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} percolator {input.pin}'
        ' -m {output.psms} --num-threads {threads}'
        ' --trainFDR {params.train_fdr} --testFDR {params.test_fdr}'
        ' --only-psms --search-input separate'
        ' >> {log.out} 2>&1'
        
        
rule check_percolator_output:
    input:
        psms=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.temp'),
    output:
        psms=temp(os.path.join(config['work_dir'], '{dataset}', '{ms_run}.percolator_psms')),
    group: 'msgfpercolator'
    log:
        out=os.path.join(config['log_dir'], '{dataset}/{ms_run}_percolator_check.log'),
    params:
        scripts=config['scripts_dir'],
        conda_wrapper=config['conda_wrapper'],
    threads: config['percolator_threads']
    resources:
        mem_mb=config['percolator_mem_gb'] * 1024,
        time_hours=config['percolator_time_hr'],
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/check_percolator_output.py'
        ' -i {input.psms} -l {log.out} -o {output.psms}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/check_percolator_output.py'
        ' -i {input.psms} -l {log.out} -o {output.psms}'
        ' >> {log.out} 2>&1'


rule filter_significant_psms:
    input:
        psms=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.percolator_psms'),
        pin=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.pin'),
    output:
        psms=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.sig.psms'),
        pin=os.path.join(config['work_dir'], '{dataset}', '{ms_run}.sig.pin'),
    group: 'msgfpercolator'
    log:
        out=os.path.join(config['log_dir'], '{dataset}/{ms_run}_filter_psms.log'),
    params: 
        fdr=config['qvalue_cutoff'],
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/filter_significant_psms.py'
        ' -i {input.psms} -o {output.psms} -f {params.fdr}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/filter_significant_psms.py'
        ' -i {input.psms} -o {output.psms} -f {params.fdr}'
        ' >> {log.out} 2>&1'
        ' && '
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/filter_significant_pin_by_psms.py'
        ' -i {output.psms} -p {input.pin} -o {output.pin}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/filter_significant_pin_by_psms.py'
        ' -i {output.psms} -p {input.pin} -o {output.pin}'
        ' >> {log.out} 2>&1'


def get_dataset_all_runs(wildcards):
    all_runs = ['%s/%s/%s.sig.psms'%(config['work_dir'], wildcards.dataset, ms_run) for ms_run in DATASET_ALL_MS_RUNS[wildcards.dataset]]
    return all_runs


rule combine_psms_for_a_dataset:
    input:
        psms_list=get_dataset_all_runs,
    output:
        psms=os.path.join(config['results_dir'], '{dataset}.sig_psms_count.txt'),
    log:
        out=os.path.join(config['log_dir'], '{dataset}_combine_psms.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/collapse_psms_in_all_samples.py'
        ' -i {input.psms_list} -o {output.psms}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/collapse_psms_in_all_samples.py'
        ' -i {input.psms_list} -o {output.psms}'
        ' >> {log.out} 2>&1'


def get_dataset_all_runs_intensity(wildcards):
    all_runs = ['%s/%s/%s.sig.pin'%(config['work_dir'], wildcards.dataset, ms_run) for ms_run in DATASET_ALL_MS_RUNS[wildcards.dataset]]
    return all_runs


rule combine_intensity_for_a_dataset:
    input:
        pin_list=get_dataset_all_runs_intensity,
    output:
        intensity=os.path.join(config['results_dir'], '{dataset}.sig_psms_intensity.txt'),
    log:
        out=os.path.join(config['log_dir'], '{dataset}_combine_intensity.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/collapse_intensity_in_all_pin_msgf.py'
        ' -i {input.pin_list} -o {output.intensity}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/collapse_intensity_in_all_pin_msgf.py'
        ' -i {input.pin_list} -o {output.intensity}'
        ' >> {log.out} 2>&1'


rule combine_all_psms:
    input:
        psms_list=expand(os.path.join(config['results_dir'], '{dataset}.sig_psms_count.txt'), dataset=DATASETS),
    output:
        psms=os.path.join(config['results_dir'], 'all_sig_psms_count.txt'),
    log:
        out=os.path.join(config['log_dir'], 'combine_psms.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/combine_psm_tables.py'
        ' -i {input.psms_list} -o {output.psms}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/combine_psm_tables.py'
        ' -i {input.psms_list} -o {output.psms}'
        ' >> {log.out} 2>&1'


rule combine_all_intensity:
    input:
        psms_list=expand(os.path.join(config['results_dir'], '{dataset}.sig_psms_intensity.txt'), dataset=DATASETS),
    output:
        psms=temp(os.path.join(config['results_dir'], 'all_sig_psms_intensity.txt')),
    log:
        out=os.path.join(config['log_dir'], 'combine_intensity.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/combine_psm_tables.py'
        ' -i {input.psms_list} -o {output.psms}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/combine_psm_tables.py'
        ' -i {input.psms_list} -o {output.psms}'
        ' >> {log.out} 2>&1'


rule filter_novel_peptides:
    input:
        peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.txt'),
        protein_fasta=os.path.join(config['reference_dir'], 'protein_db.fa'),
        index_table=os.path.join(config['reference_dir'], 'protein_db.index.txt'),
    output:
        uniprot_peptide_table=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.uniprot_tmp.txt')),
        gencode_peptide_table=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.gencode_tmp.txt')),
        other_canonical_peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.other_canonical.txt'),
        novel_peptide_table=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel_tmp.txt')),
    group: "filter_novel_peptides"
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'filter_novel_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/filter_novel_peptides.py'
        ' -pep {input.peptide_table} -ref {input.protein_fasta} -idx {input.index_table}'
        ' -o {output.novel_peptide_table} --output_gencode {output.gencode_peptide_table}'
        ' --output_uniprot {output.uniprot_peptide_table}'
        ' --output_canonical {output.other_canonical_peptide_table}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/filter_novel_peptides.py'
        ' -pep {input.peptide_table} -ref {input.protein_fasta} -idx {input.index_table}'
        ' -o {output.novel_peptide_table} --output_gencode {output.gencode_peptide_table}'
        ' --output_uniprot {output.uniprot_peptide_table}'
        ' --output_canonical {output.other_canonical_peptide_table}'
        ' >> {log.out} 2>&1'


rule search_novel_peptides_in_uniprot:
    input:
        peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.novel_tmp.txt'),
        uniprot_protein=get_uniprot_protein,
    output:
        mapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.uniprot.txt')),
        unmapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_uniprot.txt')),
    group: "filter_novel_peptides"
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'search_novel_peptides_in_uniprot.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/search_peptides_in_fasta.py'
        ' -pep {input.peptide_table}'
        ' -o1 {output.mapped} -o2 {output.unmapped}'
        ' -ref {input.uniprot_protein}'
        ' --header_delimiter \'|\' --header_index 1\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/search_peptides_in_fasta.py'
        ' -pep {input.peptide_table}'
        ' -o1 {output.mapped} -o2 {output.unmapped}'
        ' -ref {input.uniprot_protein}'
        ' --header_delimiter \'|\' --header_index 1'
        ' >> {log.out} 2>&1'


rule search_novel_peptides_in_gencode_pc:
    input:
        peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_uniprot.txt'),
        gencode_protein=get_gencode_protein,
    output:
        mapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.txt')),
        unmapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.txt')),
    group: "filter_novel_peptides"
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'search_novel_peptides_in_gencode_pc.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/search_peptides_in_fasta.py'
        ' -pep {input.peptide_table}'
        ' -o1 {output.mapped} -o2 {output.unmapped}'
        ' -ref {input.gencode_protein}'
        ' --header_delimiter \'|\' --header_index 1\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/search_peptides_in_fasta.py'
        ' -pep {input.peptide_table}'
        ' -o1 {output.mapped} -o2 {output.unmapped}'
        ' -ref {input.gencode_protein}'
        ' --header_delimiter \'|\' --header_index 1'
        ' >> {log.out} 2>&1'


rule get_fasta_novel_peptides:
    input:
        peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.txt'),
    output:
        fasta=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.fa')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'get_fasta_novel_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/convert_tabular_to_fasta.py'
        ' -i {input.peptide_table} -o {output.fasta}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/convert_tabular_to_fasta.py'
        ' -i {input.peptide_table} -o {output.fasta}'
        ' >> {log.out} 2>&1'


rule search_novel_peptides_in_gencode_gtf:
    input:
        fasta=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.fa'),
        gtf=get_genome_gtf,
        genome_fasta=get_genome_fasta,
    output:
        mapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.gencode_v39.gtf')),
        unmapped=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.fa')),
    group: "filter_novel_peptides"
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'search_novel_peptides_in_gencode_gtf.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/search_peptides.py -v 1'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}'
        ' --output_unmapped {output.unmapped}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/search_peptides.py -v 1'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}'
        ' --output_unmapped {output.unmapped}'
        ' >> {log.out} 2>&1'


rule search_non_gencode_novel_peptides:
    input:
        fasta=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.fa'),
        gtf=get_rna_gtf,
        genome_fasta=get_genome_fasta,
    output:
        mapped=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.rna_sample.gtf'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'search_non_gencode_novel_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/search_peptides.py'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/search_peptides.py'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}'
        ' >> {log.out} 2>&1'


rule get_fasta_gencode_pc_peptides:
    input:
        peptide_table=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.txt'),
    output:
        fasta=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.fa'),
    log:
        out=os.path.join(config['log_dir'], 'search_gencode_pc_peptides.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/convert_tabular_to_fasta.py'
        ' -i {input.peptide_table} -o {output.fasta}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/convert_tabular_to_fasta.py'
        ' -i {input.peptide_table} -o {output.fasta}'
        ' >> {log.out} 2>&1'


rule search_gencode_pc_peptides:
    input:
        fasta=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.fa'),
        gtf=get_genome_gtf,
        genome_fasta=get_genome_fasta,
    output:
        mapped=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.gencode_v39.gtf'),
        unmapped=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.no_v39_mapping.fa'),
    log:
        out=os.path.join(config['log_dir'], 'search_gencode_pc_peptides.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/search_peptides.py'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}'
        ' --output_unmapped {output.unmapped}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/search_peptides.py'
        ' -i {input.fasta} -c {input.gtf} -g {input.genome_fasta} -o {output.mapped}'
        ' --output_unmapped {output.unmapped}'
        ' >> {log.out} 2>&1'


rule select_novel_not_gencode_mapped_peptides:
    input:
        novel_all=os.path.join(config['results_dir'], 'all_sig_psms_count.novel_tmp.txt'),
        novel_not_gencode=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.fa'),
    output:
        novel_not_gencode=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.txt'),
    log:
        out=os.path.join(config['log_dir'], 'select_novel_not_gencode_mapped_peptides.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/select_peptides_by_fasta.py'
        ' -i {input.novel_all} -f {input.novel_not_gencode} -o {output.novel_not_gencode}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/select_peptides_by_fasta.py'
        ' -i {input.novel_all} -f {input.novel_not_gencode} -o {output.novel_not_gencode}'
        ' >> {log.out} 2>&1'


rule select_novel_gencode_mapped_peptides:
    input:
        novel_all=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.txt'),
        novel_not_gencode=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.fa'),
    output:
        novel_gencode=temp(os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.gencode_v39.txt')),
    log:
        out=os.path.join(config['log_dir'], 'select_novel_gencode_mapped_peptides.log'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/select_peptides_by_fasta.py'
        ' -i {input.novel_all} -f {input.novel_not_gencode} -o {output.novel_gencode} --retaining False\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/select_peptides_by_fasta.py'
        ' -i {input.novel_all} -f {input.novel_not_gencode} -o {output.novel_gencode} --retaining False'
        ' >> {log.out} 2>&1'


rule combine_uniprot_peptides:
    input:
        uniprot=os.path.join(config['results_dir'], 'all_sig_psms_count.uniprot_tmp.txt'),
        novel_uniprot=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.uniprot.txt'),
    output:
        uniprot=os.path.join(config['results_dir'], 'all_sig_psms_count.uniprot.txt'),
    log:
        out=os.path.join(config['log_dir'], 'combine_uniprot_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cat {input.uniprot} > {output.uniprot}\" > {log.out}'
        ' && '
        'echo \"Command: tail -n+2 {input.novel_uniprot} >> {output.uniprot}\" >> {log.out}'
        ' && '
        'cat {input.uniprot} > {output.uniprot}'
        ' && '
        'tail -n+2 {input.novel_uniprot} >> {output.uniprot}'


rule combine_gencode_peptides:
    input:
        gencode=os.path.join(config['results_dir'], 'all_sig_psms_count.gencode_tmp.txt'),
        novel_gencode_pc=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.gencode_pc_translation.txt'),
        novel_gencode_predict=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode_pc_translation.gencode_v39.txt'),
    output:
        gencode=os.path.join(config['results_dir'], 'all_sig_psms_count.gencode.txt'),
    log:
        out=os.path.join(config['log_dir'], 'combine_gencode_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: cat {input.gencode} > {output.gencode}\" > {log.out}'
        ' && '
        'echo \"Command: tail -n+2 {input.novel_gencode_pc} >> {output.gencode}\" >> {log.out}'
        ' && '
        'echo \"Command: tail -n+2 {input.novel_gencode_predict} >> {output.gencode}\" >> {log.out}'
        ' && '
        'cat {input.gencode} > {output.gencode}'
        ' && '
        'tail -n+2 {input.novel_gencode_pc} >> {output.gencode}'
        ' && '
        'tail -n+2 {input.novel_gencode_predict} >> {output.gencode}'


rule plot_novel_peptides:
    input:
        peptide_gtf=os.path.join(config['results_dir'], 'all_sig_psms_count.novel.not_gencode.rna_sample.gtf'),
        ref_gtf=get_rna_gtf,
    output:
        touch(os.path.join(config['visualization_dir'], 'plot_novel_peps.done')),
        folder=directory(config['visualization_dir']),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'plot_novel_peptides.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/plot_peptides_to_tx.py'
        ' -i {input.peptide_gtf} -r {input.ref_gtf} -o {output.folder}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/plot_peptides_to_tx.py'
        ' -i {input.peptide_gtf} -r {input.ref_gtf} -o {output.folder}'
        ' >> {log.out} 2>&1'


rule predict_gene_abun_by_psm:
    input:
        psm=os.path.join(config['results_dir'], 'all_sig_psms_count.txt'),
        index=os.path.join(config['reference_dir'], 'protein_db.index.txt'),
    output:
        gene_table=os.path.join(config['results_dir'], 'all_gene_level_count.txt'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_gene_abun_by_psm.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/get_gene_level_exp_by_EM.py'
        ' -p {input.psm} -i {input.index} -o {output.gene_table}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/get_gene_level_exp_by_EM.py'
        ' -p {input.psm} -i {input.index} -o {output.gene_table}'
        ' >> {log.out} 2>&1'


rule convert_psm_intensity_log2:
    input:
        psm=os.path.join(config['results_dir'], 'all_sig_psms_intensity.txt'),
    output:
        psm=os.path.join(config['results_dir'], 'all_sig_psms_intensity.log2.txt'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'combine_intensity.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/convert_expression_to_log2.py'
        ' -i {input.psm} -o {output.psm} --columns 1,2\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/convert_expression_to_log2.py'
        ' -i {input.psm} -o {output.psm} --columns 1,2'
        ' >> {log.out} 2>&1'


rule predict_gene_abun_by_intensity:
    input:
        psm=os.path.join(config['results_dir'], 'all_sig_psms_intensity.txt'),
        index=os.path.join(config['reference_dir'], 'protein_db.index.txt'),
    output:
        gene_table=temp(os.path.join(config['results_dir'], 'all_gene_level_intensity.txt')),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_gene_abun_by_psm.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/get_gene_level_exp_by_EM.py'
        ' -p {input.psm} -i {input.index} -o {output.gene_table}\"'
        ' > {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/get_gene_level_exp_by_EM.py'
        ' -p {input.psm} -i {input.index} -o {output.gene_table}'
        ' >> {log.out} 2>&1'


rule convert_gene_intensity_log2:
    input:
        gene_table=os.path.join(config['results_dir'], 'all_gene_level_intensity.txt'),
    output:
        gene_table=os.path.join(config['results_dir'], 'all_gene_level_intensity.log2.txt'),
    params:
        conda_wrapper=config['conda_wrapper'],
        scripts=config['scripts_dir'],
    log:
        out=os.path.join(config['log_dir'], 'predict_gene_abun_by_intensity.log'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'echo \"Command: bash {params.conda_wrapper} python {params.scripts}/convert_expression_to_log2.py'
        ' -i {input.gene_table} -o {output.gene_table} --columns 1\"'
        ' >> {log.out}'
        ' && '
        'bash {params.conda_wrapper} python {params.scripts}/convert_expression_to_log2.py'
        ' -i {input.gene_table} -o {output.gene_table} --columns 1'
        ' >> {log.out} 2>&1'

