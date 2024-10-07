# ProteoSeq
Workflow of integrative proteotranscriptomic data analysis.

## Table of Contents
- [1. Dependencies](#1-dependencies)
- [2. Test installation and run ProteoSeq](#2-test-installation-and-run-proteoseq)
- [3. Input files and pre-processing the data](#3-input-files-and-pre-processing-the-data)
  - [3.1 Mass spectrometry data](#31-mass-spectrometry-data)
  - [3.2 Use RNA-seq results to construct protein database](#32-use-rna-seq-results-to-construct-protein-database)
  - [3.3 Use protein sequences to construct protein database](#33-use-protein-sequences-to-construct-protein-database)
- [4. Configuration](#4-configuration)
  - [4.1 Input the mzML files of MS proteomics data](#41-input-the-mzml-files-of-ms-proteomics-data)
  - [4.2 RNA GTF file for sample-specific protein database](#42-rna-gtf-file-for-sample-specific-protein-database)
  - [4.3 Protein sequences for sample-specific database](#43-protein-sequences-for-sample-specific-database)
  - [4.4 Genome reference](#44-genome-reference)
  - [4.5 Canonical reference](#45-canonical-reference)
  - [4.6 Method to construct decoy database](#46-method-to-construct-decoy-database)
- [5. Outputs and visualization](#5-outputs-and-visualization)
  - [5.1 References](#51-references)
  - [5.2 Outputs on peptides level](#52-outputs-on-peptides-level)
  - [5.3 Outputs on proteins level](#53-outputs-on-proteins-level)
  - [5.4 Visualization of novel peptides](#54-visualization-of-novel-peptides)
  - [5.5 Obtain mappings and visualization of given peptides](#55-obtain-mappings-and-visualization-of-given-peptides)


## 1. Dependencies
ProteoSeq requires the following to be installed and available on $PATH:
- [Python 3](https://www.python.org/)
  - [NumPy](http://numpy.org/)
  - [SciPy](https://anaconda.org/cctbx202208/scipy) >= 1.9.0
  - [Biopython](https://anaconda.org/cctbx202112/biopython) >= 1.7.9
  - [Matplotlib](https://matplotlib.org/)
- [Percolator](https://anaconda.org/bioconda/percolator) == 3.05.0
- [MS-GF+](https://github.com/MSGFPlus/msgfplus)
- [Java](https://www.java.com)

If using the Snakemake then its installation script can install a conda environment to $HOME path with the dependencies.
```
./install
```

ProteoSeq uses [**MS-GF+**](https://msgfplus.github.io) as the search engine to identify peptides in LC-MS/MS datasets. We recommend downloading it from: https://github.com/MSGFPlus/msgfplus/releases/download/v2021.09.06/MSGFPlus_v20210906.zip. Decompress it and set up the path to `MSGFPlus.jar` in the configuration file `snakemake_config.yaml`. Java is required for running MS-GF+, for more information, please refer to the [documentation](https://msgfplus.github.io/msgfplus/) of MS-GF+.
```
mkdir msgf_plus
cd msgf_plus
curl -L 'https://github.com/MSGFPlus/msgfplus/releases/download/v2021.09.06/MSGFPlus_v20210906.zip' -O
unzip MSGFPlus_v20210906.zip
cd ..
```

ProteoSeq requires a **GTF file** (for transcripts annotation) or a **FASTA file** (for protein sequences) as input to construct the sample-specific protein database. If the data is in another format then other tools can be used to create the transcript annotation file: [**ESPRESSO**](https://github.com/Xinglab/espresso)

# 2. Test installation and run ProteoSeq
A small test data set is provided in `test_data.tar.gz`. Run `tar zxvf test_data.tar.gz` to decompress it. The unpacked files are:
- `rna_seq.gtf`: annotation file including a small list of novel transcripts discovered in long-read RNA-seq data. This file can be used as input to `RNA_gtf` to search novel peptides encoded by transcripts in the file.
- `novel_protein.fasta`: a file including a small list of novel protein sequences. This file can be used as input to `protein_fasta` to search novel peptides encoded by proteins in the file.
- `ms_dataset_1` & `ms_dataset_2`: two MS dataset including a `msgf_config.txt` for each dataset and mzML folder containing one to two mzML files.

The `snakemake_config.yaml` uses the paths for the test data by default. ProteoSeq can be run in a command line by:
```
./run
```

## 3. Input files and pre-processing the data
### 3.1 Mass spectrometry data
ProteoSeq accepts DDA-MS data. The input MS data must be centroided **mzML** files. If you have the vendor-format data files, please convert them to mzML format. [ProteoWizard msConvert](https://proteowizard.sourceforge.io/index.html) may be used to convert other formats to **mzML** format. Parameters (peptide length, protease used in the MS experiment, modifications, etc.) need to be specified in a `msgf_config.txt` file for each MS dataset. Example files can be find in test_data.
### 3.2 Use RNA-seq results to construct protein database
ProteoSeq can use RNA-seq data to guide building sample-specific protein database. It adopts the standard [**GTF**](https://useast.ensembl.org/info/website/upload/gff.html) file as input. The GTF file must contain at least "**transcripts**" and "**exons**" features. As default, the GTF file output from most RNA-seq processing tools (e.g., [ESPRESSO](https://github.com/Xinglab/espresso)) can be directly used by setting the path to `RNA_gtf` in the configuration file `snakemake_config.yaml`. We also provide some scripts that are not embedded in the `snakemake` pipeline but can be run in a command line to pre-processing the GTF to remove low-confidence transcripts.
#### 3.2.1 Filter out low-expressing transcripts
We provide a script to filter out low expressing transcripts at certain threshold. The expression table should be in tabular text format, with each row represents a transcript and each column represents a sample. Values should be numbers representing the expression level of a given transcript in each sample. Specify sample in the script to filter by the expression in a given sample. By default (--sample All), the script uses maximum expression (--value max) in all samples for a given transcript. This will select the transcripts expressing at a certain cutoff in at least one sample.
```
EXPRESSION_TABLE=rna_expression.txt
## Select transcripts present in at least one sample
# python scripts/filter_gtf_by_cpm.py -i ${SAMPLE_GTF} -e ${EXPRESSION_TABLE} -o rna_sample_0.gtf --sample ${SAMPLE} --cutoff 0 --value max --include_equal False
## Select transcripts present at certain cutoff in at least one sample
python scripts/filter_gtf_by_cpm.py -i ${SAMPLE_GTF} -e ${EXPRESSION_TABLE} -o rna_sample_001.gtf --sample ${SAMPLE} --cutoff 0.01 --value max --include_equal True
```
#### 3.2.2 Select transcripts in certain chromosomes
We provide a script that can select transcripts in interested chromosomes from the GTF. This step is recommended because the default output GTF of many RNA-seq processing tools may contain transcripts in scaffolds and unlocalized/unplaced sequences. By default, the script selects transcripts in chr1-22, X&Y. <br />
```
python scripts/filter_chromosomes_in_gtf.py -i ${SAMPLE_GTF} -o rna_sample.primary_chrom.gtf`
```
### 3.3 Use protein sequences to construct protein database
ProteoSeq can directly accept protein sequences (`protein_fasta`) to construct sample-specific protein database. It adopts the standard **FASTA** file as input. Peptides that do not map to canonical protein sequences (see **4.5** for the definition of canonical sequences) and map to protein sequences in this file will be reported as novel peptides.

## 4. Configuration
### 4.1 Input the mzML files of MS proteomics data
ProteoSeq can analyze various DDA-MS datasets at a same time. Each dataset should have a `msgf_config` file to specify its spectra type, MS instrument and experimental protocols such as the proteases (e.g., trypsin) used in the MS experiments. Raw MS files need to be converted to centroided `mzML` format and put in a separate folder for each dataset.
```
ms_datasets:
  dataset_1_name:
    mzml_dir: 'test_data/ms_dataset_1/mzML'
    msgf_config: 'test_data/ms_dataset_1/msgf_config.txt'
  dataset_2_name:
    mzml_dir: 'test_data/ms_dataset_2/mzML'
    msgf_config: 'test_data/ms_dataset_2/msgf_config.txt'
```
### 4.2 RNA GTF file for sample-specific protein database
ProteoSeq can use a single GTF file processed from RNA-seq data as the input to guide the construction of sample-specific searching space. The GTF file should include transcript annotations including **transcript** and **exon** coordinates. The absolute path to the GTF file needs to be set up in `snakemake_config.yaml`. An example file can be found in `test_data`.
<br/>
If CDS is not annotated for a transcript in the GTF file and if `predict_ORF` is on, ProteoSeq will look for the longest in-frame ORF in the transcript sequence and translate protein sequence. 
<br/>
ProteoSeq also predicts if a transcript is targeted to **NMD pathway** with three [rules](https://pubmed.ncbi.nlm.nih.gov/27618451/): 1. > one exon; 2. CDS >= 150nt; 3. lastEJC - stop_codon >= 50nt. By default, NMD transcripts are discarded. If `discard_NMD` is set to false, translation product of NMD transcripts will be included in the protein database.
```
RNA_gtf: '/path/to/rna_seq.gtf'
predict_ORF: true
discard_NMD: true
```

### 4.3 Protein sequences for sample-specific database
ProteoSeq can accept protein sequences to construct sample-specific protein database. The absolute path to the FASTA file needs to be set up in `snakemake_config.yaml`. An example file can be found in `test_data`.
```
protein_fasta: '/path/to/novel_protein.fasta'
```

### 4.4 Genome reference
A FASTA file of genome sequences and a GTF file of GENCODE annotation is required for ProteoSeq. If the paths are not set up in the configuration file, ProteoSeq will automatically download `GRCh38.primary_assembly.genome.fa` and `gencode.v39.annotation.gtf`.
```
GENCODE_gtf_path: '/path/to/GENCODE.gtf' #Default if empty
genome_fasta_path: '/path/to/genome.fa' #Default if empty
```
Default downloads include:
```
GENCODE_gtf_url: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz'
genome_fasta_url: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz'
```

### 4.5 Canonical reference
An important feature of ProteoSeq is that for all identified peptides, ProteoSeq scans their mapping to parental protein and transcript isoforms at a single amino acid resolution. Therefore, it can rapidly classify peptides into different categories and uncover high-confidence novel peptides. 
<br/>
First, ProteoSeq will scan all peptides in UniProt sequences. If `UniProt_protein_path` is not set, ProteoSeq will automatically download the UniProt reviewed sequences plus isoforms from `UniProt_protein_url`. Human protein sequences (OX=9606) will be selected and sequences with ambiguous amino acid ("X") will be discarded. 
<br/>
By default, the UniProt protein sequences will be combined to proteins in sample-specific protein database for analyzing the MS data, unless `merge_UniProt_protein` is set to false.
```
merge_UniProt_protein: true
UniProt_protein_path: ''
UniProt_protein_url: 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2021_04/knowledgebase/uniprot_sprot-only2021_04.tar.gz'
```
After filtering out peptides mapping to UniProt sequences, ProteoSeq will scan other peptides in GENCODE annotated sequences. If `GENCODE_protein_path` is not set, ProteoSeq will automatically download the GENCODE proteins from `GENCODE_protein_url`. 
By default, the GENCODE proteins are not included in sample-specific protein database, to avoid inflating the searching space. Users can turn on `merge_GENCODE_protein` to include them in the protein database.
```
merge_GENCODE_protein: false
GENCODE_protein_path: ''
GENCODE_protein_url: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.pc_transcripts.fa.gz'
```

Users can specify their own canonical protein sequences. Peptides mapping to these provided canonical protein sequences will be reported as a separate file and will not be considered as novel peptides.
```
customized_canonical_protein_path: ''
```
### 4.6 Method to construct decoy database
The target-decoy search strategy is commonly used in proteomics to increase confidence in peptide identification. Each mass spectrum is searched against both a target protein database and a decoy database, generating two lists of peptide-spectrum matches (PSMs). These PSMs are then analyzed using [Percolator](http://percolator.ms)  to statistically estimate confidence levels. Common methods for constructing decoy databases include reversing or shuffling the protein sequences. ProteoSeq offers multiple options for decoy database construction and can run Percolator to identify significant peptides based on user-defined significance thresholds.
```
# Decoy database construction method
# If both are on, ProteoSeq will first reverse each sequence and then shuffle the amino acids in each sequence
reverse: false
shuffle: true
# Cutoffs for Percolator
train_fdr: 0.05
test_fdr: 0.05
qvalue_cutoff: 0.05
```
## 5. Outputs and visualization
### 5.1 References
Reference files are retained in the `reference_dir`. 

### 5.2 Outputs on peptides level
A few files will be output to `results_dir`.
- **all_sig_psms_count.txt**: All identified peptides and the count of significant PSMs with given FDR cutoff. Peptides are classified into three groups: UniProt, GENCODE and Novel.
- **all_sig_psms_count.uniprot.txt**: Peptides map to UniProt protein sequences.
- **all_sig_psms_count.gencode.txt**: Peptides which do not map to UniProt sequences and are encoded by GENCODE annotated proteins or transcripts.
- **all_sig_psms_count.other_canonical.txt**: Peptides which do not fall into the previous two categories and map to proteins provided in `customized_canonical_protein_path`. 
- **all_sig_psms_count.novel.txt**: Peptides which do not fall into the previous categories and map to transcripts provided in `RNA_gtf` or proteins provided in `protein_fasta`.

### 5.3 Outputs on proteins level
ProteoSeq by default runs an EM algorithm with spectral counts or total MS2 ion intensities  to predict the abundance of proteins. Outputs include:
- **all_gene_level_count.txt**: Gene-level protein abundance in each MS sample, predicted from spectral counts (set `predict_gene_by_PSM_count: true` to enable this prediction).
- **all_gene_level_intensity.log2.txt**: Gene-level protein abundance in each MS sample, predicted from log2 transformed total MS2 ion intensities (set `predict_gene_by_intensity: true` to enable this prediction).

### 5.4 Visualization of novel peptides
By default, ProteoSeq will plot the mappings of all novel peptides (defined as in **5.2**) to the transcripts provided in the 'RNA_gtf' file. For each of the genes with novel peptides, a PDF image will be output to `visualization_dir`. Set `enable_visualization:` to false in the configuration file to skip this step.

### 5.5 Obtain mappings and visualization of given peptides
Visualization can also be run manually with built-in python scripts for any given peptide sequences. It includes two steps: (1) given peptides are searched to obtain their mappings to a given list of transcripts. (2) the mappings are plotted to each gene.

In the first step, given peptides need to be stored in a FASTA format file, each peptide needs to have an ID line starting with ">" and a sequence line. The searching space can either be the sample-specific `rna_seq.gtf` generated by ProteoSeq or the whole annotated transcriptome provided by GENCODE (e.g., `gencode.v39.annotation.gtf`. The GTF file must contain at least the coordinates of transcripts and exons. If the CDS coordinates are not found in the GTF file, ProteoSeq will automatically predict the longest ORF in the transcripts. A FASTA file with genome sequences is also required.

```
# python scripts/search_peptides.py -i peptides.fa -c ${GENCODE_GTF} -g ${GENOME_FA} -o peptides_in_gencode_tx.gtf 
python scripts/search_peptides.py -i peptides.fa -c rna_seq.gtf -g ${GENOME_FA} -o ${WORK_DIR}/peptides_in_sample_specific_tx.gtf
```

In the second step, the mapping of a peptide to transcripts are plotted. For each gene, a high-resolution PDF image will be output. The coordinates of peptides (these can be predicted by running search_peptides.py), coordinates of exons and transcripts (the same GTF file for running search_peptides.py) are required for visualization. 

```
# python scripts/plot_peptides_to_tx.py -i peptides_in_gencode_tx.gtf -r ${GENCODE_GTF} -o visualization_dir
python scripts/plot_peptides_to_tx.py -i peptides_in_sample_specific_tx.gtf -r rna_seq.gtf -o visualization_dir
```
