### rMAP-Candida: Rapid Mycology Analysis Pipeline for *Candida*

**rMAP-Candida** is a reproducible, Dockerized bioinformatics workflow for the analysis of paired-end Illumina sequencing data from clinical *Candida* isolates. The pipeline is designed to support fungal genomic surveillance, species confirmation, quality control, genome assembly, annotation, antifungal resistance characterization, virulence assessment & integrated HTML reporting.

The workflow is implemented in **WDL** & is intended to run with **Cromwell**, making it portable across local workstations, institutional servers & cloud-based environments.

---

## Overview

Clinical *Candida* infections are increasingly important in hospital & public health settings because of rising antifungal resistance, species diversity & the emergence of multidrug-resistant pathogens such as *Candida auris*. Accurate genomic characterization of *Candida* isolates can support outbreak investigation, antifungal resistance surveillance, infection prevention & research into pathogen evolution.

**rMAP-Candida** provides an end-to-end analysis framework for fungal paired-end genome sequencing data, with emphasis on:

- Reproducibility
- Clinical isolate analysis
- Species-level identification
- Antifungal resistance gene & mutation detection
- Assembly quality assessment
- Genome annotation
- Summary reporting for surveillance

---

## Key features

- Paired-end Illumina FASTQ input support
- Adapter trimming & read preprocessing
- Read-level quality control using FastQC
- Aggregated QC reporting using MultiQC
- *Candida* species typing using Kraken2/Bracken or a fungal-focused database
- Genome assembly using SPAdes or a fungal genome assembler
- Assembly quality assessment using QUAST
- Fungal genome completeness assessment using BUSCO
- Genome annotation using tools such as Funannotate, Prokka, or related fungal annotation approaches
- Antifungal resistance characterization
- Virulence-associated gene screening
- Optional phylogenomic analysis for closely related isolates
- Integrated sample-level & cohort-level HTML report
- Dockerized execution for reproducibility
- WDL/Cromwell workflow structure for portability & scalability

---

## Intended use

This pipeline is intended for research & genomic surveillance applications involving *Candida* spp., including but not limited to:

- *Candida albicans*
- *Candida auris*
- *Candida glabrata*
- *Candida parapsilosis*
- *Candida tropicalis*
- *Candida krusei*
- *Candida dubliniensis*
- Other clinically relevant *Candida* species represented in the reference databases

---

## Example applications

The pipeline can be used for:

1. **Clinical isolate genomic characterization**
   - Identify the likely *Candida* species
   - Assess read & assembly quality
   - Generate annotated genome outputs

2. **Antifungal resistance surveillance**
   - Screen for known resistance-associated mutations or genes
   - Summarize predicted resistance markers
   - Compare resistance profiles across samples

3. **Outbreak and transmission investigation**
   - Compare genomes from multiple isolates
   - Identify clustering patterns
   - Support phylogenomic interpretation

4. **Public health mycology research**
   - Characterize *Candida* diversity
   - Generate preliminary data for manuscripts & grants
   - Support local fungal genomic surveillance capacity

---

## Workflow summary

The pipeline performs the following major steps:

```text
Input paired-end FASTQ files
        |
        v
Read quality control
        |
        v
Adapter trimming & read preprocessing
        |
        v
Post-trimming QC & MultiQC report
        |
        v
Candida species typing
        |
        v
Genome assembly
        |
        v
Assembly quality assessment
        |
        v
BUSCO completeness assessment
        |
        v
Genome annotation
        |
        v
Antifungal resistance & virulence screening
        |
        v
Optional phylogenomic analysis
        |
        v
Integrated HTML report
```

---

## Pipeline modules

### 1. Read preprocessing

Raw paired-end FASTQ reads are processed to remove adapters, trim low-quality bases & discard poor-quality reads.

Typical tools:

- Trimmomatic
- fastp
- Cutadapt

Outputs include:

- Trimmed FASTQ files
- Read count summaries
- Trimming logs

---

### 2. Quality control

The pipeline evaluates sequencing quality before & after trimming.

Typical tools:

- FastQC
- MultiQC

Outputs include:

- Per-sample FastQC HTML reports
- Aggregated MultiQC report
- QC summary table

---

### 3. Species identification

The pipeline performs taxonomic classification to identify the most likely *Candida* species represented in each sample.

Typical tools:

- Kraken2
- Bracken
- Custom fungal or *Candida*-focused Kraken2 database

Outputs include:

- Kraken2 classification report
- Bracken abundance table
- Most probable species per sample
- Species typing summary table

Example summary:

| Sample ID | Most probable species | Evidence | Interpretation |
|---|---|---|---|
| Sample_001 | *Candida albicans* | High read support | Species supported |
| Sample_002 | *Candida auris* | High read support | Priority pathogen |
| Sample_003 | *Candida glabrata* | Moderate read support | Species supported |

---

### 4. Genome assembly

High-quality reads are assembled into draft genomes.

Typical tools:

- SPAdes
- SKESA
- Shovill
- Unicycler, where appropriate

Outputs include:

- Draft assembly FASTA
- Assembly graph files, if available
- Assembly logs

---

### 5. Assembly quality assessment

The draft assembly is assessed for contiguity, length, GC content & other assembly-level quality metrics.

Typical tools:

- QUAST

Outputs include:

- N50
- Number of contigs
- Total assembly length
- Largest contig
- GC content
- QUAST HTML report

Example output table:

| Sample ID | Contigs | N50 | Total length | GC % |
|---|---:|---:|---:|---:|
| Sample_001 | 126 | 245,000 | 14.2 Mb | 33.5 |
| Sample_002 | 214 | 118,000 | 12.4 Mb | 44.8 |

---

### 6. Genome completeness assessment

The pipeline assesses fungal genome completeness using BUSCO.

Typical tools:

- BUSCO
- Fungal lineage datasets such as `saccharomycetes_odb10`, depending on the target species

Outputs include:

- BUSCO short summary
- BUSCO full table
- Completeness score
- Fragmented & missing gene counts

Example interpretation:

| Sample ID | Complete BUSCOs | Fragmented | Missing | Interpretation |
|---|---:|---:|---:|---|
| Sample_001 | 96.3% | 1.2% | 2.5% | Good-quality assembly |
| Sample_002 | 84.5% | 5.5% | 10.0% | Review assembly quality |

---

### 7. Genome annotation

The assembled genomes can be annotated to identify genes, coding sequences & relevant genomic features.

Possible tools:

- Funannotate
- Prokka
- Augustus
- GeneMark-ES
- InterProScan, if available

Outputs include:

- Annotated GFF
- Protein FASTA
- CDS FASTA
- GenBank file
- Annotation summary table

---

### 8. Antifungal resistance characterization

The pipeline screens for known antifungal resistance-associated genes, mutations, or genomic markers.

Potential targets include genes associated with resistance to:

- Azoles
- Echinocandins
- Polyenes
- Flucytosine

Important resistance-associated genes include:

| Gene | Relevance |
|---|---|
| `ERG11` | Azole resistance |
| `TAC1` | Azole resistance regulation |
| `MRR1` | Multidrug resistance regulation |
| `CDR1` | Efflux-mediated azole resistance |
| `CDR2` | Efflux-mediated azole resistance |
| `MDR1` | Multidrug transporter |
| `FKS1` | Echinocandin resistance |
| `FKS2` | Echinocandin resistance, especially in *C. glabrata* |
| `ERG3` | Azole & polyene-related pathways |
| `FCY1` | Flucytosine resistance |
| `FUR1` | Flucytosine resistance |

Outputs include:

- Resistance gene screening results
- Mutation summary table
- Antifungal resistance marker report
- Per-sample predicted resistance interpretation, where supported

Example table:

| Sample ID | Species | Gene | Mutation / marker | Drug class | Interpretation |
|---|---|---|---|---|---|
| Sample_001 | *C. albicans* | ERG11 | Y132F | Azole | Resistance-associated marker detected |
| Sample_002 | *C. auris* | FKS1 | S639F | Echinocandin | Resistance-associated marker detected |
| Sample_003 | *C. glabrata* | FKS2 | None detected | Echinocandin | No known marker detected |

> Note: Resistance prediction should be interpreted cautiously & ideally confirmed with phenotypic antifungal susceptibility testing where available.

---

### 9. Virulence-associated gene screening

The pipeline optionally screen for virulence-associated genes or gene families relevant to *Candida* pathogenicity.

Potential targets include:

| Gene family | Relevance |
|---|---|
| `ALS` family | Adhesion |
| `SAP` family | Secreted aspartyl proteases |
| `HWP1` | Hyphal wall protein |
| `PLB` family | Phospholipases |
| `LIP` family | Lipases |
| Biofilm-associated genes | Biofilm formation & persistence |

Outputs include:

- Virulence gene presence/absence table
- Per-sample virulence summary
- Heatmap-ready matrix

---

### 10. Optional phylogenomic analysis

For datasets containing multiple isolates of the same or closely related species, the pipeline may support phylogenomic comparison.

Possible tools:

- Snippy
- MUMmer
- Parsnp
- IQ-TREE2
- FastTree
- MAFFT
- SNP-sites

Outputs include:

- Core genome alignment
- SNP matrix
- Maximum-likelihood phylogenetic tree
- Tree image
- Pairwise SNP distance matrix

---

### 11. Integrated HTML report

The final output is an integrated HTML report that summarizes the full analysis in a user-friendly format.

The report include:

1. Run summary
2. Sample metadata
3. Read QC & trimming summary
4. Species typing results
5. Assembly quality metrics
6. BUSCO completeness results
7. Genome annotation summary
8. Antifungal resistance markers
9. Virulence-associated gene screening
10. Optional phylogenomic outputs
11. Pipeline provenance and software versions

---

## Repository structure

Recommended repository layout:

```text
rMAP_Candida/
├── README.md
├── LICENSE
├── rMAP_Candida.wdl
├── inputs/
│   └── rMAP_Candida.inputs.example.json
├── metadata/
│   └── samples.tsv
├── docker/
│   ├── Dockerfile.candida_kraken2_bracken
│   ├── Dockerfile.amr
│   └── Dockerfile.report
├── scripts/
│   ├── parse_fastqc.py
│   ├── parse_bracken.py
│   ├── parse_busco.py
│   ├── parse_amr_results.py
│   ├── merge_candida_reports.py
│   └── generate_html_report.py
├── databases/
│   └── README.md
├── docs/
│   ├── index.html
│   ├── assets/
│   └── reports/
├── test_data/
│   └── README.md
└── examples/
    ├── example_samples.tsv
    └── example_inputs.json
```

---

## Requirements

### Software

The pipeline requires:

- Java 8 or later
- Cromwell
- Docker
- Git

Recommended:

- At least 8 CPU cores
- At least 16 GB RAM
- At least 100 GB free disk space, depending on number of samples & database size

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/gmboowa/rMAP-Candida.git
cd rMAP-Candida
```

### 2. Download Cromwell

```bash
wget https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar
```

Or place your existing Cromwell JAR file somewhere accessible, for example:

```bash
~/cromwell<version>.jar
```

### 3. Confirm Docker is running

```bash
docker --version
docker ps
```

On macOS with Colima:

```bash
colima start --cpu 8 --memory 16 --disk 100
docker ps
```

---

## Input requirements

### FASTQ files

The pipeline expects paired-end FASTQ files.

Accepted formats:

```text
.fastq
.fastq.gz
.fq
.fq.gz
```

Example:

```text
Sample_001_R1.fastq.gz
Sample_001_R2.fastq.gz
Sample_002_R1.fastq.gz
Sample_002_R2.fastq.gz
```

---

## Sample metadata file

A metadata file is recommended for tracking samples.

Example `samples.tsv`:

```tsv
sample_id	read1	read2	collection_date	source
Sample_001	~/Sample_001_R1.fastq.gz	~/Sample_001_R2.fastq.gz	2025-06-01	Blood culture
Sample_002	~/Sample_002_R1.fastq.gz	~/Sample_002_R2.fastq.gz	2025-06-03	Urine
Sample_003	~/Sample_003_R1.fastq.gz	~/Sample_003_R2.fastq.gz	2025-06-05	Blood culture
```

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `read1` | Path to forward FASTQ file |
| `read2` | Path to reverse FASTQ file |

Optional columns:

| Column | Description |
|---|---|
| `location` | Sampling location |
| `collection_date` | Date of sample collection |
| `source` | Clinical source or specimen type |
| `species_expected` | Expected species, if known |
| `phenotype` | Antifungal susceptibility result, if available |

---

## Example input JSON

Example `inputs/rMAP-Candida.inputs.example.json`:

```json
{
  "rMAP_Candida.samples_tsv": "~/samples.tsv",
  "rMAP_Candida.output_prefix": "candida_analysis",
  "rMAP_Candida.outdir": "~/output",

  "rMAP_Candida.threads": 8,
  "rMAP_Candida.memory_gb": 16,

  "rMAP_Candida.kraken2_db": "~/candida_kraken2_db",
  "rMAP_Candida.busco_lineage": "saccharomycetes_odb10",

  "rMAP_Candida.run_trimming": true,
  "rMAP_Candida.run_species_typing": true,
  "rMAP_Candida.run_assembly": true,
  "rMAP_Candida.run_busco": true,
  "rMAP_Candida.run_annotation": true,
  "rMAP_Candida.run_amr": true,
  "rMAP_Candida.run_virulence": true,
  "rMAP_Candida.run_phylogeny": false
}
```

---

## Running the pipeline

### Basic run

```bash
java -jar ~/cromwell-<version>.jar run rMAP_Candida.wdl \
  --inputs inputs/rMAP-Candida.inputs.example.json
```

Example:

```bash
java -jar ~/cromwell-92.jar run ~/Candida/rMAP-Candida.wdl \
  --inputs ~/rMAP-Candida.inputs.example.json
```

---

## Output directory

The final output directory contain:

```text
output/
├── qc/
│   ├── fastqc/
│   ├── multiqc/
│   └── trimming_summary.tsv
├── species_typing/
│   ├── kraken2/
│   ├── bracken/
│   └── species_summary.tsv
├── assembly/
│   ├── assemblies/
│   └── assembly_summary.tsv
├── quast/
│   ├── reports/
│   └── quast_summary.tsv
├── busco/
│   ├── reports/
│   └── busco_summary.tsv
├── annotation/
│   ├── gff/
│   ├── proteins/
│   ├── cds/
│   └── annotation_summary.tsv
├── antifungal_resistance/
│   ├── resistance_markers.tsv
│   ├── resistance_summary.tsv
│   └── resistance_interpretation.tsv
├── virulence/
│   ├── virulence_genes.tsv
│   └── virulence_summary.tsv
├── phylogeny/
│   ├── core_alignment.fasta
│   ├── tree.nwk
│   ├── tree.png
│   └── snp_distance_matrix.tsv
├── report/
│   ├── rMAP_Myc_Candida_report.html
│   └── assets/
└── provenance/
    ├── run_metadata.txt
    └── software_versions.tsv
```

---

## Main output files

| Output | Description |
|---|---|
| `rMAP_Myc_Candida_report.html` | Integrated final HTML report |
| `multiqc_report.html` | Aggregated sequencing QC report |
| `species_summary.tsv` | Most likely species per sample |
| `assembly_summary.tsv` | Assembly statistics |
| `busco_summary.tsv` | Genome completeness results |
| `resistance_summary.tsv` | Antifungal resistance marker summary |
| `virulence_summary.tsv` | Virulence gene summary |
| `software_versions.tsv` | Software & container versions |
| `run_metadata.txt` | Pipeline run information |

---

## Example final report sections

The integrated HTML report is designed to include the following sections:

### 1. Run summary

Includes:

- Workflow name
- Run date
- Number of samples
- Output prefix
- Pipeline version
- Analysis status

### 2. Sample overview

Includes:

- Sample ID
- FASTQ files
- Collection metadata, where available
- Analysis completion status

### 3. Read QC & trimming

Includes:

- Raw read counts
- Trimmed read counts
- Percentage reads retained
- FastQC status
- MultiQC link

### 4. Species typing

Includes:

- Most probable *Candida* species
- Kraken2/Bracken evidence
- Classification confidence
- Species-level summary plot

### 5. Genome assembly

Includes:

- Assembly length
- Number of contigs
- N50
- GC content
- Assembly interpretation

### 6. BUSCO completeness

Includes:

- Complete BUSCO percentage
- Fragmented BUSCO percentage
- Missing BUSCO percentage
- Completeness interpretation

### 7. Genome annotation

Includes:

- Number of predicted genes
- Number of coding sequences
- Protein FASTA outputs
- GFF annotation outputs

### 8. Antifungal resistance markers

Includes:

- Resistance-associated genes
- Resistance-associated mutations
- Drug classes affected
- Interpretation notes

### 9. Virulence-associated genes

Includes:

- Adhesion-related genes
- Biofilm-associated genes
- Secreted enzyme genes
- Presence/absence summaries

### 10. Phylogenomics, where enabled

Includes:

- SNP alignment
- Phylogenetic tree
- Pairwise SNP matrix
- Cluster interpretation

### 11. Provenance

Includes:

- Docker images used
- Software versions
- Databases used
- Runtime parameters

---

## Antifungal resistance interpretation

Antifungal resistance prediction in this pipeline is based on genomic markers reported in the literature or curated databases. However, genomic prediction may not fully explain antifungal susceptibility.

For clinical interpretation, genomic findings should be considered alongside:

- Phenotypic antifungal susceptibility testing
- Species identification
- Clinical history
- Treatment exposure
- Local epidemiology
- Laboratory validation

This pipeline is intended for research & surveillance use unless clinically validated in a certified diagnostic setting.

---

## Suggested databases

Depending on availability, the following database types may be used or integrated:

| Database type | Purpose |
|---|---|
| Custom *Candida* Kraken2 database | Species identification |
| NCBI RefSeq fungal genomes | Reference genome source |
| BUSCO fungal lineage datasets | Genome completeness |
| CARD or custom resistance database | Resistance gene screening |
| Custom antifungal resistance mutation database | Resistance mutation detection |
| VFDB-like or custom fungal virulence database | Virulence gene screening |

---

## Building a custom *Candida* Kraken2/Bracken database

A custom fungal or *Candida*-focused database can improve classification speed & reduce irrelevant taxonomic assignments.

Example conceptual steps:

```bash
kraken2-build --download-taxonomy --db candida_db

kraken2-build --download-library fungi --db candida_db

kraken2-build --build --db candida_db --threads 8

bracken-build -d candida_db -t 8 -k 35 -l 150
```

For a more targeted database, curated *Candida* reference genomes can be added manually.

Example:

```bash
kraken2-build --add-to-library Candida_albicans.fna --db candida_db
kraken2-build --add-to-library Candida_auris.fna --db candida_db
kraken2-build --add-to-library Candida_glabrata.fna --db candida_db
kraken2-build --build --db candida_db --threads 8
```

---

## Docker images

The pipeline is designed to use Docker containers for reproducibility.

Example image categories:

| Module | Example Docker image |
|---|---|
| Trimming | `staphb/trimmomatic:0.39` |
| FastQC | `staphb/fastqc:0.11.9` |
| MultiQC | `multiqc/multiqc` |
| Kraken2/Bracken | Custom *Candida* Kraken2/Bracken image |
| Assembly | SPAdes or Shovill image |
| QUAST | QUAST image |
| BUSCO | BUSCO image |
| Annotation | Funannotate or Prokka image |
| AMR screening | Custom antifungal resistance image |
| Reporting | Python/R reporting image |

---

## Example local execution on macOS

If using macOS with Colima:

```bash
colima start --cpu 8 --memory 16 --disk 150
```

Confirm Docker:

```bash
docker ps
```

Run workflow:

```bash
java -jar ~/cromwell-92.jar run rMAP_Candida.wdl \
  --inputs ~/rMAP-Candida.inputs.example.json
```

---

## Troubleshooting

### Docker is not running

Error example:

```text
Cannot connect to the Docker daemon
```

Fix:

```bash
docker ps
```

If using Colima:

```bash
colima start --cpu 8 --memory 16 --disk 150
```

---

### Cromwell cannot localize files

Check that all file paths in the input JSON are absolute paths.

Good:

```json
"~/Sample_001_R1.fastq.gz"
```

---

### BUSCO fails

Possible causes:

- Incorrect lineage dataset
- Poor assembly quality
- Missing BUSCO database
- Insufficient memory
- Docker image issue

Check:

```bash
find cromwell-executions -path "*/call-BUSCO/*/execution/*" -type f
```

Review BUSCO logs:

```bash
cat busco*.log
```

---

### Kraken2 database errors

Possible causes:

- Database path is incorrect
- Database was not fully built
- Database files are missing
- Memory is insufficient

Check database directory:

```bash
ls -lh ~/candida_kraken2_db
```

Expected files may include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

---

### AMR task fails

Possible causes:

- Resistance database path is missing
- Input assembly file is empty
- Gene database is not indexed
- Tool expects nucleotide input but received protein input, or vice versa
- Output filename expected by WDL does not match the actual tool output

Check AMR execution directory:

```bash
find cromwell-executions -path "*/call-AMR*/*/execution/*" -type f
```

Review logs:

```bash
cat stderr
cat stdout
```

---

## Recommended quality thresholds

Suggested initial thresholds:

| Metric | Suggested threshold |
|---|---|
| Raw reads | >500,000 paired reads |
| Trimmed reads retained | >70% |
| Assembly size | Species-dependent |
| Contigs | Lower is better |
| N50 | Higher is better |
| BUSCO completeness | >90% preferred |
| Species typing | Dominant species supported by high read proportion |

These thresholds should be adapted based on sequencing depth, species, study design & laboratory context.

---

## Limitations

- Genomic prediction of antifungal resistance may not fully match phenotypic susceptibility.
- Resistance interpretation depends on database completeness & species-specific knowledge.
- Species identification depends on the completeness & quality of the reference database.
- Mixed infections or contamination may complicate interpretation.
- Fungal genomes can vary substantially in size, ploidy, repeat content & assembly complexity.
- The pipeline requires validation before clinical diagnostic use.

---

## Recommended citation

If you use this pipeline, please cite the repository & any tools used by the workflow.

Suggested citation format:

```text
rMAP-Candida: Rapid Mycology Analysis Pipeline for Candida genomic surveillance. 
```
---

## Contributing

Contributions are welcome.

Potential areas for contribution include:

- Improved antifungal resistance databases
- Species-specific resistance interpretation rules
- Additional fungal pathogens
- Improved HTML reporting
- Phylogenomic modules
- Benchmarking datasets
- Cloud execution support
- Documentation & tutorials

To contribute:

1. Fork the repository
2. Create a new branch
3. Make your changes
4. Test the workflow
5. Submit a pull request

---

## Development roadmap

Planned improvements include:

- Expanded support for non-*Candida* fungal pathogens
- Species-specific resistance interpretation
- Integration of curated *Candida auris* resistance markers
- Improved outbreak investigation module
- Pairwise SNP clustering
- Interactive HTML visualizations
- GitHub Pages example reports
- Terra/Google Cloud execution support
- Benchmarking with public clinical datasets

---
## License

This project is released under the terms specified in the `LICENSE` file.

Recommended license options:

- MIT License

---

## Disclaimer

This pipeline is intended for research, training & public health genomic surveillance. It is not intended for direct clinical decision-making unless appropriately validated and implemented under relevant regulatory and laboratory quality frameworks.
