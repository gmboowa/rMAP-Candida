## rMAP-Candida: Rapid Mycological Analysis Pipeline for *Candida spp* Genomic Surveillance

<p align="center">
  <img src="/docs/assets/Project.jpg" alt="rMAP-Candida workflow illustration" width="100%">
</p>

<p align="center">
  <em>rMAP-Candida workflow overview for Candida spp. whole-genome sequencing, species typing, assembly assessment, antifungal-resistance marker screening, phylogenomics & surveillance reporting.</em>
</p>

---


**rMAP-Candida** is a Dockerized WDL/Cromwell workflow for paired-end *Candida* whole-genome sequencing analysis. It integrates read QC and trimming, *Candida*-focused species typing, de novo assembly, assembly-contiguity assessment, optional genome-completeness assessment, curated antifungal-resistance marker screening, species-aware phylogenomics, pairwise SNP-distance summaries, and an integrated HTML surveillance report.

The workflow is intended for **research, training, method benchmarking, outbreak-investigation support, and public-health genomic surveillance**. It is **not a standalone clinical diagnostic system**. Antifungal-resistance results are reported as genomic screening evidence and should be interpreted with species identity, sample metadata, validated catalogues, and phenotypic antifungal susceptibility testing where required.


---

## Quick start: local Cromwell test


### 1. Clone the repository

```bash
git clone https://github.com/gmboowa/rMAP-Candida.git
cd rMAP-Candida
```

### 2. Prepare the two-sample test FASTQs

The example JSON expects these files:

```text
example/fastq/ERR263534_1.fastq.gz
example/fastq/ERR263534_2.fastq.gz
example/fastq/ERR331060_1.fastq.gz
example/fastq/ERR331060_2.fastq.gz
```
Prepare the two-sample test FASTQs

The example JSON expects these files:

```text
example/fastq/ERR263534_1.fastq.gz
example/fastq/ERR263534_2.fastq.gz
example/fastq/ERR331060_1.fastq.gz
example/fastq/ERR331060_2.fastq.gz
```

Create a FASTQ directory:

```bash
mkdir -p example/fastq
cd example/fastq
fastq-dump --split-3 --gzip ERR263534
fastq-dump --split-3 --gzip ERR331060
cd ../..
```

Download the two test samples from SRA using `fastq-dump`:

```bash
fastq-dump --split-3 --gzip ERR263534
fastq-dump --split-3 --gzip ERR331060
```

Confirm that the expected paired-end FASTQ files are downloaded:


If `fastq-dump` is not installed, install the SRA Toolkit first.

On macOS with Homebrew:

```bash
brew install sratoolkit
```

Then confirm that `fastq-dump` is available:

```bash
fastq-dump --version
```

Install Docker and confirm it is running

rMAP-Candida uses Docker containers for the workflow tools and databases. Docker must be installed and running before launching Cromwell.

#### macOS

Install Docker Desktop from:

```text
https://www.docker.com/products/docker-desktop/
```

After installation, open **Docker Desktop** from Applications and wait until it says Docker is running.

Confirm Docker works:

```bash
docker --version
docker info
```

You can also test Docker with:

```bash
docker run hello-world
```

#### Linux

Install Docker using your system package manager or the official Docker installation guide:

```text
https://docs.docker.com/engine/install/
```

After installation, start Docker:

```bash
sudo systemctl start docker
sudo systemctl enable docker
```

Confirm Docker works:

```bash
docker --version
docker info
docker run hello-world
```

If Docker requires `sudo`, either run Docker commands with `sudo` or add your user to the Docker group according to your institution’s system policy.

#### Check rMAP-Candida Docker images

Before running Cromwell, confirm that key workflow containers can be pulled:

```bash
docker pull gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
docker pull gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9
docker pull gmboowa/rmap-candida-refs:2026.05
```

Check that the Candida Kraken2/Bracken database is visible inside the container:

```bash
docker run --rm gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db \
  bash -lc 'ls -lh /opt/kraken2_db/candida && ls /opt/kraken2_db/candida | head'
```

Expected database files include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

### 3. Run the workflow without a custom Cromwell config

From the repository directory:

```bash
cd ~/rMAP-Candida
```

Run the two-sample test workflow:

```bash
java -jar ~/cromwell-91.jar run \
  rMAP_Candida.wdl \
  --inputs example/rMAP-Candida.inputs.example.two_samples.json
```

If your local Cromwell/Docker setup produces Docker hash-lookup or Docker metadata errors, use the repository-provided no-Docker-hash-lookup Cromwell config instead.

### 4. Confirm success

```bash
grep -i "workflow.*succeeded" cromwell-*.log || true
find cromwell-executions/rMAP_Candida -name "rMAP_Candida_report.html" -print
find cromwell-executions/rMAP_Candida \( \
  -name "rMAP_Candida_summary.tsv" -o \
  -name "rMAP_Candida_surveillance_summary.tsv" -o \
  -name "rMAP_Candida_pairwise_snp_distances.tsv" \
\) -print
```

The successful test report should summarize samples, top species groups, one curated AMR marker hit, and a median assembly N50 in bp. Species-aware phylogeny may be skipped for the two-sample test dataset if there are fewer than three same-species samples.

### 3. Run the workflow without a custom Cromwell config


```bash
java -jar ~/cromwell-91.jar run rMAP_Candida.wdl --inputs ~/rMAP-Candida.inputs.json
```

### 4. Confirm success

```bash
grep -i "workflow.*succeeded" cromwell-*.log || true
find cromwell-executions/rMAP_Candida -name "rMAP_Candida_report.html" -print
find cromwell-executions/rMAP_Candida \( \
  -name "rMAP_Candida_summary.tsv" -o \
  -name "rMAP_Candida_surveillance_summary.tsv" -o \
  -name "rMAP_Candida_pairwise_snp_distances.tsv" \
\) -print
```

The successful test report should summarize samples, top species groups, one curated AMR marker hit, and a median assembly N50 in bp. Species-aware phylogeny.

---



## Repository layout

```text
rMAP-Candida/
├── README.md
├── rMAP_Candida.wdl
├── docker/
├
├── example/
│   ├── rMAP_Candida.inputs.example.json
│   ├── rMAP_Candida_accessions.tsv
│   └── surveillance_metadata.tsv
├
└── docs/
    ├── assets/
    ├── reports/
    └── index.html

```

---

## Main workflow modules

| Module | WDL file | Main tools | Purpose | Key outputs |
|---|---|---|---|---|
| Read QC and trimming | `wdl/modules/read_qc.wdl` | fastp, FastQC, MultiQC | Adapter/quality trimming & read-level QC | `*.fastp.html`, `*.fastp.json`, MultiQC HTML |
| Species typing | `wdl/modules/species_typing.wdl` | Kraken2, Bracken | *Candida*-focused species classification | `*.kraken2.report`, `*.bracken.tsv`, `*.top_species.tsv` |
| Assembly and assembly QC | `wdl/modules/assembly_qc.wdl` | MEGAHIT, QUAST, Compleasm/BUSCO | De novo assembly, contiguity metrics, optional completeness | contigs FASTA, QUAST TSV, Compleasm/BUSCO summaries |
| AMR screening | `wdl/modules/amr.wdl` | ChroQueTas/FungAMR-derived parsing | Curated genomic antifungal-resistance marker screening | AMR summary/raw TSV and HTML |
| Species-aware phylogeny | `wdl/modules/phylogeny.wdl` | Snippy, BCFtools, IQ-TREE2, ETE3 | Per-species core-SNP analysis & tree rendering | core SNP FASTA, Newick tree, PNG tree, group summary |
| Reporting | `wdl/modules/reporting.wdl` | built-in report logic | Integrated surveillance report & TSV summaries | `rMAP_Candida_report.html`, summary TSVs |

---

## Required inputs

The workflow expects paired-end FASTQ files & matching sample names.

| Input | Type | Description |
|---|---|---|
| `rMAP_Candida.sample_names` | `Array[String]+` | Sample IDs. Must be in the same order as R1 & R2 files. |
| `rMAP_Candida.read1s` | `Array[File]+` | Read 1 FASTQ files. |
| `rMAP_Candida.read2s` | `Array[File]+` | Read 2 FASTQ files. |

Minimal two-sample example:

```json
{
  "rMAP_Candida.sample_names": ["ERR263534", "ERR331060"],
  "rMAP_Candida.read1s": [
    "example/fastq/ERR263534_1.fastq.gz",
    "example/fastq/ERR331060_1.fastq.gz"
  ],
  "rMAP_Candida.read2s": [
    "example/fastq/ERR263534_2.fastq.gz",
    "example/fastq/ERR331060_2.fastq.gz"
  ],
  "rMAP_Candida.do_trimming": true,
  "rMAP_Candida.do_quality_control": true,
  "rMAP_Candida.do_species_typing": true,
  "rMAP_Candida.do_assembly": true,
  "rMAP_Candida.do_assembly_qc": true,
  "rMAP_Candida.do_compleasm": false,
  "rMAP_Candida.do_busco": false,
  "rMAP_Candida.do_fungal_amr": false,
  "rMAP_Candida.do_phylogeny": false,
  "rMAP_Candida.render_phylogeny_tree": false,
  "rMAP_Candida.max_cpus": 2,
  "rMAP_Candida.max_memory_gb": 8,
  "rMAP_Candida.max_disk_gb": 100
}
```

---

## Optional metadata

The workflow can include optional surveillance metadata in the HTML report.

```json
"rMAP_Candida.surveillance_metadata_tsv": "example/surveillance_metadata.example.tsv"
```

Recommended columns:

```text
sample_id	country	site	collection_date	specimen_type	patient_group	ward_or_facility	sequencing_platform	species	metadata_note
```

If no metadata file is supplied, the workflow still runs & the report states that no surveillance metadata was provided.

---

## Containers, databases, and public access

All major steps are Dockerized. Users should pull and inspect the containers before running the workflow.

### Species typing database

The custom *Candida*-focused Kraken2/Bracken database is bundled in:

```text
gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
```

Inside the container, the database path is:

```text
/opt/kraken2_db/candida
```

Check public access:

```bash
docker pull gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
docker run --rm gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db \
  bash -lc 'ls -lh /opt/kraken2_db/candida && ls /opt/kraken2_db/candida | head'
```

Expected database files include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

The database build recipe and source-list notes should be maintained in:

```text
docker/kraken2_bracken/README.md
```

### AMR container

```bash
docker pull gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9
```

### Reference-genome container

```bash
docker pull gmboowa/rmap-candida-refs:2026.05
docker run --rm gmboowa/rmap-candida-refs:2026.05
bash -lc 'cat /opt/rmap_candida_refs/references.tsv'
```

---

## Expected outputs

Cromwell stores outputs under:

```text
cromwell-executions/rMAP_Candida/<workflow-id>/
```

Main report outputs:

```text
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_report.html
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_summary.tsv
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_surveillance_summary.tsv
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_pairwise_snp_distances.tsv
```

Additional outputs may include:

```text
*.fastp.html
*.fastp.json
*.kraken2.report
*.bracken.tsv
*.top_species.tsv
*.contigs.fasta
*.assembly_summary.tsv
*.quast.report.tsv
*.compleasm.summary.tsv
*.fungal_amr.summary.tsv
*.fungal_amr.raw.tsv
*.fungal_amr.html
*.core_snp_alignment.fasta
*.treefile
*.nwk
*.png
```

---

## Species-aware phylogenomics & diploidy/polymorphism handling

rMAP-Candida does **not** combine different *Candida* species into one core-SNP tree. It first groups isolates by species and only builds a tree when enough same-species samples and a configured reference are available.

This design is deliberate because *Candida* and related yeasts differ in ploidy, genome plasticity, heterozygosity, copy-number variation, aneuploidy, loss of heterozygosity, and recombination behavior. The workflow therefore reports phylogeny eligibility explicitly and treats SNP distances as descriptive genomic relatedness metrics, not direct transmission proof.

Default phylogeny settings:

```json
"rMAP_Candida.min_species_samples_for_tree": 3,
"rMAP_Candida.snippy_phylogeny_species": ["Candidozyma auris", "Candida albicans"],
"rMAP_Candida.haploid_phylogeny_species": ["Nakaseomyces glabratus"],
"rMAP_Candida.core_site_min_fraction": 0.95
```

If a species has fewer than three eligible samples, the report will state that phylogeny was skipped because there were too few same-species samples. This is expected behavior for the two-sample test dataset.

---

## Antifungal-resistance interpretation

AMR results are genomic screening outputs. They do **not** automatically indicate clinical susceptibility or resistance.

The report distinguishes:

```text
curated genomic AMR marker detected
```

from:

```text
no curated genomic AMR marker detected - susceptibility not inferred
```

Resistance may involve mechanisms that are not fully captured by short-read marker screening or the configured catalogue, including:

- point mutations in known loci such as `ERG11`, `FKS1`, `FKS2`, or species-specific orthologues;
- promoter or regulatory changes such as `TAC1`, `UPC2`, `MRR1`, or `PDR1`-associated mechanisms;
- copy-number changes;
- aneuploidy;
- loss of heterozygosity;
- structural variants;
- species-specific mechanisms not included in the installed database.

Clinically important findings should be interpreted with isolate metadata & phenotypic antifungal susceptibility testing where available.

---

## Why WDL/Cromwell?

WDL/Cromwell was selected because it provides a portable workflow description, explicit typed inputs and outputs, scatter-based parallelism for sample-level tasks, Docker runtime declarations, and execution provenance. The same workflow can be run locally for a two-sample reviewer test, on a workstation or HPC environment for larger batches, or on WDL-compatible cloud platforms.

This implementation is designed to reduce installation barriers: users do not need to install each bioinformatics tool manually, because task-level dependencies are containerized and versioned.

---

## Troubleshooting

### The workflow works without a config file

Use the no-config command first:

```bash
java -jar ~/cromwell-91.jar run rMAP_Candida.wdl --inputs example/rMAP-Candida.inputs.json
```

Only use a custom Cromwell config if your local Docker setup produces hash-lookup or Docker metadata errors.

### Docker image cannot be pulled

Check the exact tag:

```bash
docker pull gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
docker pull gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9
docker pull gmboowa/rmap-candida-refs:2026.05
```

### FASTQ file paths fail

Check the files in the JSON exist from the directory where you run Cromwell:

```bash
jq -r '."rMAP_Candida.read1s"[], ."rMAP_Candida.read2s"[]' rMAP-Candida.inputs.json | xargs -I{} ls -lh {}
```

### Phylogeny is skipped

This is expected when there are fewer than three same-species samples. The workflow reports this rather than silently failing.

---

## Interpretation & limitations

- rMAP-Candida is not a standalone clinical diagnostic system.
- Species calls depend on the quality & taxonomic coverage of the Kraken2/Bracken database.
- Mixed cultures, contamination, low sequencing depth, or incomplete databases can reduce species-confidence scores.
- QUAST reports assembly contiguity, not gene-level completeness.
- Compleasm/BUSCO completeness depends on assembly quality & lineage selection.
- AMR marker absence does not imply susceptibility.
- Copy-number variation, aneuploidy, loss of heterozygosity, promoter changes, regulatory mechanisms & structural variation may not be fully captured.
- Species-aware trees & low SNP distances should not be interpreted as direct transmission without epidemiological metadata & species-appropriate validation.
- Public SRA/ENA/BioSample metadata may be incomplete & should be verified before surveillance interpretation.

---

## Suggested citation

If you use this workflow, cite the repository and the major tools used in the workflow, including WDL/Cromwell, Docker, fastp, FastQC, MultiQC, Kraken2, Bracken, MEGAHIT, QUAST, Compleasm/BUSCO where enabled, ChroQueTas/FungAMR-derived resources, Snippy, IQ-TREE2, and ETE3.

```text
rMAP-Candida: Rapid Mycological Analysis Pipeline for Candida genomic surveillance.
https://github.com/gmboowa/rMAP-Candida
```
