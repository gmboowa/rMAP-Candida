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

### 3. Run the complete workflow locally from raw paired-end reads

Local Docker/Colima execution of rMAP-Candida has been validated with **Cromwell 91** and the supplied local backend configuration:

```text
example/cromwell.local.fifo_portable.conf
```

For the validated local setup, do not run Cromwell without the configuration file. The configuration:

- selects Cromwell's Local/Shared File System backend;
- limits concurrent tasks to reduce CPU and memory oversubscription on a laptop;
- converts WDL `cpu` and `memory` runtime attributes into Docker `--cpus` and `--memory` limits;
- passes the container-visible `${docker_script}` path to Docker;
- places Cromwell and tool temporary FIFOs under native `/tmp`, avoiding named-pipe failures on macOS/Colima bind-mounted directories; and
- disables Docker hash lookup for this validated local environment.

> **Local Cromwell compatibility:** Use `cromwell-91.jar` with this configuration. During local validation, Cromwell 92 produced a Local/SFS container-runtime coercion error before a downstream Docker task could launch. This was an execution-engine compatibility issue, not a failure of the rMAP-Candida task logic.

#### 3.1 Prerequisites

Confirm that Java and Docker are available:

```bash
java -version
docker info
```

When using Colima on macOS:

```bash
colima status
```

Place Cromwell 91 somewhere accessible, for example:

```text
~/cromwell-91.jar
```

Confirm the version:

```bash
java -jar ~/cromwell-91.jar --version
```

Expected:

```text
cromwell 91
```

#### 3.2 Prepare the raw-read input JSON

Copy the supplied example:

```bash
cp example/rMAP_Candida.inputs.example.json \
   example/rMAP_Candida.local_test.inputs.json
```

Edit `example/rMAP_Candida.local_test.inputs.json` so that:

- `sample_names`, `read1s`, and `read2s` have the same number of entries;
- each R1 file is paired with the corresponding R2 file;
- FASTQ paths are absolute local paths or otherwise visible to Cromwell;
- the full raw-read workflow modules are enabled; and
- preassembled-contig resume mode is disabled.

A minimal two-sample raw-read configuration has the following form:

```json
{
  "rMAP_Candida.sample_names": [
    "sample_01",
    "sample_02"
  ],
  "rMAP_Candida.read1s": [
    "/absolute/path/to/sample_01_R1.fastq.gz",
    "/absolute/path/to/sample_02_R1.fastq.gz"
  ],
  "rMAP_Candida.read2s": [
    "/absolute/path/to/sample_01_R2.fastq.gz",
    "/absolute/path/to/sample_02_R2.fastq.gz"
  ],

  "rMAP_Candida.do_trimming": true,
  "rMAP_Candida.do_quality_control": true,
  "rMAP_Candida.do_species_typing": true,
  "rMAP_Candida.do_assembly": true,
  "rMAP_Candida.do_assembly_qc": true,
  "rMAP_Candida.do_compleasm": true,
  "rMAP_Candida.do_busco": false,
  "rMAP_Candida.do_fungal_amr": true,
  "rMAP_Candida.do_phylogeny": true,
  "rMAP_Candida.render_phylogeny_tree": true,

  "rMAP_Candida.use_preassembled_contigs": false,

  "rMAP_Candida.fungal_kraken2_bracken_docker": "gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db",
  "rMAP_Candida.kraken_db_path": "/opt/kraken2_db/candida",
  "rMAP_Candida.fungamr_docker": "gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9",
  "rMAP_Candida.compleasm_docker": "huangnengcsu/compleasm:v0.2.7",
  "rMAP_Candida.candida_refs_docker": "gmboowa/rmap-candida-refs:2026.05",
  "rMAP_Candida.candida_refs_manifest": "/opt/rmap_candida_refs/references.tsv",

  "rMAP_Candida.megahit_disable_hw_acceleration": true,
  "rMAP_Candida.megahit_memory_percent": 65,
  "rMAP_Candida.megahit_mem_flag": 0,
  "rMAP_Candida.fastqc_memory_mb": 1024,

  "rMAP_Candida.min_species_samples_for_tree": 3,
  "rMAP_Candida.iqtree2_bootstraps": 1000,

  "rMAP_Candida.max_cpus": 2,
  "rMAP_Candida.max_memory_gb": 8
}
```

For a species-specific phylogeny, at least `min_species_samples_for_tree` samples must have the same species call and a matching reference in the bundled Candida reference manifest. Mixed-species trees are intentionally not generated.

#### 3.3 Run the raw-read workflow in the foreground

From the repository root:

```bash
CROMWELL_JAR="$HOME/cromwell-91.jar"
CROMWELL_CONF="$PWD/example/cromwell.local.fifo_portable.conf"
WDL="$PWD/rMAP_Candida.wdl"
INPUTS="$PWD/example/rMAP_Candida.local_test.inputs.json"

java \
  -Xms1g \
  -Xmx3g \
  -XX:+UseG1GC \
  -Dconfig.file="$CROMWELL_CONF" \
  -jar "$CROMWELL_JAR" \
  run "$WDL" \
  --inputs "$INPUTS"
```

This command starts the complete workflow from the raw paired-end FASTQ files. When enabled in the JSON, the workflow performs:

1. read trimming with fastp;
2. raw/trimmed-read quality assessment;
3. Candida-focused Kraken2/Bracken species typing;
4. de novo MEGAHIT assembly;
5. QUAST assembly-contiguity assessment;
6. optional Compleasm/BUSCO completeness analysis;
7. fungal AMR characterization;
8. species-aware core-SNP phylogeny;
9. IQ-TREE inference and tree rendering; and
10. generation of the final integrated HTML and TSV reports.

#### 3.4 Run the raw-read workflow in the background on macOS

For a long local run, use `nohup` to detach the process from the terminal and `caffeinate` to prevent idle sleep while macOS remains awake:

```bash
cd /path/to/rMAP-Candida

CROMWELL_JAR="$HOME/cromwell-91.jar"
CROMWELL_CONF="$PWD/example/cromwell.local.fifo_portable.conf"
WDL="$PWD/rMAP_Candida.wdl"
INPUTS="$PWD/example/rMAP_Candida.local_test.inputs.json"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG="$PWD/rMAP_Candida_raw_reads_${RUN_ID}.log"
PIDFILE="$PWD/rMAP_Candida_raw_reads_${RUN_ID}.pid"

nohup caffeinate -dimsu \
  java \
  -Xms1g \
  -Xmx3g \
  -XX:+UseG1GC \
  -Dconfig.file="$CROMWELL_CONF" \
  -jar "$CROMWELL_JAR" \
  run "$WDL" \
  --inputs "$INPUTS" \
  > "$LOG" 2>&1 < /dev/null &

PID=$!
echo "$PID" > "$PIDFILE"
disown "$PID"

echo "rMAP-Candida raw-read workflow started"
echo "PID: $PID"
echo "PID file: $PIDFILE"
echo "Log file: $LOG"
```

Closing the terminal does not stop the detached process. Closing a MacBook lid normally suspends macOS and pauses Docker/Colima, so keep the machine awake and connected to power for an uninterrupted run.

Monitor the log:

```bash
tail -f "$LOG"
```

Press `Ctrl-C` to stop viewing the log; this does not stop Cromwell.

Check whether Cromwell is still running:

```bash
PID="$(cat "$PIDFILE")"

if kill -0 "$PID" 2>/dev/null; then
  echo "Cromwell is running: PID $PID"
  ps -p "$PID" -o pid,ppid,state,etime,%cpu,%mem,command
else
  echo "Cromwell is no longer running"
fi
```

View the active Docker task:

```bash
docker ps --format \
'table {{.ID}}\t{{.Status}}\t{{.Image}}'
```

View current container resource use:

```bash
docker stats --no-stream
```

### 4. Confirm successful completion and final report generation

#### 4.1 Confirm the workflow succeeded

For a background run:

```bash
grep -E \
"Workflow actor.*completed with status 'Succeeded'|workflow finished with status 'Succeeded'|transitioned to state Succeeded" \
"$LOG"
```

A successful run should include a message similar to:

```text
SingleWorkflowRunnerActor workflow finished with status 'Succeeded'
```

Also check for explicit failures:

```bash
grep -E \
'exited with return code|completed with status .Failed.|workflow finished with status .Failed.|transitioned to state Failed' \
"$LOG" || true
```

#### 4.2 Locate the final integrated HTML report

```bash
REPORT="$(find cromwell-executions/rMAP_Candida \
  -path '*/call-MERGE_MYC_REPORTS/execution/rMAP_Candida_report.html' \
  -type f -print0 |
  xargs -0 ls -t 2>/dev/null |
  head -n 1)"

if [ -n "$REPORT" ] && [ -s "$REPORT" ]; then
  echo "Final HTML report: $REPORT"
  ls -lh "$REPORT"
else
  echo "Final HTML report was not found or is empty"
  exit 1
fi
```

Open it on macOS:

```bash
open "$REPORT"
```

#### 4.3 Confirm the AMR and phylogeny sections are present

```bash
grep -q 'id="amr"' "$REPORT" &&
  echo "AMR section present"

grep -q 'id="phylogeny"' "$REPORT" &&
  echo "Phylogeny section present"
```

For a phylogeny-eligible dataset, confirm that the report contains embedded tree content:

```bash
grep -Eq \
'class="tree-img"|data:image/png;base64|<svg' \
"$REPORT" &&
  echo "Rendered phylogeny content present"
```

If no species group meets the sample-count/reference requirements, the phylogeny section remains in the final report and explains why no tree was generated.

#### 4.4 Confirm that sample-level outputs were produced

Count successful assembly tasks:

```bash
find cromwell-executions/rMAP_Candida \
  -path '*/call-ASSEMBLY/shard-*/execution/rc' \
  -exec sh -c '
    for rc do
      [ "$(cat "$rc")" = "0" ] && echo "$rc"
    done
  ' sh {} + |
wc -l
```

Locate final contigs:

```bash
find cromwell-executions/rMAP_Candida \
  -path '*/call-ASSEMBLY/shard-*/execution/*.contigs.fasta' \
  -type f -size +0c -print
```

Locate AMR summaries:

```bash
find cromwell-executions/rMAP_Candida \
  -path '*/call-AMR/shard-*/execution/*.fungal_amr.summary.tsv' \
  -type f -size +0c -print
```

Locate species-aware tree outputs:

```bash
find cromwell-executions/rMAP_Candida \
  -type f \
  \( -name '*.treefile' -o -name '*.nwk' -o -name '*.png' \) \
  -size +0c -print
```

Locate the main tabular outputs:

```bash
find cromwell-executions/rMAP_Candida \( \
  -name "rMAP_Candida_summary.tsv" -o \
  -name "rMAP_Candida_surveillance_summary.tsv" -o \
  -name "rMAP_Candida_pairwise_snp_distances.tsv" \
\) -print
```

A complete raw-read validation should demonstrate that rMAP-Candida starts with paired FASTQ files, completes assembly and all enabled downstream analyses, and produces the non-empty integrated HTML report. AMR findings are genomic screening evidence and should be interpreted with species identity, validated marker catalogues, clinical metadata, and phenotypic antifungal susceptibility testing where appropriate.
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
