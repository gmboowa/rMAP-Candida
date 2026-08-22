## rMAP-Candida: Rapid Mycological Analysis Pipeline for *Candida* spp. Genomic Surveillance

<p align="center">
  <img src="/docs/assets/Project.jpg" alt="rMAP-Candida workflow illustration" width="100%">
</p>

<p align="center">
  <em>rMAP-Candida workflow overview for *Candida* spp. whole-genome sequencing, species typing, assembly assessment, antifungal-resistance marker screening, phylogenomics & surveillance reporting.</em>
</p>

---


**rMAP-Candida** is a Dockerized WDL/Cromwell workflow for paired-end *Candida* whole-genome sequencing analysis. It integrates read QC & trimming, *Candida*-focused species typing, *de novo* assembly, assembly-contiguity assessment, optional genome-completeness assessment, curated antifungal-resistance marker screening, species-aware phylogenomics, pairwise SNP-distance summaries & an integrated HTML surveillance report.

The workflow is intended for **research, training, method benchmarking, outbreak-investigation support & public-health genomic surveillance**. It is **not a standalone clinical diagnostic system**. Antifungal-resistance results are reported as genomic screening evidence & should be interpreted with species identity, sample metadata, validated catalogues & phenotypic antifungal susceptibility testing where required.


---

## Quick start: local Cromwell test


### 1. Clone the repository

```bash
git clone https://github.com/gmboowa/rMAP-Candida.git
cd rMAP-Candida
```

### 2. Prepare the two-sample raw-read test data

The example starts from paired-end raw FASTQ files & runs the complete workflow, including trimming, quality control, species typing, *de novo* assembly, assembly assessment, antifungal-resistance screening, phylogenomics where eligible & final HTML reporting.

The example input JSON expects:

```text
example/fastq/ERR263534_1.fastq.gz
example/fastq/ERR263534_2.fastq.gz
example/fastq/ERR331060_1.fastq.gz
example/fastq/ERR331060_2.fastq.gz
```

Create the FASTQ directory & download the two public test runs:

```bash
mkdir -p example/fastq
cd example/fastq

fastq-dump --split-3 --gzip ERR263534
fastq-dump --split-3 --gzip ERR331060

cd ../..
```

Confirm that all four files are present & non-empty:

```bash
ls -lh \
  example/fastq/ERR263534_1.fastq.gz \
  example/fastq/ERR263534_2.fastq.gz \
  example/fastq/ERR331060_1.fastq.gz \
  example/fastq/ERR331060_2.fastq.gz
```

If `fastq-dump` is not installed, install the SRA Toolkit first. On macOS with Homebrew:

```bash
brew install sratoolkit
fastq-dump --version
```

### 3. Confirm Java, Docker & Cromwell 91

rMAP-Candida uses Dockerized tasks. Confirm that Java & Docker are available:

```bash
java -version
docker --version
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

Expected output should identify Cromwell version 91.

Before launching the workflow, confirm that the main custom images are publicly accessible:

```bash
docker pull gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
docker pull gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9
docker pull gmboowa/rmap-candida-refs:2026.05
```

Confirm that the bundled Candida Kraken2/Bracken database is available inside its container:

```bash
docker run --rm \
  gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db \
  bash -lc 'ls -lh /opt/kraken2_db/candida && ls /opt/kraken2_db/candida | head'
```

Expected database files include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

### 4. Use the supplied local Cromwell configuration

Validated local Docker/Colima execution uses:

```text
example/cromwell.local.fifo_portable.conf
```

The configuration file:

- selects Cromwell's Local/Shared File System backend;
- limits concurrent tasks to reduce CPU & memory oversubscription on a local machine;
- maps WDL `cpu` & `memory` runtime attributes to Docker `--cpus` & `--memory` limits;
- launches the container-visible `${docker_script}` path;
- places Cromwell temporary files & FIFOs under native `/tmp`, avoiding named-pipe failures on macOS/Colima bind-mounted directories; &
- disables Docker hash lookup for the validated local execution environment.

> **Local compatibility note:** Use `cromwell-91.jar` together with this configuration for the validated local Docker/Colima route. Cromwell 92 produced a Local/SFS container-runtime coercion error during local testing before a downstream Docker task could launch.

### 5. Confirm that the input JSON runs from raw reads

The supplied input should be:

```text
example/rMAP-Candida.inputs.example.two_samples.json
```

Confirm that it enables *de novo* assembly & does not use preassembled contigs:

```bash
jq '{
  do_trimming: ."rMAP_Candida.do_trimming",
  do_quality_control: ."rMAP_Candida.do_quality_control",
  do_species_typing: ."rMAP_Candida.do_species_typing",
  do_assembly: ."rMAP_Candida.do_assembly",
  do_assembly_qc: ."rMAP_Candida.do_assembly_qc",
  do_fungal_amr: ."rMAP_Candida.do_fungal_amr",
  do_phylogeny: ."rMAP_Candida.do_phylogeny",
  use_preassembled_contigs: ."rMAP_Candida.use_preassembled_contigs"
}' example/rMAP-Candida.inputs.example.two_samples.json
```

For a complete raw-read run, the relevant values should include:

```json
{
  "do_trimming": true,
  "do_quality_control": true,
  "do_species_typing": true,
  "do_assembly": true,
  "do_assembly_qc": true,
  "do_fungal_amr": true,
  "do_phylogeny": true,
  "use_preassembled_contigs": false
}
```

The arrays `sample_names`, `read1s`, & `read2s` must have the same number of entries & must use the same sample order.

### 6. Run the complete raw-read workflow

#### Foreground execution

From the repository root:

```bash
CROMWELL_JAR="$HOME/cromwell-91.jar"
CROMWELL_CONF="$PWD/example/cromwell.local.fifo_portable.conf"
WDL="$PWD/rMAP_Candida.wdl"
INPUTS="$PWD/example/rMAP-Candida.inputs.example.two_samples.json"

java \
  -Xms1g \
  -Xmx3g \
  -XX:+UseG1GC \
  -Dconfig.file="$CROMWELL_CONF" \
  -jar "$CROMWELL_JAR" \
  run "$WDL" \
  --inputs "$INPUTS"
```

#### Background execution on macOS

For a long local run, use `nohup` to detach the process from the terminal & `caffeinate` to prevent idle sleep while laptop remains awake:

```bash
cd ~/rMAP-Candida

CROMWELL_JAR="$HOME/cromwell-91.jar"
CROMWELL_CONF="$PWD/example/cromwell.local.fifo_portable.conf"
WDL="$PWD/rMAP_Candida.wdl"
INPUTS="$PWD/example/rMAP-Candida.inputs.example.two_samples.json"

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

Closing the terminal does not stop the detached process. Closing the laptop lid normally suspends the run & pauses Docker/Colima, so keep the machine awake & connected to power for an uninterrupted run.

Monitor the run:

```bash
tail -f "$LOG"
```

Press `Ctrl-C` to stop viewing the log; this does not stop Cromwell.

Check whether the process is still running:

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
docker ps --format 'table {{.ID}}\t{{.Status}}\t{{.Image}}'
```

### 7. Confirm successful completion & final report generation

For a background run, confirm that Cromwell completed successfully:

```bash
grep -E \
"Workflow actor.*completed with status 'Succeeded'|workflow finished with status 'Succeeded'|transitioned to state Succeeded" \
"$LOG"
```

Check for explicit failures:

```bash
grep -E \
'exited with return code|completed with status .Failed.|workflow finished with status .Failed.|transitioned to state Failed' \
"$LOG" || true
```

Locate the most recently generated integrated report:

```bash
REPORT="$(find cromwell-executions/rMAP_Candida \
  -path '*/call-MERGE_MYC_REPORTS/execution/rMAP_Candida_report.html' \
  -type f -print0 | \
  xargs -0 ls -t 2>/dev/null | \
  head -n 1)"

if [ -n "$REPORT" ] && [ -s "$REPORT" ]; then
  echo "Final HTML report: $REPORT"
  ls -lh "$REPORT"
else
  echo "Final HTML report was not found or is empty"
  exit 1
fi
```

Confirm that the final HTML contains the AMR & phylogeny sections:

```bash
grep -q 'id="amr"' "$REPORT" && \
  echo "AMR section present"

grep -q 'id="phylogeny"' "$REPORT" && \
  echo "Phylogeny section present"
```

For a phylogeny-eligible dataset, check for embedded tree content:

```bash
grep -Eq \
'class="tree-img"|data:image/png;base64|<svg' \
"$REPORT" && \
  echo "Rendered phylogeny content present"
```

If no species group meets the minimum same-species sample & reference requirements, the phylogeny section remains in the report & states why no tree was generated.

Open the report on macOS:

```bash
open "$REPORT"
```

Locate sample-level & final outputs:

```bash
find cromwell-executions/rMAP_Candida \
  -path '*/call-ASSEMBLY/shard-*/execution/*.contigs.fasta' \
  -type f -size +0c -print

find cromwell-executions/rMAP_Candida \
  -path '*/call-AMR/shard-*/execution/*.fungal_amr.summary.tsv' \
  -type f -size +0c -print

find cromwell-executions/rMAP_Candida \
  -type f \
  \( -name '*.treefile' -o -name '*.nwk' -o -name '*.png' \) \
  -size +0c -print

find cromwell-executions/rMAP_Candida \( \
  -name 'rMAP_Candida_summary.tsv' -o \
  -name 'rMAP_Candida_surveillance_summary.tsv' -o \
  -name 'rMAP_Candida_pairwise_snp_distances.tsv' \
\) -print
```

A complete validation should demonstrate that rMAP-Candida starts with paired raw FASTQ files, completes *de novo* assembly & all enabled downstream analyses & produces a non-empty integrated HTML report. AMR findings are genomic screening evidence & should be interpreted with species identity, validated marker catalogues, clinical metadata & phenotypic antifungal susceptibility testing where appropriate.

---

## Repository layout

```text
rMAP-Candida/
├── .gitignore
├── README.md
├── rMAP_Candida.wdl
│
├── docker/
│   ├── amr_chroquetas/
│   ├── candida_refs/
│   └── kraken2_bracken/
│
├── example/
│   ├── cromwell.local.fifo_portable.conf
│   ├── cromwell.options.no_docker_hash_lookup.json
│   ├── download_two_sample_test_data.sh
│   ├── rMAP_Candida.inputs.example.json
│   ├── rMAP_Candida_accessions.tsv
│   └── surveillance_metadata.tsv
│
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
| Assembly and assembly QC | `wdl/modules/assembly_qc.wdl` | MEGAHIT, QUAST, Compleasm/BUSCO | *De novo* assembly, contiguity metrics, optional completeness | contigs FASTA, QUAST TSV, Compleasm/BUSCO summaries |
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

## Containers, databases & public access

All major steps are Dockerized. Users should pull & inspect the containers before running the workflow.

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

The database build recipe & source-list notes should be maintained in:

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

rMAP-Candida does **not** combine different *Candida* spp. into one core-SNP tree. It first groups isolates by species & only builds a tree when enough same-species samples & a configured reference are available.

This design is deliberate because *Candida* spp. & related yeasts differ in ploidy, genome plasticity, heterozygosity, copy-number variation, aneuploidy, loss of heterozygosity & recombination behavior. The workflow therefore reports phylogeny eligibility explicitly & treats SNP distances as descriptive genomic relatedness metrics, not direct transmission proof.

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

WDL/Cromwell was selected because it provides a portable workflow description, explicit typed inputs & outputs, scatter-based parallelism for sample-level tasks, Docker runtime declarations & execution provenance. The same workflow can be run locally for a two-sample test, on a workstation or HPC environment for larger batches, or on WDL-compatible cloud platforms.

This implementation is designed to reduce installation barriers: users do not need to install each bioinformatics tool manually, because task-level dependencies are containerized & versioned.

---

## Troubleshooting

### Local execution requires the supplied validated configuration

For the validated local Docker/Colima route, use Cromwell 91 together with:

```text
example/cromwell.local.fifo_portable.conf
```

Run:

```bash
java   -Xms1g   -Xmx3g   -XX:+UseG1GC   -Dconfig.file=example/cromwell.local.fifo_portable.conf   -jar ~/cromwell-91.jar   run rMAP_Candida.wdl   --inputs example/rMAP-Candida.inputs.example.two_samples.json
```

Do not omit `-Dconfig.file` when following the validated local instructions. The supplied configuration handles the Local/SFS backend, task concurrency, Docker resource limits, container-visible task scripts & native `/tmp` FIFO handling.

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

If you use **rMAP-Candida** or any component of this workflow in your research, please cite the publication describing the workflow, the software repository & the major tools used.

### Primary publication

Mboowa G, Sserwadda I, Kanyerezi S, Kidenya BR, Bwambale J, Musinguzi B. **rMAP-Candida: a modular Dockerized WDL/Cromwell workflow for reproducible Candida species typing, assembly-contiguity assessment, antifungal-resistance marker screening, and phylogenomic surveillance.** *Frontiers in Bioinformatics*. 2026;6. doi:10.3389/fbinf.2026.1896572.

Full article: https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1896572/full

### Software repository

**rMAP-Candida: Rapid Mycological Analysis Pipeline for Candida genomic surveillance.**
https://github.com/gmboowa/rMAP-Candida

### Major tools

Please also cite the major software and resources used by the workflow, as applicable to the modules you enable:

* WDL / Cromwell
* Docker
* fastp
* FastQC
* MultiQC
* Kraken2
* Bracken
* MEGAHIT
* QUAST
* Compleasm and/or BUSCO
* ChroQueTas / FungAMR-derived resources
* Snippy
* IQ-TREE2
* ETE3

When citing results generated by a specific module, please cite the corresponding software/resource in addition to the primary **rMAP-Candida** publication.
