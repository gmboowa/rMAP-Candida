**rMAP-Candida: Rapid Mycological Analysis Pipeline for *Candida* genomic surveillance**

**rMAP-Candida** is a Dockerized WDL/Cromwell workflow for paired-end *Candida* whole-genome sequencing analysis. The current working workflow supports read QC, optional trimming, *Candida*-focused species typing, MEGAHIT genome assembly, QUAST assembly assessment, BUSCO fungal completeness assessment, ChroQueTas/FungAMR-based antifungal-resistance marker screening & an integrated HTML surveillance report.

The workflow is designed for **research, training & public health genomic surveillance**. Antifungal-resistance results should be interpreted as **genomic screening evidence**, not as direct clinical susceptibility calls, unless validated against local phenotypic antifungal susceptibility testing & appropriate clinical-laboratory standards.

---

## Current stable workflow

Use the current patched workflow:

```text
rMAP_Candida.wdl
```

This version preserves the working AMR detection & reporting logic while adding hardened BUSCO parsing & reporting.

### Major fixes in this version

- Preserves the working ChroQueTas/FungAMR marker parser.
- Correctly reports curated AMR marker hits such as `Cyp51 Y132F` & `Cyp51 V125A`.
- Keeps `No marker detected` separate from susceptibility interpretation.
- Displays detected `Gene / status` entries using a deep-red visual emphasis in the integrated HTML report.
- Hardens BUSCO execution & parsing using a *Candida*-appropriate lineage strategy.
- Reports BUSCO status & notes instead of silently masking BUSCO failures as missing values.
- Generates a clean integrated HTML report with summary cards, tables & provenance-oriented output sections.

---

## Workflow overview

```text
Paired-end FASTQ files
        │
        ├── 1. Read trimming & QC
        │       └── fastp HTML/JSON summaries
        │
        ├── 2. Species typing
        │       └── Kraken2 + Bracken with custom Candida-focused database
        │
        ├── 3. Genome assembly
        │       └── MEGAHIT contigs & assembly summary
        │
        ├── 4. Assembly quality assessment
        │       └── QUAST transposed per-sample report
        │
        ├── 5. Genome completeness assessment
        │       └── BUSCO fungal completeness table
        │
        ├── 6. Antifungal-resistance marker screening
        │       └── ChroQueTas/FungAMR marker-level evidence
        │
        └── 7. Integrated report
                ├── HTML report
                └── summary TSV
```

---

## Main workflow modules

| Module | Tool / method | Main purpose | Key outputs |
|---|---|---|---|
| Read trimming & QC | `fastp` | Adapter/quality trimming & read-level QC | `*.fastp.html`, `*.fastp.json`, trimmed FASTQs |
| Species typing | Kraken2 + Bracken | *Candida*-focused taxonomic classification | `*.kraken2.report`, `*.bracken.tsv`, `*.top_species.tsv` |
| Assembly | MEGAHIT | De novo assembly from paired-end reads | `*.contigs.fasta`, `*.assembly_summary.tsv` |
| Assembly QC | QUAST | Assembly continuity & GC metrics | `*.quast.report.tsv`, QUAST report folder |
| Completeness | BUSCO | Fungal genome completeness | `*.busco.summary.tsv`, BUSCO status & note columns |
| AMR marker screening | ChroQueTas/FungAMR | Curated antifungal-resistance marker detection | `*.fungal_amr.summary.tsv`, `*.fungal_amr.raw.tsv`, `*.fungal_amr.html` |
| Reporting | Built-in WDL report task | Integrated surveillance report | `rMAP_Candida_report.html`, `rMAP_Candida_summary.tsv` |

---

## Inputs

The workflow expects paired-end FASTQ files & matching sample names.

### Required inputs

| Input | Type | Description |
|---|---|---|
| `rMAP_Candida.sample_names` | `Array[String]` | Sample IDs. Must match the order of `read1s` & `read2s`. |
| `rMAP_Candida.read1s` | `Array[File]` | R1 FASTQ files, usually `*_1.fastq.gz` or `*_R1.fastq.gz`. |
| `rMAP_Candida.read2s` | `Array[File]` | R2 FASTQ files, usually `*_2.fastq.gz` or `*_R2.fastq.gz`. |

### Recommended toggles

| Input | Recommended value | Description |
|---|---:|---|
| `rMAP_Candida.do_trimming` | `true` | Run fastp trimming & QC. |
| `rMAP_Candida.do_quality_control` | `true` | Keep read QC enabled. |
| `rMAP_Candida.do_species_typing` | `true` | Run Kraken2/Bracken species typing. |
| `rMAP_Candida.do_assembly` | `true` | Run MEGAHIT assembly. |
| `rMAP_Candida.do_assembly_qc` | `true` | Run QUAST assembly assessment. |
| `rMAP_Candida.do_busco` | `true` | Run BUSCO completeness assessment. |
| `rMAP_Candida.do_fungal_amr` | `true` | Run ChroQueTas/FungAMR AMR screening. |

### Container & database inputs

| Input | Recommended value | Description |
|---|---|---|
| `rMAP_Candida.fungal_kraken2_bracken_docker` | `gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db` | Custom Kraken2/Bracken container with *Candida*-focused database. |
| `rMAP_Candida.kraken_db_path` | `/opt/kraken2_db/candida` | Path to the bundled Kraken2 database inside the container. |
| `rMAP_Candida.fungamr_docker` | `gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed` | AMR container with ChroQueTas/FungAMR support & required dependencies. |
| `rMAP_Candida.bracken_level` | `S` | Species-level Bracken abundance estimation. |
| `rMAP_Candida.max_cpus` | `8` | Maximum CPUs to use locally. |
| `rMAP_Candida.max_memory_gb` | `14` | Maximum memory in GB for local Cromwell execution. |

---

## Example input JSON

```json
{
  "rMAP_Candida.sample_names": [
    "SRR21675877",
    "SRR12073482",
    "SRR21675900"
  ],
  "rMAP_Candida.read1s": [
    "~/SRR21675877_1.fastq.gz",
    "~/SRR12073482_1.fastq.gz",
    "~/SRR21675900_1.fastq.gz"
  ],
  "rMAP_Candida.read2s": [
    "~/SRR21675877_2.fastq.gz",
    "~/SRR12073482_2.fastq.gz",
    "~/SRR21675900_2.fastq.gz"
  ],
  "rMAP_Candida.do_trimming": true,
  "rMAP_Candida.do_quality_control": true,
  "rMAP_Candida.do_species_typing": true,
  "rMAP_Candida.do_assembly": true,
  "rMAP_Candida.do_assembly_qc": true,
  "rMAP_Candida.do_busco": true,
  "rMAP_Candida.do_fungal_amr": true,
  "rMAP_Candida.fungal_kraken2_bracken_docker": "gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db",
  "rMAP_Candida.kraken_db_path": "/opt/kraken2_db/candida",
  "rMAP_Candida.fungamr_docker": "gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed",
  "rMAP_Candida.bracken_level": "S",
  "rMAP_Candida.max_cpus": 8,
  "rMAP_Candida.max_memory_gb": 14
}
```

> Important: Ensure every R1 file has the correct matching R2 file. A common mistake is accidentally entering an R1 file in the `read2s` list.

---

## Running the workflow locally with Cromwell

From the directory containing your input JSON:

```bash
java -jar ~/cromwell-92.jar run  rMAP_Candida.wdl --inputs rMAP-Candida.inputs.json
```

For long local runs, use `tmux` so that the workflow continues even if your terminal disconnects:

```bash
tmux new -s candida_run
java -jar ~/cromwell-92.jar run ~/rMAP_Candida.wdl --inputs rMAP-Candida.inputs.json
```
---

## Expected final outputs

Cromwell stores task outputs under:

```text
cromwell-executions/rMAP_Candida/<workflow-id>/
```

The final integrated report is produced by the report merge task:

```text
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_report.html
```

The final summary table is:

```text
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_summary.tsv
```

---

## Integrated HTML report sections

The current working report contains the following sections:

1. **Executive summary**
   - Number of samples analyzed
   - Top species groups
   - Total AMR hits
   - Median N50
   - Interpretation note for genomic AMR screening

2. **Sample-level surveillance summary**
   - Species assignment
   - Species read percentage
   - Contig count
   - N50
   - BUSCO completeness
   - AMR hit count

3. **Candida species typing using Kraken2/Bracken**
   - Top species call
   - Read percentage
   - Clade reads
   - Taxon reads
   - TaxID
   - Evidence source

4. **MEGAHIT assembly summary**
   - Number of contigs
   - Total assembly length
   - N50
   - Largest contig

5. **Assembly quality assessment with QUAST**
   - Contig count
   - Largest contig
   - Total length
   - GC percentage
   - N50

6. **BUSCO fungal completeness**
   - Complete BUSCO percentage
   - Single-copy BUSCO percentage
   - Duplicated BUSCO percentage
   - Fragmented BUSCO percentage
   - Missing BUSCO percentage
   - Number of BUSCO markers
   - BUSCO status
   - BUSCO note

7. **Fungal antifungal-resistance characterization**
   - Sample
   - Species
   - Drug class
   - Drug
   - Gene / status
   - Mutation
   - Effect
   - Evidence level
   - Interpretation

8. **Output navigation & provenance**
   - Notes on where per-sample QC, assembly, QUAST, BUSCO, AMR & report files are available.

---

## BUSCO logic

The workflow runs BUSCO per assembly & uses a *Candida*-appropriate fungal lineage strategy.

Current BUSCO behavior:

1. Uses `saccharomycetes_odb10` as the preferred default for *Candida* assemblies.
2. Retries using broader fungal lineages such as `ascomycota_odb10` or `fungi_odb10` when needed.
3. Parses BUSCO output into a standardized per-sample summary table.
4. Reports real BUSCO status & notes when results are missing or incomplete.

The integrated report should show values such as:

```text
Complete_BUSCO_%
Single_copy_BUSCO_%
Duplicated_BUSCO_%
Fragmented_BUSCO_%
Missing_BUSCO_%
BUSCO markers
Status
BUSCO note
```

Recommended interpretation:

| BUSCO completeness | Suggested interpretation |
|---:|---|
| `>= 95%` | Strong completeness for many *Candida* WGS assemblies |
| `90–94.9%` | Generally acceptable, but inspect assembly metrics |
| `< 90%` | Review read quality, assembly fragmentation, species identity & lineage selection |
| `NA` | Inspect BUSCO logs, status & note columns |

---

## Antifungal-resistance marker logic

The workflow uses ChroQueTas/FungAMR-derived outputs to identify curated antifungal-resistance markers. The current parser is designed to capture exact marker rows from ChroQueTas output files & summarize them into the integrated report.

### Important interpretation rule

`No marker detected` means only that the configured AMR scanner did not detect a curated marker in the harvested output. It **does not** mean that the isolate is susceptible.

Resistance can arise through mechanisms such as:

- `ERG11` / `Cyp51` alterations
- `TAC1`, `UPC2`, `MRR1`, or `PDR1` regulatory changes
- Efflux-mediated mechanisms
- `FKS` alterations
- Flucytosine-associated markers such as `FCY1`, `FCY2`, or `FUR1`
- Copy-number changes
- Aneuploidy or loss of heterozygosity
- Species-specific mechanisms not represented in the installed marker database

### Example validated positive-control behavior

The current workflow has successfully reported marker-level AMR evidence in *Candidozyma auris* validation samples, including:

| Sample | Species | Gene / status | Mutation | Effect | Evidence level |
|---|---|---|---|---|---|
| `SRR21675877` | *Candidozyma auris* | `Cyp51` | `Y132F` | `FungAMR MUTATION` | `FungAMR curated marker detected by ChroQueTas` |
| `SRR21675900` | *Candidozyma auris* | `Cyp51` | `V125A` | `FungAMR MUTATION` | `FungAMR curated marker detected by ChroQueTas` |
| `SRR12073482` | *Candidozyma auris* | `Cyp51` | `V125A` | `FungAMR MUTATION` | `FungAMR curated marker detected by ChroQueTas` |

These are useful workflow-validation controls because they test whether the AMR module can detect & report curated ChroQueTas/FungAMR markers in the integrated HTML report.

---

## Recommended validation checks after a run

### 1. Confirm Cromwell completed successfully

```bash
grep -i "workflow.*succeeded" cromwell-*.log
```

Or inspect the end of the terminal log for a succeeded status.

### 2. Confirm the integrated report exists

```bash
find cromwell-executions/rMAP_Candida \
  -name "rMAP_Candida_report.html" \
  -print
```

### 3. Confirm the summary TSV exists

```bash
find cromwell-executions/rMAP_Candida \
  -name "rMAP_Candida_summary.tsv" \
  -print
```

### 4. Check AMR summaries

```bash
find cromwell-executions/rMAP_Candida \
  -name "*.fungal_amr.summary.tsv" \
  -print
```

Inspect a positive-control sample:

```bash
cat cromwell-executions/rMAP_Candida/<workflow-id>/call-AMR/shard-0/execution/*.fungal_amr.summary.tsv
```

Expected columns:

```text
sample_id
species
drug_class
drug
gene_or_status
mutation
effect
evidence_level
interpretation
```

### 5. Check BUSCO summaries

```bash
find cromwell-executions/rMAP_Candida \
  -name "*.busco.summary.tsv" \
  -print
```

If BUSCO values are `NA`, inspect the corresponding BUSCO task execution folder & logs.

---

## Troubleshooting

### AMR results show `No marker detected`

This may be expected for some isolates. It does not mean susceptible.

Check whether ChroQueTas actually ran & whether marker files were generated:

```bash
find cromwell-executions/rMAP_Candida \
  -path "*/call-AMR/*/execution/amr_out/*" \
  -type f | sort
```

Review scanner logs:

```bash
cat cromwell-executions/rMAP_Candida/<workflow-id>/call-AMR/shard-0/execution/SAMPLE.fungal_amr.log
cat cromwell-executions/rMAP_Candida/<workflow-id>/call-AMR/shard-0/execution/amr_out/SAMPLE.ChroQueTas.stdout
cat cromwell-executions/rMAP_Candida/<workflow-id>/call-AMR/shard-0/execution/amr_out/SAMPLE.ChroQueTas.stderr
```

### AMR task reports stale output directory

This usually means a ChroQueTas output folder already exists inside the task execution directory.

The patched workflow is designed to remove or isolate stale AMR output folders during execution. If this still occurs during manual testing, remove the manual output directory before rerunning:

```bash
rm -rf manual_amr_test
```

### `miniprot is required & not installed`

Use the fixed AMR container:

```text
gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed
```

Do not use the older AMR container for the current workflow unless it has been rebuilt with `miniprot` & the required ChroQueTas dependencies.

### BUSCO values remain `NA`

Check the BUSCO status & note columns in the report & summary TSV. Then inspect the task logs:

```bash
find cromwell-executions/rMAP_Candida \
  -path "*/call-BUSCO/*/execution/*" \
  -type f | sort
```

Common causes include:

- Incomplete assembly
- Wrong or unavailable BUSCO lineage database
- Memory limit too low
- Docker image issue
- Interrupted BUSCO download or local cache issue

### Kraken2/Bracken database errors

Confirm the expected bundled database path:

```text
~/kraken2_db/candida
```

Expected files inside the Kraken2 database include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

### Paired FASTQ input errors

Confirm that each sample has a matched R1 and R2:

```bash
ls -lh *_1.fastq.gz *_2.fastq.gz
```

For each sample:

```text
SAMPLE_1.fastq.gz  -> read1s
SAMPLE_2.fastq.gz  -> read2s
```

---

## Recommended quality thresholds

These are practical starting thresholds & should be adapted to the study design, sequencing depth, species & laboratory context.

| Metric | Threshold |
|---|---:|
| Species dominance by Kraken2/Bracken | Preferably `>80%` for single-isolate WGS |
| Assembly size | Species-dependent; inspect for major deviations |
| Contig count | Lower is generally better |
| N50 | Higher is generally better |
| BUSCO completeness | Preferably `>90%`; `>95%` is stronger |
| AMR marker detection | Treat as screening evidence, not a phenotype |

---

## Repository structure

Repository layout:

```text
rMAP_Candida/
├── README.md
├── LICENSE
├── workflows/
│   └── rMAP_Candida.wdl
├── inputs/
│   └── rMAP-Candida.inputs.example.json
├── docs/
│   ├── index.html
│   └── reports/
│       └── example_report.html
├── docker/
│   ├── kraken2_bracken/
│   └── amr_chroquetas/
├── test_data/
│   └── README.md
└── examples/
    ├── positive_controls.md
    └── run_local_cromwell.sh
```

---

## Limitations

- The pipeline is not a standalone clinical diagnostic system.
- Genomic AMR markers may not fully predict phenotypic antifungal susceptibility.
- AMR interpretation depends on the completeness & accuracy of the installed marker database.
- Copy-number variation, aneuploidy, loss of heterozygosity, promoter changes & regulatory mechanisms may not be fully captured.
- Species identification depends on the completeness & quality of the Kraken2/Bracken database.
- Mixed cultures, contamination, or low sequencing depth can affect interpretation.
- BUSCO completeness depends on assembly quality & lineage selection.
- Results should be validated before use in clinical decision-making.

---

## Suggested citation

If you use this workflow, cite the repository & the major tools used by the workflow.


```text
rMAP-Candida: Rapid Mycological Analysis Pipeline for Candida genomic surveillance.
```

Also cite the relevant tools used in the workflow, including Cromwell/WDL, fastp, Kraken2, Bracken, MEGAHIT, QUAST, BUSCO, ChroQueTas & FungAMR/CARD resources where applicable.

---

## Disclaimer

`rMAP-Candida` is intended for research, training & public health genomic surveillance. It is not intended for direct clinical decision-making unless the workflow, databases, reporting logic & interpretation rules have been validated under appropriate clinical laboratory, regulatory, biosafety, and quality-management frameworks.

