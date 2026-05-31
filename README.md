## rMAP-Candida: Rapid Mycological Analysis Pipeline for *Candida spp* Genomic Surveillance

**rMAP-Candida** is a Dockerized WDL/Cromwell workflow for paired-end *Candida spp.* whole-genome sequencing analysis & genomic surveillance. The current surveillance-ready workflow supports read trimming & QC, *Candida*-focused species typing, MEGAHIT assembly, QUAST assembly assessment, fast genome completeness assessment using **Compleasm**, ChroQueTas/FungAMR-based antifungal-resistance marker screening, species-aware core-SNP phylogeny, SNP-distance summaries, optional surveillance metadata integration & an integrated HTML surveillance report.

The workflow is designed for **research, training, outbreak investigation support & public health genomic surveillance**. Antifungal-resistance results should be interpreted as **genomic screening evidence**, not as direct clinical susceptibility calls, unless validated against local phenotypic antifungal susceptibility testing & appropriate clinical-laboratory standards.

---

## Current stable workflow

Use the current surveillance-ready workflow:

```text
rMAP_Candida.wdl
```
---

BUSCO is informative but slow & can take several hours for only a few *Candida spp* samples.

The current workflow now uses:

```text
Genome completeness assessment using Compleasm
```

Compleasm provides BUSCO-like conserved ortholog completeness metrics but is much faster & better suited for routine surveillance-style runs.

Default completeness settings:

```json
"rMAP_Candida.do_busco": false,
"rMAP_Candida.do_compleasm": true,
"rMAP_Candida.compleasm_docker": "huangnengcsu/compleasm:v0.2.7",
"rMAP_Candida.compleasm_lineage": "saccharomycetes",
"rMAP_Candida.compleasm_odb": "odb12"
```

### Surveillance-ready reporting

The integrated HTML report now includes:

- A **Surveillance Readiness Dashboard**
- Species confidence & mixed-sample warnings
- Assembly QC interpretation
- Compleasm completeness interpretation
- AMR marker interpretation
- Phylogeny eligibility & inclusion status
- Species-aware SNP distance & closest-neighbor summaries
- Optional surveillance metadata table
- Sticky table of contents
- Collapsible sample-level cards
- Downloadable TSV outputs for surveillance interpretation

### Optional metadata support

The workflow can now accept an optional surveillance metadata TSV:

```wdl
File? surveillance_metadata_tsv
```

Recommended metadata columns:

```text
sample_id
country
site
collection_date
specimen_type
patient_group
ward_or_facility
sequencing_platform
species
metadata_note
```

If no metadata file is provided, the workflow still runs & the report notes that no surveillance metadata was supplied.

---

## Workflow overview

```text
Paired-end FASTQ files
        │
        ├── 1. Read trimming & QC
        │       └── fastp HTML/JSON summaries & trimmed FASTQs
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
        │       └── Compleasm completeness table
        │
        ├── 6. Antifungal-resistance marker screening
        │       └── ChroQueTas/FungAMR marker-level evidence
        │
        ├── 7. Species-aware core-SNP phylogeny
        │       ├── Reference-based alignment
        │       ├── IQ-TREE phylogeny
        │       ├── ETE3 rendered tree
        │       └── Newick tree outputs
        │
        ├── 8. SNP-distance & cluster-support summaries
        │       └── Pairwise SNP distance matrix & closest-neighbor table
        │
        └── 9. Integrated surveillance report
                ├── rMAP_Candida_report.html
                ├── rMAP_Candida_summary.tsv
                ├── rMAP_Candida_surveillance_summary.tsv
                └── rMAP_Candida_pairwise_snp_distances.tsv
```

---

## Main workflow modules

| Module | Tool / method | Main purpose | Key outputs |
|---|---|---|---|
| Read trimming & QC | `fastp` | Adapter/quality trimming & read-level QC | `*.fastp.html`, `*.fastp.json`, trimmed FASTQs |
| Species typing | Kraken2 + Bracken | *Candida*-focused taxonomic classification | `*.kraken2.report`, `*.bracken.tsv`, `*.top_species.tsv` |
| Assembly | MEGAHIT | De novo assembly from paired-end reads | `*.contigs.fasta`, `*.assembly_summary.tsv` |
| Assembly QC | QUAST | Assembly continuity & GC metrics | `*.quast.report.tsv`, QUAST report folder |
| Completeness | Compleasm | Fast conserved ortholog genome completeness assessment | `*.compleasm.summary.tsv`, status & note columns |
| AMR marker screening | ChroQueTas/FungAMR | Curated antifungal-resistance marker detection | `*.fungal_amr.summary.tsv`, `*.fungal_amr.raw.tsv`, `*.fungal_amr.html` |
| Phylogeny | Snippy/BCFtools/IQ-TREE/ETE3 logic in WDL | Species-aware core-SNP phylogenetic analysis | core SNP alignment, Newick tree, rendered tree image |
| SNP distances | Built-in report logic | Pairwise SNP distance & closest-neighbor summary | `rMAP_Candida_pairwise_snp_distances.tsv` |
| Reporting | Built-in WDL report task | Integrated genomic surveillance report | `rMAP_Candida_report.html`, surveillance summary TSVs |

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
| `rMAP_Candida.do_busco` | `false` | Legacy BUSCO module is disabled by default. |
| `rMAP_Candida.do_compleasm` | `true` | Run Compleasm genome completeness assessment. |
| `rMAP_Candida.do_fungal_amr` | `true` | Run ChroQueTas/FungAMR AMR screening. |
| `rMAP_Candida.do_phylogeny` | `true` | Run species-aware phylogenetic analysis where enough samples are available. |
| `rMAP_Candida.render_phylogeny_tree` | `true` | Render publication-style tree image for the HTML report. |

### Optional surveillance metadata input

| Input | Type | Description |
|---|---|---|
| `rMAP_Candida.surveillance_metadata_tsv` | `File?` | Optional metadata table with sample-level country, site, date, specimen type & other surveillance fields. |

### Container & database inputs

| Input | Recommended value | Description |
|---|---|---|
| `rMAP_Candida.fungal_kraken2_bracken_docker` | `gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db` | Custom Kraken2/Bracken container with *Candida*-focused database. |
| `rMAP_Candida.kraken_db_path` | `/opt/kraken2_db/candida` | Path to the bundled Kraken2 database inside the container. |
| `rMAP_Candida.fungamr_docker` | `gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed` | AMR container with ChroQueTas/FungAMR support & required dependencies. |
| `rMAP_Candida.compleasm_docker` | `huangnengcsu/compleasm:v0.2.7` | Compleasm container. |
| `rMAP_Candida.compleasm_lineage` | `saccharomycetes` | Recommended Compleasm lineage for most *Candida* assemblies. |
| `rMAP_Candida.compleasm_odb` | `odb12` | Compleasm ortholog database version. |
| `rMAP_Candida.candida_refs_docker` | `gmboowa/rmap-candida-refs:2026.05` | Container containing curated species reference genomes for phylogeny. |
| `rMAP_Candida.candida_refs_manifest` | `/opt/rmap_candida_refs/references.tsv` | Manifest of available reference genomes in the reference container. |
| `rMAP_Candida.tree_visualization_docker` | `gmboowa/ete3-render:1.18` | ETE3 tree rendering container. |
| `rMAP_Candida.bracken_level` | `S` | Species-level Bracken abundance estimation. |
| `rMAP_Candida.max_cpus` | `8` | Maximum CPUs to use locally. |
| `rMAP_Candida.max_memory_gb` | `14` | Maximum memory in GB for local Cromwell execution. |

---

## Example input JSON

```json
{
  "rMAP_Candida.sample_names": [
    "SRR30894132",
    "SRR30894141",
    "SRR36120464",
    "SRR30894125"
  ],
  "rMAP_Candida.read1s": [
    "~/SRR30894132_1.fastq.gz",
    "~/SRR30894141_1.fastq.gz",
    "~/SRR36120464_1.fastq.gz",
    "~/SRR30894125_1.fastq.gz"
  ],
  "rMAP_Candida.read2s": [
    "~/SRR30894132_2.fastq.gz",
    "~/SRR30894141_2.fastq.gz",
    "~/SRR36120464_2.fastq.gz",
    "~/SRR30894125_2.fastq.gz"
  ],
  "rMAP_Candida.do_trimming": true,
  "rMAP_Candida.do_quality_control": true,
  "rMAP_Candida.do_species_typing": true,
  "rMAP_Candida.do_assembly": true,
  "rMAP_Candida.do_assembly_qc": true,
  "rMAP_Candida.do_busco": false,
  "rMAP_Candida.do_compleasm": true,
  "rMAP_Candida.do_fungal_amr": true,
  "rMAP_Candida.do_phylogeny": true,
  "rMAP_Candida.render_phylogeny_tree": true,
  "rMAP_Candida.fungal_kraken2_bracken_docker": "gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db",
  "rMAP_Candida.kraken_db_path": "/opt/kraken2_db/candida",
  "rMAP_Candida.fungamr_docker": "gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed",
  "rMAP_Candida.compleasm_docker": "huangnengcsu/compleasm:v0.2.7",
  "rMAP_Candida.compleasm_lineage": "saccharomycetes",
  "rMAP_Candida.compleasm_odb": "odb12",
  "rMAP_Candida.candida_refs_docker": "gmboowa/rmap-candida-refs:2026.05",
  "rMAP_Candida.candida_refs_manifest": "/opt/rmap_candida_refs/references.tsv",
  "rMAP_Candida.bracken_level": "S",
  "rMAP_Candida.min_species_samples_for_tree": 3,
  "rMAP_Candida.snippy_phylogeny_species": [
    "Candida albicans"
  ],
  "rMAP_Candida.haploid_phylogeny_species": [
    "Nakaseomyces glabratus"
  ],
  "rMAP_Candida.min_depth_for_phylogeny": 10,
  "rMAP_Candida.min_mapping_quality": 20,
  "rMAP_Candida.min_base_quality_for_phylogeny": 20,
  "rMAP_Candida.min_variant_quality_for_phylogeny": 20,
  "rMAP_Candida.core_site_min_fraction": 0.95,
  "rMAP_Candida.iqtree2_model": "GTR+G",
  "rMAP_Candida.iqtree2_bootstraps": 1000,
  "rMAP_Candida.max_cpus": 8,
  "rMAP_Candida.max_memory_gb": 14,
  "rMAP_Candida.tree_visualization_docker": "gmboowa/ete3-render:1.18"
}
```

To include optional metadata, add:

```json
"rMAP_Candida.surveillance_metadata_tsv": "~/surveillance_metadata.tsv"
```

---

## Example surveillance metadata TSV

```text
sample_id	country	site	collection_date	specimen_type	patient_group	ward_or_facility	sequencing_platform	species	metadata_note
SRR30894132	not_available	not_available	not_available	not_available	not_available	not_available	not_available	Candida albicans	Verify public BioSample/SRA/ENA metadata before final interpretation.
SRR30894141	not_available	not_available	not_available	not_available	not_available	not_available	not_available	Candida albicans	Verify public BioSample/SRA/ENA metadata before final interpretation.
SRR36120464	not_available	not_available	not_available	not_available	not_available	not_available	not_available	Candida albicans	Verify public BioSample/SRA/ENA metadata before final interpretation.
SRR30894125	not_available	not_available	not_available	not_available	not_available	not_available	not_available	Candida albicans	Verify public BioSample/SRA/ENA metadata before final interpretation.
```

---

## Running the workflow locally with Cromwell

From the directory containing your WDL & input JSON:

```bash
java -jar ~/cromwell-92.jar \
  run ~/rMAP_Candida.wdl \
  --inputs ~/rMAP-Candida.inputs.json \
  --options ~/cromwell.options.no_docker_hash_lookup.json
```

For long local runs, use `tmux` so that the workflow continues even if your terminal disconnects:

```bash
tmux new -s candida_run
java -jar ~/cromwell-92.jar \
  run ~/rMAP_Candida.wdl \
  --inputs ~/rMAP-Candida.inputs.json \
  --options ~/cromwell.options.no_docker_hash_lookup.json
```

Detach safely with:

```bash
Ctrl-b d
```

Check the session later:

```bash
tmux ls
tmux attach -t candida_run
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

The main summary tables are:

```text
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_summary.tsv
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_surveillance_summary.tsv
call-MERGE_MYC_REPORTS/execution/rMAP_Candida_pairwise_snp_distances.tsv
```

Depending on enabled modules & available samples, additional outputs may include:

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

## Integrated HTML report sections

The current surveillance-ready report contains the following sections:

1. **Executive summary**
   - Number of samples analyzed
   - Top species groups
   - Total AMR hits
   - Median N50
   - Compleasm completeness overview
   - Interpretation note for genomic AMR screening

2. **Surveillance readiness dashboard**
   - Species confidence
   - Assembly QC status
   - Compleasm status
   - AMR marker status
   - Phylogeny eligibility
   - Overall surveillance interpretation

3. **Sample-level surveillance summary**
   - Species assignment
   - Species read percentage
   - Contig count
   - N50
   - Compleasm completeness
   - AMR hit count
   - Collapsible detailed interpretation

4. **Surveillance metadata**
   - Country
   - Site
   - Collection date
   - Specimen type
   - Patient group
   - Ward/facility
   - Sequencing platform
   - Metadata notes

5. ***Candida* species typing using Kraken2/Bracken**
   - Top species call
   - Species confidence category
   - Read percentage
   - Clade reads
   - Taxon reads
   - TaxID
   - Evidence source

6. **MEGAHIT assembly summary**
   - Number of contigs
   - Total assembly length
   - N50
   - Largest contig

7. **Assembly quality assessment with QUAST**
   - Contig count
   - Largest contig
   - Total length
   - GC percentage
   - N50
   - Assembly QC interpretation

8. **Genome completeness assessment using Compleasm**
   - Complete percentage
   - Single-copy percentage
   - Duplicated percentage
   - Fragmented percentage
   - Missing percentage
   - Number of markers
   - Compleasm status
   - Compleasm note

9. **Fungal antifungal-resistance characterization**
   - Sample
   - Species
   - Drug class
   - Drug
   - Gene/status
   - Mutation
   - Effect
   - Evidence level
   - Interpretation

10. **Species-aware core-SNP phylogeny**
    - Tree image
    - Newick output
    - Species group
    - Included samples
    - Excluded samples
    - Reference used
    - Notes on tree interpretation

11. **Phylogeny eligibility table**
    - Species
    - Species group size
    - Reference availability
    - Core genome usability
    - Included in tree
    - Inclusion/exclusion reason

12. **Species-aware SNP distance & cluster summary**
    - Pairwise SNP distances
    - Closest genetic neighbor
    - Conservative cluster-support flag
    - Interpretation note: low SNP distance suggests close relatedness but should be interpreted with epidemiological metadata

13. **Output navigation & provenance**
    - Tool/module notes
    - Output file locations
    - Parameters & containers used

---

## Compleasm logic

The workflow runs Compleasm per assembly using a *Candida*-appropriate lineage.

Current default behavior:

1. Uses `saccharomycetes` as the preferred lineage for most *Candida* assemblies.
2. Uses `odb12` as the ortholog database setting.
3. Parses Compleasm output into standardized per-sample summary columns.
4. Reports status & notes when results are missing or incomplete.
5. Keeps the HTML section title as:

```text
Genome completeness assessment using Compleasm
```

The integrated report should show values such as:

```text
Complete_%
Single_copy_%
Duplicated_%
Fragmented_%
Missing_%
Markers
Status
Compleasm note
```

Recommended interpretation:

| Compleasm completeness | Suggested interpretation |
|---:|---|
| `>= 95%` | Strong completeness for many *Candida* WGS assemblies |
| `90–94.9%` | Generally acceptable, but inspect assembly metrics |
| `< 90%` | Review read quality, assembly fragmentation, species identity & lineage selection |
| `NA` | Inspect Compleasm logs, status & note columns |

---

## Species confidence & mixed-sample warning

The report classifies species calls using the top-species read proportion from Kraken2/Bracken.

| Top-species abundance | Report interpretation |
|---:|---|
| `>= 95%` | High-confidence species assignment |
| `80–94.9%` | Moderate-confidence assignment; review for mixed content or lower purity |
| `< 80%` | Low-confidence or possible mixed sample; review before surveillance interpretation |

This helps flag samples that may represent mixed cultures, contamination, insufficient database coverage, or lower-confidence classification.

---

## Surveillance readiness interpretation

The report integrates several metrics into a sample-level surveillance status.

Example status logic:

| Metric | PASS / Ready | Review |
|---|---|---|
| Species confidence | Top species `>= 95%` | Top species `< 95%` |
| Assembly quality | Acceptable contig count and N50 | Very fragmented assembly or low N50 |
| Compleasm completeness | `>= 95%` | `< 95%` or `NA` |
| AMR screen | Marker result produced | AMR task failed or missing |
| Phylogeny | Included or clearly not applicable | Expected but excluded without clear reason |

Suggested report labels:

```text
Ready for surveillance interpretation
Review species assignment
Review assembly quality
Review completeness
AMR marker detected
Phylogeny excluded
```

---

## Species-aware phylogeny

The workflow supports species-aware core-SNP phylogeny. It should not combine unrelated *Candida* species into a single tree. It only attempts tree construction for species groups with enough samples.

Important inputs:

| Input | Purpose |
|---|---|
| `rMAP_Candida.min_species_samples_for_tree` | Minimum number of same-species samples required to build a tree. |
| `rMAP_Candida.snippy_phylogeny_species` | Species treated with diploid-aware or Snippy-based reference workflow. |
| `rMAP_Candida.haploid_phylogeny_species` | Species treated with haploid-style phylogeny assumptions where appropriate. |
| `rMAP_Candida.candida_refs_docker` | Container with curated species references. |
| `rMAP_Candida.candida_refs_manifest` | Reference manifest inside the reference container. |
| `rMAP_Candida.iqtree2_model` | IQ-TREE model, e.g. `GTR+G`. |
| `rMAP_Candida.iqtree2_bootstraps` | Bootstrap replicates. |
| `rMAP_Candida.render_phylogeny_tree` | Whether to render a PNG tree for the HTML report. |

Tree interpretation should be conservative:

```text
Low SNP distance or close phylogenetic clustering suggests genetic relatedness, but outbreak or transmission interpretation requires epidemiological metadata, sampling dates, facility information & validated species-specific thresholds.
```

---

## SNP distance output

The report includes a species-aware pairwise SNP distance & closest-neighbor summary.

Expected final TSV:

```text
rMAP_Candida_pairwise_snp_distances.tsv
```

Suggested interpretation:

| SNP-distance pattern | Interpretation |
|---|---|
| Low pairwise SNP distance | Possible close genetic relatedness; interpret with metadata |
| High pairwise SNP distance | Samples are likely more genetically distinct |
| Missing SNP distance | Tree/alignment was not generated or sample was excluded |

Avoid hard outbreak claims unless thresholds are validated for the species, sequencing approach & analysis workflow.

---

## Antifungal-resistance marker logic

The workflow uses ChroQueTas/FungAMR-derived outputs to identify curated antifungal-resistance markers. The current parser is designed to capture exact marker rows from ChroQueTas output files & summarize them into the integrated report.

### Important interpretation rule

`No marker detected` means only that the configured AMR scanner did not detect a curated marker in the harvested output. It **does not** mean that the isolate is susceptible.

Resistance can arise through mechanisms such as:

- `ERG11` / `Cyp51` alterations
- `TAC1`, `UPC2`, `MRR1`, or `PDR1` regulatory changes
- Efflux-mediated mechanisms
- `FKS1` / `FKS2` alterations
- Flucytosine-associated markers such as `FCY1`, `FCY2`, or `FUR1`
- Copy-number changes
- Aneuploidy or loss of heterozygosity
- Species-specific mechanisms not represented in the installed marker database

### Species-specific AMR interpretation examples

| Species | Important loci/mechanisms to consider |
|---|---|
| *Candida albicans* | `ERG11`, `TAC1`, `UPC2`, `MRR1`, `MDR1`, `CDR1/CDR2`, `FKS1` |
| *Candidozyma auris* | `ERG11`, `TAC1B`, `FKS1`, copy-number changes, aneuploidy/LOH |
| *Nakaseomyces glabratus* | `PDR1`, `CDR1/CDR2/SNQ2`, `FKS1/FKS2` |
| *Candida tropicalis* | `ERG11`, efflux-associated genes, `FKS1`, species-specific azole mechanisms |
| *Pichia kudriavzevii* | Intrinsic & acquired azole-related mechanisms require species-aware interpretation |

### Example positive-control behavior

Earlier workflow validation successfully reported marker-level AMR evidence in *Candidozyma auris* validation samples, including:

| Sample | Species | Gene/status | Mutation | Effect | Evidence level |
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

### 3. Confirm the summary TSVs exist

```bash
find cromwell-executions/rMAP_Candida \
  -name "rMAP_Candida_summary.tsv" \
  -o -name "rMAP_Candida_surveillance_summary.tsv" \
  -o -name "rMAP_Candida_pairwise_snp_distances.tsv" \
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

### 5. Check Compleasm summaries

```bash
find cromwell-executions/rMAP_Candida \
  -name "*.compleasm.summary.tsv" \
  -print
```

If Compleasm values are `NA`, inspect the corresponding Compleasm task execution folder & logs.

### 6. Check tree & SNP-distance outputs

```bash
find cromwell-executions/rMAP_Candida \
  \( -name "*.treefile" -o -name "*.nwk" -o -name "*tree*.png" -o -name "*snp*distance*.tsv" \) \
  -print
```
---
## Running rMAP-Candida Locally with Controlled Cromwell Shards

When running `rMAP-Candida` locally on a laptop or workstation, it is useful to limit the number of Cromwell jobs/shards running at the same time. This prevents too many Docker containers from running in parallel & reduces the risk of memory pressure, hanging tasks, or system slowdowns.

This workflow has been tested locally with Cromwell 92 using a reusable local backend configuration that limits concurrent jobs while allowing input files from any location under the user home directory.

### 1. Create a reusable Cromwell local configuration

Create a file called:

```bash
~/cromwell.local.reusable.shard2.no_memory_attr.conf

```
```bash
include required(classpath("application"))

system {
  max-concurrent-workflows = 1
}

backend {
  default = "Local"

  providers {
    Local {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

      config {
        # Controls the maximum number of Cromwell jobs/shards
        # running at the same time.
        concurrent-job-limit = 2

        runtime-attributes = """
        String? docker
        String? docker_user
        String? dockerWorkingDir
        Int cpu = 1
        Float memory_gb = 4
        String? disks
        """

        submit = """
        /bin/bash ${script}
        """

        submit-docker = """
        docker run --rm \
          -v ${cwd}:${docker_cwd} \
          -v /Users/gerald:/Users/gerald \
          -w ${docker_cwd} \
          ${docker} \
          /bin/bash ${script}
        """
      }
    }
  }
}
```
---
2. Why this configuration is useful

The key line is:

```bash

concurrent-job-limit = 2
```
This limits Cromwell to running only two jobs or shards at the same time. For example, if the workflow scatters over many samples, Cromwell will queue the remaining jobs and only execute two concurrently.

For very memory-limited runs, this can be reduced to:
```bash
concurrent-job-limit = 1
```
For larger machines, it can be increased, for example:
```bash
concurrent-job-limit = 3
```

or:

```bash
concurrent-job-limit = 4
```

3. Why the Docker mount is reusable

The Docker configuration uses:

-v /Users/gerald:/Users/gerald

This makes the configuration reusable across different projects, because any input files stored under the user home directory are visible inside Docker containers.

For example, the following paths will all be accessible:

~/Desktop/
~/Documents/
~/Downloads/

This avoids hard-coding project-specific folders such as:

~/folder1
~/folder2

4. Important note about the memory runtime attribute

For local Cromwell execution with this Docker configuration, the WDL should not include unresolved runtime declarations such as:

memory: memory

or:

memory: "~{memory_gb} GB"

if these values are not explicitly declared and passed correctly in each task.

Otherwise, Cromwell may fail during initialization with an error such as:

Task FASTQC has an invalid runtime attribute memory = !! NOT FOUND !!

For the local reusable configuration above, the runtime memory attribute is intentionally omitted from the backend configuration. Resource control is mainly handled by limiting parallel jobs using:

```bash
concurrent-job-limit = 2
```

5. Example command

Run the workflow using the Cromwell configuration file with -Dconfig.file before -jar:

```bash
java \
  -Dconfig.file=~/cromwell.local.reusable.shard2.no_memory_attr.conf \
  -jar ~/cromwell-92.jar \
  run ~/rMAP_Candida.no_runtime_memory.wdl \
  --inputs ~/rMAP-Candida.inputs_samples.json \
  --options ~/cromwell.options.no_docker_hash_lookup.json
```
---

## Troubleshooting

### Compleasm values remain `NA`

Check the Compleasm status & note columns in the report & summary TSV. Then inspect the task logs:

```bash
find cromwell-executions/rMAP_Candida \
  -path "*/call-COMPLEASM*/*/execution/*" \
  -type f | sort
```

Common causes include:

- Compleasm lineage download failed
- No internet access inside Docker
- Local lineage cache unavailable
- Assembly invalid or empty
- Memory limit too low
- Docker image issue

If running without internet, pre-download or mount the required Compleasm lineage resources if supported by your local container setup.

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

### `miniprot is required and not installed`

Use the fixed AMR container:

```text
gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed
```

Do not use the older AMR container for the current workflow unless it has been rebuilt with `miniprot` & the required ChroQueTas dependencies.

### Kraken2/Bracken database errors

Confirm the expected bundled database path inside the container:

```text
/opt/kraken2_db/candida
```

Expected files inside the Kraken2 database include:

```text
hash.k2d
opts.k2d
taxo.k2d
```

### Paired FASTQ input errors

Confirm that each sample has a matched R1 & R2:

```bash
ls -lh *_1.fastq.gz *_2.fastq.gz
```

For each sample:

```text
SAMPLE_1.fastq.gz  -> read1s
SAMPLE_2.fastq.gz  -> read2s
```

### Tree does not appear in the report

Check:

```bash
grep -n "do_phylogeny" rMAP-Candida.inputs*.json
grep -n "render_phylogeny_tree" rMAP-Candida.inputs*.json
```

Confirm that enough same-species samples were included:

```json
"rMAP_Candida.min_species_samples_for_tree": 3
```

For a 4-sample *Candida albicans* test run, all four samples should be eligible if species typing confirms *Candida albicans* and alignment/QC thresholds pass.

---

## Recommended quality thresholds

These are practical starting thresholds & should be adapted to the study design, sequencing depth, species & laboratory context.

| Metric | Threshold |
|---|---:|
| Species dominance by Kraken2/Bracken | High confidence: `>=95%`; review: `<95%` |
| Assembly size | Species-dependent; inspect for major deviations |
| Contig count | Lower is generally better |
| N50 | Higher is generally better |
| Compleasm completeness | Preferably `>90%`; `>95%` is stronger |
| AMR marker detection | Treat as screening evidence, not a phenotype |
| Core-site fraction | Default: `0.95` |
| Minimum species samples for tree | Default: `3` |

---

## Repository structure

Suggested repository layout:

```text
rMAP-Candida/
├── README.md
├── LICENSE
├── workflows/
│   └── rMAP_Candida.wdl
├── inputs/
│   ├── rMAP-Candida.inputs.example.json
│  
├── metadata/
│   └── surveillance_metadata.tsv
├── docs/
│   ├── index.html
│   └── reports/
│       └── example_report.html
├── docker/
│   ├── kraken2_bracken/
│   ├── amr_chroquetas/
│   └── candida_refs/
├── test_data/
    └── README.md

```

---

## Limitations

- The pipeline is not a standalone clinical diagnostic system.
- Genomic AMR markers may not fully predict phenotypic antifungal susceptibility.
- AMR interpretation depends on the completeness and accuracy of the installed marker database.
- Copy-number variation, aneuploidy, loss of heterozygosity, promoter changes & regulatory mechanisms may not be fully captured.
- Species identification depends on the completeness & quality of the Kraken2/Bracken database.
- Mixed cultures, contamination or low sequencing depth can affect interpretation.
- Compleasm completeness depends on assembly quality & lineage selection.
- Phylogenetic clustering & low SNP distances should not be interpreted as direct transmission without epidemiological metadata.
- Public SRA/ENA/BioSample metadata may be incomplete or missing; verify sample origin & collection dates before final surveillance interpretation.
- Results should be validated before use in clinical decision-making.

---

## Suggested citation

If you use this workflow, cite the repository & the major tools used by the workflow.

```text
rMAP-Candida: Rapid Mycological Analysis Pipeline for Candida genomic surveillance.
```

Also cite the relevant tools used in the workflow, including Cromwell/WDL, fastp, Kraken2, Bracken, MEGAHIT, QUAST, Compleasm, ChroQueTas, FungAMR, IQ-TREE, ETE3 & relevant AMR/database resources where applicable.

---

## Disclaimer

`rMAP-Candida` is intended for research, training & public health genomic surveillance. It is not intended for direct clinical decision-making unless the workflow, databases, reporting logic & interpretation rules have been validated under appropriate clinical laboratory, regulatory, biosafety & quality-management frameworks.
