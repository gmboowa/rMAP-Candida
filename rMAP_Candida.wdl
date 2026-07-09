version 1.0


workflow rMAP_Candida {
  input {
    Array[String]+ sample_names
    Array[File]+ read1s
    Array[File]+ read2s

    Boolean do_trimming = true
    Boolean do_quality_control = true
    Boolean do_species_typing = true
    Boolean do_assembly = true
    Boolean do_assembly_qc = true
    Boolean do_compleasm = false
    Boolean do_busco = false  # legacy/deprecated; retained only for backward-compatible JSONs
    Boolean do_fungal_amr = true
    Boolean do_phylogeny = false
    Boolean render_phylogeny_tree = true

    # Optional surveillance metadata TSV for interpretation in the HTML report.
    # Recommended columns: sample_id, country, site, collection_date, specimen_type, patient_group, ward_or_facility, sequencing_platform.
    File? surveillance_metadata_tsv

    # Dockerized Candida reference bundle.
    # Build/push this image separately and keep the references.tsv manifest inside it.
    # This lets users run phylogeny without providing local reference FASTA paths in JSON.
    String candida_refs_docker = "gmboowa/rmap-candida-refs:2026.05"
    String candida_refs_manifest = "/opt/rmap_candida_refs/references.tsv"
    String tree_visualization_docker = "gmboowa/ete3-render:1.18"

    # Species-aware phylogeny inputs.
    # Provide one reference per species that you want to include in phylogeny.
    # Example:
    #   phylogeny_reference_species = ["Candidozyma auris", "Candida albicans"]
    #   phylogeny_reference_fastas  = ["/path/C_auris_B8441.fasta", "/path/C_albicans_SC5314.fasta"]
    Array[String] phylogeny_reference_species = []
    Array[File] phylogeny_reference_fastas = []
    String snippy_docker = "staphb/snippy:4.6.0"
    String iqtree2_docker = "gmboowa/iqtree2-python:2.3.4"
    String iqtree2_model = "GTR+G"
    Int iqtree2_bootstraps = 1000
    Int min_species_samples_for_tree = 3
    Int min_mapping_quality = 20
    Int min_base_quality_for_phylogeny = 20
    Int min_depth_for_phylogeny = 10
    Int min_variant_quality_for_phylogeny = 20
    Float core_site_min_fraction = 0.95

    # Species-aware phylogeny branching.
    # rc170 borrows the working rMAP-TB logic: eligible species groups use
    # Snippy -> snippy-core -> IQ-TREE2 whenever listed here. This avoids the
    # fragile diploid bcftools consensus branch that caused Candida albicans
    # groups to be reported as variant_calling_failed before any tree could be built.
    # Keep the array editable in JSON if you want to add/remove species later.
    Array[String] snippy_phylogeny_species = ["Candidozyma auris", "Candida albicans"]
    Array[String] haploid_phylogeny_species = ["Nakaseomyces glabratus"]

    String fungal_kraken2_bracken_docker = "gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db"
    String fungamr_docker = "gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9"
    String megahit_docker = "quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6"
    String quast_docker = "staphb/quast:5.2.0"
    String compleasm_docker = "huangnengcsu/compleasm:v0.2.7"
    String busco_docker = "ezlabgva/busco:v5.7.1_cv1"  # legacy/deprecated
    String fastp_docker = "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    String fastqc_docker = "staphb/fastqc:0.12.1"
    String multiqc_docker = "multiqc/multiqc:v1.24"

    Int max_cpus = 8
    Int max_memory_gb = 32
    Int max_disk_gb = 1000
    Int min_read_length = 50
    Int bracken_read_length = 150
    String compleasm_lineage = "saccharomycetes"
    String compleasm_odb = "odb12"
    String busco_lineage = "saccharomycetes_odb10"  # legacy/deprecated
    String kraken_db_path = "/opt/kraken2_db/candida"
    String bracken_level = "S"
  }

  Int n = length(sample_names)
  Int cpu_4 = if max_cpus < 4 then max_cpus else 4
  Int cpu_8 = if max_cpus < 8 then max_cpus else 8

  # rc170 safety: even if an older JSON explicitly sets
  # rMAP_Candida.snippy_phylogeny_species = ["Candidozyma auris"], force
  # Candida albicans into the Snippy/core branch so the current 4-sample
  # albicans group can produce a species-specific tree.
  Array[String] effective_snippy_phylogeny_species = flatten([snippy_phylogeny_species, ["Candida albicans"]])

  #
  # Stage 1: per-sample read-level work and assembly.
  #
  # The first scatter deliberately contains only tasks that should run directly from
  # the paired FASTQ inputs.  Assembly-dependent tasks are launched in a second,
  # explicit scatter below over the complete Array[File] of assembly outputs.  This
  # prevents QUAST and AMR from accidentally receiving only shard-0 or a
  # single selected assembly FASTA.
  #
  scatter (i in range(n)) {
    if (do_trimming) {
      call FASTP_TRIMMING as TRIM {
        input:
          sample_name = sample_names[i],
          read1 = read1s[i],
          read2 = read2s[i],
          min_read_length = min_read_length,
          docker_image = fastp_docker,
          cpu = max_cpus,
          memory_gb = max_memory_gb,
          disk_gb = max_disk_gb
      }
    }

    File analysis_read1 = select_first([TRIM.trimmed_read1, read1s[i]])
    File analysis_read2 = select_first([TRIM.trimmed_read2, read2s[i]])

    if (do_quality_control) {
      call FASTQC as QC {
        input:
          sample_name = sample_names[i],
          read1 = analysis_read1,
          read2 = analysis_read2,
          docker_image = fastqc_docker,
          cpu = max_cpus,
          memory_gb = max_memory_gb,
          disk_gb = max_disk_gb
      }
    }

    if (do_species_typing) {
      call FUNGAL_SPECIES_TYPING as SPECIES {
        input:
          sample_name = sample_names[i],
          read1 = analysis_read1,
          read2 = analysis_read2,
          kraken_db_path = kraken_db_path,
          bracken_read_length = bracken_read_length,
          bracken_level = bracken_level,
          docker_image = fungal_kraken2_bracken_docker,
          cpu = max_cpus,
          memory_gb = max_memory_gb,
          disk_gb = max_disk_gb
      }
    }

    if (do_assembly) {
      call FUNGAL_ASSEMBLY as ASSEMBLY {
        input:
          sample_name = sample_names[i],
          read1 = analysis_read1,
          read2 = analysis_read2,
          docker_image = megahit_docker,
          cpu = max_cpus,
          memory_gb = max_memory_gb,
          disk_gb = max_disk_gb
      }
    }
  }

  #
  # Stage 2: per-assembly downstream work.
  #
  # This second scatter is intentionally outside the assembly scatter.  Cromwell
  # must first collect all ASSEMBLY.contigs_fasta outputs, then QUAST, optional BUSCO/Compleasm, and
  # AMR are scattered over the full array.  With two samples, this creates:
  #   call-QUAST/shard-0 and call-QUAST/shard-1
  #   call-BUSCO/shard-0     and call-BUSCO/shard-1       only if do_busco=true
  #   call-COMPLEASM/shard-0 and call-COMPLEASM/shard-1   only if do_compleasm=true
  #   call-AMR/shard-0       and call-AMR/shard-1
  #
  if (do_assembly) {
    Array[File] assembled_contigs = select_all(ASSEMBLY.contigs_fasta)

    if (do_assembly_qc) {
      scatter (q in range(length(assembled_contigs))) {
        call ASSEMBLY_QC as QUAST {
          input:
            sample_name = sample_names[q],
            contigs = assembled_contigs[q],
            docker_image = quast_docker,
            cpu = max_cpus,
            memory_gb = max_memory_gb,
            disk_gb = max_disk_gb
        }
      }
    }

    if (do_busco) {
      scatter (b in range(length(assembled_contigs))) {
        call BUSCO_FUNGAL as BUSCO {
          input:
            sample_name = sample_names[b],
            assembly_fasta = assembled_contigs[b],
            busco_lineage = busco_lineage,
            busco_docker = busco_docker,
            cpu = max_cpus,
            memory_gb = max_memory_gb,
            disk_gb = max_disk_gb
        }
      }
    }

    if (do_compleasm) {
      scatter (b in range(length(assembled_contigs))) {
        call COMPLEASM_FUNGAL as COMPLEASM {
          input:
            sample_name = sample_names[b],
            assembly_fasta = assembled_contigs[b],
            compleasm_lineage = compleasm_lineage,
            compleasm_odb = compleasm_odb,
            compleasm_docker = compleasm_docker,
            cpu = max_cpus,
            memory_gb = max_memory_gb,
            disk_gb = max_disk_gb
        }
      }
    }

    if (do_fungal_amr && do_species_typing) {
      Array[File] species_summary_files_for_amr = select_all(SPECIES.top_species_tsv)

      scatter (a in range(length(assembled_contigs))) {
        call FUNGAL_AMR_CHARACTERIZATION as AMR {
          input:
            sample_id = sample_names[a],
            assembly_fasta = assembled_contigs[a],
            species_summary_tsv = species_summary_files_for_amr[a],
            fungal_amr_docker_image = fungamr_docker,
            threads = max_cpus,
            memory_gb = max_memory_gb,
            disk_gb = max_disk_gb
        }
      }
    }
  }


  #
  # Stage 3: optional species-aware Candida core-SNP phylogeny.
  #
  # IMPORTANT:
  # Phylogeny is performed separately per top species call. Do not mix Candida species
  # in one SNP tree. A species-specific reference FASTA must be supplied/exported for
  # each species to be included. rc170 uses the rMAP-TB style Snippy/snippy-core
  # branch for species listed in snippy_phylogeny_species, now including Candida albicans
  # by default. Other species still have the bcftools consensus fallback.
  # rc172 patch: invalid/unclassified groups such as NA, Unknown, Unclassified,
  # and No_call are skipped before phylogeny; valid species groups continue.
  # Recombinant regions are not explicitly filtered in this implementation; resulting
  # trees should be interpreted as broad genomic relatedness/lineage structure rather than
  # definitive transmission inference.
  #
  if (do_phylogeny && do_species_typing) {
    call CANDIDA_EXPORT_REFERENCE_FASTAS {
      input:
        refs_manifest = candida_refs_manifest,
        docker_image = candida_refs_docker,
        cpu = max_cpus,
        memory_gb = max_memory_gb,
        disk_gb = max_disk_gb
    }

    call CANDIDA_SNIPPY_CORE_BY_SPECIES {
      input:
        sample_names = sample_names,
        read1s = analysis_read1,
        read2s = analysis_read2,
        species_top_tsvs = select_all(SPECIES.top_species_tsv),
        reference_species = CANDIDA_EXPORT_REFERENCE_FASTAS.reference_species,
        reference_fastas = CANDIDA_EXPORT_REFERENCE_FASTAS.reference_fastas,
        min_species_samples_for_tree = min_species_samples_for_tree,
        docker_image = snippy_docker,
        cpu = max_cpus,
        memory_gb = max_memory_gb,
        disk_gb = max_disk_gb,
        min_quality = min_mapping_quality,
        min_base_quality = min_base_quality_for_phylogeny,
        min_depth = min_depth_for_phylogeny,
        min_variant_quality = min_variant_quality_for_phylogeny,
        core_site_min_fraction = core_site_min_fraction,
        haploid_species = haploid_phylogeny_species,
        snippy_species = effective_snippy_phylogeny_species
    }

    Array[String] phylogeny_group_labels = read_lines(CANDIDA_SNIPPY_CORE_BY_SPECIES.group_labels_txt)

    scatter (pt in range(length(CANDIDA_SNIPPY_CORE_BY_SPECIES.core_full_alignments))) {
      call CANDIDA_IQTREE2_PHYLOGENY as CANDIDA_IQTREE {
        input:
          alignment = CANDIDA_SNIPPY_CORE_BY_SPECIES.core_full_alignments[pt],
          species_label = phylogeny_group_labels[pt],
          model = iqtree2_model,
          bootstrap_replicates = iqtree2_bootstraps,
          docker_image = iqtree2_docker,
          cpu = max_cpus,
          memory_gb = max_memory_gb,
          disk_gb = max_disk_gb
      }

      if (render_phylogeny_tree) {
        call CANDIDA_TREE_VISUALIZATION as CANDIDA_TREE {
          input:
            input_tree = CANDIDA_IQTREE.final_tree,
            species_label = phylogeny_group_labels[pt],
            docker_image = tree_visualization_docker,
            width = 2600,
            height = 1800,
            image_format = "png",
            cpu = max_cpus,
            memory_gb = max_memory_gb,
            disk_gb = max_disk_gb
        }
      }
    }
  }

  if (do_quality_control) {
    call MULTIQC_REPORT {
      input:
        fastqc_reports = flatten(select_all(QC.fastqc_reports)),
        trimming_json = select_all(TRIM.fastp_json),
        trimming_html = select_all(TRIM.fastp_html),
        docker_image = multiqc_docker,
        cpu = max_cpus,
        memory_gb = max_memory_gb,
        disk_gb = max_disk_gb
    }
  }

  call MERGE_MYC_REPORTS {
    input:
      sample_names = sample_names,
      multiqc_report = MULTIQC_REPORT.multiqc_report,
      trimming_htmls = select_all(TRIM.fastp_html),
      species_top_tsvs = select_all(SPECIES.top_species_tsv),
      kraken_reports = select_all(SPECIES.kraken_report),
      bracken_reports = select_all(SPECIES.bracken_report),
      assembly_summaries = select_all(ASSEMBLY.assembly_summary_tsv),
      quast_reports = select_first([QUAST.quast_report_tsv, []]),
      busco_summaries = BUSCO.busco_short_summary_txt,
      busco_summary_tsvs = BUSCO.busco_summary_tsv,
      compleasm_summaries = COMPLEASM.busco_short_summary_txt,
      compleasm_summary_tsvs = COMPLEASM.busco_summary_tsv,
      amr_summaries = select_first([AMR.amr_summary_tsv, []]),
      amr_htmls = select_first([AMR.amr_report_html, []]),
      phylogeny_group_summary = CANDIDA_SNIPPY_CORE_BY_SPECIES.phylogeny_group_summary_tsv,
      phylogeny_core_alignments = select_first([CANDIDA_SNIPPY_CORE_BY_SPECIES.core_full_alignments, []]),
      phylogeny_newick_trees = select_first([CANDIDA_IQTREE.final_tree, []]),
      phylogeny_iqtree_reports = select_first([CANDIDA_IQTREE.iqtree_report, []]),
      phylogeny_tree_images = select_all(select_first([CANDIDA_TREE.tree_image, []])),
      surveillance_metadata_tsv = surveillance_metadata_tsv,
      cpu = max_cpus,
      memory_gb = max_memory_gb,
      disk_gb = max_disk_gb
  }

  output {
    File rmap_myc_html_report = MERGE_MYC_REPORTS.html_report
    File rmap_myc_summary_tsv = MERGE_MYC_REPORTS.summary_tsv
    File rmap_candida_surveillance_summary_tsv = MERGE_MYC_REPORTS.surveillance_summary_tsv
    File rmap_candida_pairwise_snp_distances_tsv = MERGE_MYC_REPORTS.pairwise_snp_distances_tsv
    Array[File] trimmed_reads_1 = select_all(TRIM.trimmed_read1)
    Array[File] trimmed_reads_2 = select_all(TRIM.trimmed_read2)
    File? multiqc_report = MULTIQC_REPORT.multiqc_report
    Array[File] fungal_species_summaries = select_all(SPECIES.top_species_tsv)
    Array[File] fungal_assemblies = select_all(ASSEMBLY.contigs_fasta)
    Array[File] fungal_amr_summaries = select_first([AMR.amr_summary_tsv, []])
    File? candida_phylogeny_group_summary = CANDIDA_SNIPPY_CORE_BY_SPECIES.phylogeny_group_summary_tsv
    Array[File] candida_core_snp_alignments = select_first([CANDIDA_SNIPPY_CORE_BY_SPECIES.core_full_alignments, []])
    Array[File] candida_phylogeny_newick_trees = select_first([CANDIDA_IQTREE.final_tree, []])
    Array[File] candida_iqtree_reports = select_first([CANDIDA_IQTREE.iqtree_report, []])
    Array[File] candida_phylogeny_tree_images = select_all(select_first([CANDIDA_TREE.tree_image, []]))
  }
}

task FASTP_TRIMMING {
  input {
    String sample_name
    File read1
    File read2
    Int min_read_length
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    fastp \
      --in1 ~{read1} \
      --in2 ~{read2} \
      --out1 ~{sample_name}_R1.trimmed.fastq.gz \
      --out2 ~{sample_name}_R2.trimmed.fastq.gz \
      --length_required ~{min_read_length} \
      --detect_adapter_for_pe \
      --thread ~{cpu} \
      --html ~{sample_name}.fastp.html \
      --json ~{sample_name}.fastp.json
  >>>

  output {
    File trimmed_read1 = "~{sample_name}_R1.trimmed.fastq.gz"
    File trimmed_read2 = "~{sample_name}_R2.trimmed.fastq.gz"
    File fastp_html = "~{sample_name}.fastp.html"
    File fastp_json = "~{sample_name}.fastp.json"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task FASTQC {
  input {
    String sample_name
    File read1
    File read2
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    mkdir -p fastqc_out
    fastqc -t ~{cpu} -o fastqc_out ~{read1} ~{read2}
  >>>

  output {
    Array[File] fastqc_reports = glob("fastqc_out/*_fastqc.html")
    Array[File] fastqc_zips = glob("fastqc_out/*_fastqc.zip")
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task MULTIQC_REPORT {
  input {
    Array[File] fastqc_reports
    Array[File] trimming_json
    Array[File] trimming_html
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    mkdir -p multiqc_inputs
    cp ~{sep=' ' fastqc_reports} multiqc_inputs/ || true
    cp ~{sep=' ' trimming_json} multiqc_inputs/ || true
    multiqc multiqc_inputs -o . -n rMAP-Myc-Candida-Candida_MultiQC.html --force
  >>>

  output {
    File multiqc_report = "rMAP-Myc-Candida-Candida_MultiQC.html"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task FUNGAL_SPECIES_TYPING {
  input {
    String sample_name
    File read1
    File read2
    String kraken_db_path
    Int bracken_read_length
    String bracken_level
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    echo "=============================================="
    echo "rMAP-Myc-Candida-Candida species typing"
    echo "Sample: ~{sample_name}"
    echo "Kraken2 DB: ~{kraken_db_path}"
    echo "=============================================="

    if [ ! -d "~{kraken_db_path}" ]; then
      echo "ERROR: Kraken2 database directory not found: ~{kraken_db_path}" >&2
      exit 1
    fi

    for f in hash.k2d opts.k2d taxo.k2d; do
      if [ ! -s "~{kraken_db_path}/$f" ]; then
        echo "ERROR: Required Kraken2 DB file missing or empty: ~{kraken_db_path}/$f" >&2
        echo "Database contents:" >&2
        ls -lh "~{kraken_db_path}" >&2 || true
        exit 1
      fi
    done

    echo "Database contents:"
    ls -lh "~{kraken_db_path}"

    echo "Running Kraken2..."
    kraken2 \
      --db "~{kraken_db_path}" \
      --threads ~{cpu} \
      --paired \
      --gzip-compressed \
      --memory-mapping \
      --use-names \
      --report "~{sample_name}.kraken2.report" \
      --output "~{sample_name}.kraken2.output" \
      "~{read1}" "~{read2}"

    echo "Running Bracken..."
    if command -v bracken >/dev/null 2>&1; then
      set +e
      bracken \
        -d "~{kraken_db_path}" \
        -i "~{sample_name}.kraken2.report" \
        -o "~{sample_name}.bracken.tsv" \
        -r ~{bracken_read_length} \
        -l "~{bracken_level}"
      BRACKEN_RC=$?
      set -e

      if [ "${BRACKEN_RC}" != "0" ] || [ ! -s "~{sample_name}.bracken.tsv" ]; then
        echo "WARNING: Bracken failed or produced no output. Creating placeholder Bracken TSV." >&2
        printf "name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
" > "~{sample_name}.bracken.tsv"
      fi
    else
      echo "WARNING: bracken command not found. Creating placeholder Bracken TSV." >&2
      printf "name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
" > "~{sample_name}.bracken.tsv"
    fi

    echo "Creating top species summary..."
    python3 <<'PY'
import csv
from pathlib import Path

sample = "~{sample_name}"
bracken_path = Path(f"{sample}.bracken.tsv")
kraken_path = Path(f"{sample}.kraken2.report")
out_path = Path(f"{sample}.top_species.tsv")

def clean_name(x):
    return " ".join(str(x).strip().split())

top = None

# Prefer Bracken species-level abundance if available.
if bracken_path.exists() and bracken_path.stat().st_size > 0:
    with bracken_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="	")
        rows = []
        for row in reader:
            name = row.get("name") or row.get("taxonomy_name") or row.get("scientific_name") or ""
            taxid = row.get("taxonomy_id") or row.get("taxid") or ""
            reads = row.get("new_est_reads") or row.get("fraction_total_reads") or row.get("kraken_assigned_reads") or "0"
            fraction = row.get("fraction_total_reads") or "0"
            try:
                score = float(row.get("new_est_reads") or 0)
            except Exception:
                score = 0.0
            if name.strip():
                rows.append((score, clean_name(name), taxid, reads, fraction))
        if rows:
            rows.sort(reverse=True, key=lambda x: x[0])
            score, name, taxid, reads, fraction = rows[0]
            top = {
                "sample_id": sample,
                "top_species": name,
                "taxid": taxid,
                "estimated_reads": reads,
                "fraction_total_reads": fraction,
                "evidence": "Bracken species-level abundance"
            }

# Fall back to Kraken2 report if Bracken has no species rows.
if top is None and kraken_path.exists() and kraken_path.stat().st_size > 0:
    best = None
    with kraken_path.open() as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            try:
                pct = float(parts[0].strip())
                clade_reads = int(parts[1].strip())
            except Exception:
                continue
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name = clean_name(parts[5])
            if rank == "S":
                candidate = (clade_reads, pct, name, taxid)
                if best is None or candidate[0] > best[0]:
                    best = candidate
    if best is not None:
        clade_reads, pct, name, taxid = best
        top = {
            "sample_id": sample,
            "top_species": name,
            "taxid": taxid,
            "estimated_reads": str(clade_reads),
            "fraction_total_reads": str(pct / 100.0),
            "evidence": "Kraken2 species-level report"
        }

if top is None:
    top = {
        "sample_id": sample,
        "top_species": "No species-level Candida call",
        "taxid": "NA",
        "estimated_reads": "0",
        "fraction_total_reads": "0",
        "evidence": "No species-level call in Bracken or Kraken2 report"
    }

with out_path.open("w", newline="") as out:
    fields = ["sample_id", "top_species", "taxid", "estimated_reads", "fraction_total_reads", "evidence"]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="	")
    writer.writeheader()
    writer.writerow(top)

print(out_path.read_text())
PY

    echo "Species typing completed for ~{sample_name}"
  >>>

  output {
    File kraken_report = "~{sample_name}.kraken2.report"
    File kraken_output = "~{sample_name}.kraken2.output"
    File bracken_report = "~{sample_name}.bracken.tsv"
    File top_species_tsv = "~{sample_name}.top_species.tsv"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


task FUNGAL_ASSEMBLY {
  input {
    String sample_name
    File read1
    File read2
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    echo "=============================================="
    echo "rMAP-Myc-Candida assembly using MEGAHIT"
    echo "Sample: ~{sample_name}"
    echo "Read 1: ~{read1}"
    echo "Read 2: ~{read2}"
    echo "Threads: ~{cpu}"
    echo "=============================================="

    rm -rf megahit_~{sample_name}

    megahit \
      -1 ~{read1} \
      -2 ~{read2} \
      -o megahit_~{sample_name} \
      -t ~{cpu} \
      --min-contig-len 500

    if [ ! -s megahit_~{sample_name}/final.contigs.fa ]; then
      echo "ERROR: MEGAHIT did not produce final.contigs.fa for ~{sample_name}" >&2
      exit 1
    fi

    cp megahit_~{sample_name}/final.contigs.fa ~{sample_name}.contigs.fasta

    # Compute assembly summary without requiring python inside the MEGAHIT image.
    awk '
      /^>/ {
        if (seqlen > 0) { print seqlen }
        seqlen = 0
        next
      }
      {
        gsub(/[[:space:]]/, "", $0)
        seqlen += length($0)
      }
      END {
        if (seqlen > 0) { print seqlen }
      }
    ' ~{sample_name}.contigs.fasta > contig_lengths.txt

    if [ ! -s contig_lengths.txt ]; then
      printf "sample_id\tassembler\tcontigs\ttotal_bp\tn50\tlargest_contig\n" > ~{sample_name}.assembly_summary.tsv
      printf "~{sample_name}\tMEGAHIT\t0\t0\t0\t0\n" >> ~{sample_name}.assembly_summary.tsv
    else
      CONTIGS=$(wc -l < contig_lengths.txt | tr -d ' ')
      TOTAL_BP=$(awk '{s += $1} END {print s + 0}' contig_lengths.txt)

      # Avoid pipefail/SIGPIPE failures from patterns such as: sort | head
      # This keeps successful MEGAHIT assemblies from being marked failed with rc=141.
      sort -nr contig_lengths.txt > contig_lengths.sorted.txt

      LARGEST_CONTIG=$(awk 'NR == 1 {print $1}' contig_lengths.sorted.txt)
      HALF_BP=$(( (TOTAL_BP + 1) / 2 ))
      N50=$(awk -v half="${HALF_BP}" '
        BEGIN {s = 0; n50 = 0}
        {
          s += $1
          if (n50 == 0 && s >= half) {
            n50 = $1
          }
        }
        END {print n50 + 0}
      ' contig_lengths.sorted.txt)

      printf "sample_id\tassembler\tcontigs\ttotal_bp\tn50\tlargest_contig\n" > ~{sample_name}.assembly_summary.tsv
      printf "~{sample_name}\tMEGAHIT\t%s\t%s\t%s\t%s\n" "${CONTIGS}" "${TOTAL_BP}" "${N50:-0}" "${LARGEST_CONTIG:-0}" >> ~{sample_name}.assembly_summary.tsv
    fi

    echo "Assembly completed for ~{sample_name}"
    cat ~{sample_name}.assembly_summary.tsv
  >>>

  output {
    File contigs_fasta = "~{sample_name}.contigs.fasta"
    File assembly_summary_tsv = "~{sample_name}.assembly_summary.tsv"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    continueOnReturnCode: [0, 141]
  }
}

task ASSEMBLY_QC {
  input {
    String sample_name
    File contigs
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    RAW_DIR="quast_~{sample_name}"
    STD_REPORT="~{sample_name}.quast.report.tsv"

    quast.py ~{contigs} -o "$RAW_DIR" -t ~{cpu} --min-contig 500 || true

    # Build a standardized one-row QUAST TSV without requiring python3 in the QUAST image.
    # Some staphb/quast images include quast.py but not a python3 executable in PATH.
    RAW_REPORT="$RAW_DIR/report.tsv"
    printf "sample_id\t# contigs\tLargest contig\tTotal length\tGC (%%)\tN50\tN90\tL50\tL90\t# N's per 100 kbp\tquast_status\n" > "$STD_REPORT"

    get_metric() {
      metric="$1"
      if [ -s "$RAW_REPORT" ]; then
        awk -F '\t' -v m="$metric" '$1 == m {print $2; found=1; exit} END {if (!found) print "NA"}' "$RAW_REPORT"
      else
        printf "NA"
      fi
    }

    if [ -s "$RAW_REPORT" ]; then
      STATUS="QUAST parsed"
    else
      STATUS="QUAST report unavailable"
    fi

    CONTIGS=$(get_metric "# contigs")
    if [ "$CONTIGS" = "NA" ]; then CONTIGS=$(get_metric "# contigs (>= 0 bp)"); fi
    if [ "$CONTIGS" = "NA" ]; then CONTIGS=$(get_metric "# contigs (>= 500 bp)"); fi

    LARGEST=$(get_metric "Largest contig")

    TOTAL=$(get_metric "Total length")
    if [ "$TOTAL" = "NA" ]; then TOTAL=$(get_metric "Total length (>= 0 bp)"); fi
    if [ "$TOTAL" = "NA" ]; then TOTAL=$(get_metric "Total length (>= 500 bp)"); fi

    GC=$(get_metric "GC (%)")
    N50=$(get_metric "N50")
    N90=$(get_metric "N90")
    L50=$(get_metric "L50")
    L90=$(get_metric "L90")
    NS=$(get_metric "# N's per 100 kbp")

    printf "~{sample_name}\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$CONTIGS" "$LARGEST" "$TOTAL" "$GC" "$N50" "$N90" "$L50" "$L90" "$NS" "$STATUS" >> "$STD_REPORT"

    test -s "$STD_REPORT"
  >>>

  output {
    File quast_report_tsv = "~{sample_name}.quast.report.tsv"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task BUSCO_FUNGAL {
  input {
    String sample_name
    File assembly_fasta
    String busco_docker
    String busco_lineage
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set +e
    set +u
    set +o pipefail

    SAMPLE="~{sample_name}"
    ASM="~{assembly_fasta}"
    PRIMARY_LINEAGE="~{busco_lineage}"
    SUMMARY_TSV="${SAMPLE}.busco_summary.tsv"
    SHORT_TXT="${SAMPLE}.busco.short_summary.txt"
    LOG="${SAMPLE}.busco.log"
    DETAILS="${SAMPLE}.busco.details.txt"

    HEADER="sample_id	complete_busco_pct	single_copy_busco_pct	duplicated_busco_pct	fragmented_busco_pct	missing_busco_pct	busco_markers	busco_lineage	busco_status	busco_note"
    printf "%b\n" "${HEADER}" > "${SUMMARY_TSV}"
    : > "${LOG}"
    : > "${DETAILS}"

    echo "[$(date)] BUSCO task started for ${SAMPLE}" >> "${LOG}"
    echo "Assembly FASTA: ${ASM}" >> "${LOG}"
    echo "Primary lineage: ${PRIMARY_LINEAGE}" >> "${LOG}"
    echo "CPU: ~{cpu}; memory_gb: ~{memory_gb}" >> "${LOG}"
    command -v busco >> "${LOG}" 2>&1 || true
    busco --version >> "${LOG}" 2>&1 || true
    python3 --version >> "${LOG}" 2>&1 || true

    # Candida assemblies often parse best against saccharomycetes_odb10, but some
    # BUSCO containers/environments only have ascomycota/fungi datasets cached.
    # This task therefore tries the requested lineage first, then the Candida-safe
    # fallbacks. A BUSCO summary is accepted if it exists, even when BUSCO exits
    # non-zero after writing results; this avoids losing valid BUSCO metrics.
    LINEAGES="${PRIMARY_LINEAGE} saccharomycetes_odb10 ascomycota_odb10 fungi_odb10"
    BUSCO_OK=0
    USED_LINEAGES=""
    SUMMARY_FILE=""
    FINAL_LINEAGE="${PRIMARY_LINEAGE}"
    FINAL_RC=127

    run_and_find_busco_summary () {
      L="$1"
      OUTDIR="$2"
      MODE_TAG="$3"
      echo "[$(date)] Trying BUSCO lineage=${L}; mode=${MODE_TAG}" >> "${LOG}"
      EXTRA_ARGS=""
      if [ "${MODE_TAG}" = "offline" ]; then
        EXTRA_ARGS="--offline"
      fi
      busco \
        -i "${ASM}" \
        -o "run_${SAMPLE}_${L}_${MODE_TAG}" \
        -m genome \
        -l "${L}" \
        --cpu "~{cpu}" \
        --out_path "${OUTDIR}" \
        --download_path busco_downloads \
        --force ${EXTRA_ARGS} >> "${LOG}" 2>&1
      RC=$?
      FOUND="$(find "${OUTDIR}" -type f \( -name 'short_summary*.txt' -o -name 'short_summary*.json' \) 2>/dev/null | sort | head -n 1)"
      echo "[$(date)] BUSCO lineage=${L}; mode=${MODE_TAG}; rc=${RC}; summary=${FOUND}" >> "${LOG}"
      printf "%s\t%s\n" "${RC}" "${FOUND}"
    }

    for L in ${LINEAGES}; do
      [ -z "${L}" ] && continue
      if echo " ${USED_LINEAGES} " | grep -q " ${L} "; then
        continue
      fi
      USED_LINEAGES="${USED_LINEAGES} ${L}"

      OUTDIR="busco_${SAMPLE}_${L}"
      RESULT="$(run_and_find_busco_summary "${L}" "${OUTDIR}" "online")"
      FINAL_RC="$(printf '%s' "${RESULT}" | awk -F '\t' '{print $1}')"
      FOUND="$(printf '%s' "${RESULT}" | awk -F '\t' '{print $2}')"

      # If online download failed, try offline with any cached lineages inside the container/workdir.
      if [ -z "${FOUND}" ]; then
        OUTDIR_OFF="busco_${SAMPLE}_${L}_offline"
        RESULT="$(run_and_find_busco_summary "${L}" "${OUTDIR_OFF}" "offline")"
        FINAL_RC="$(printf '%s' "${RESULT}" | awk -F '\t' '{print $1}')"
        FOUND="$(printf '%s' "${RESULT}" | awk -F '\t' '{print $2}')"
      fi

      if [ -n "${FOUND}" ]; then
        BUSCO_OK=1
        SUMMARY_FILE="${FOUND}"
        FINAL_LINEAGE="${L}"
        break
      fi
    done

    # Last-resort auto-lineage attempt. This is slower but helps when the requested
    # lineage name is not available or when BUSCO can infer a better fungal lineage.
    if [ "${BUSCO_OK}" -ne 1 ]; then
      OUTDIR="busco_${SAMPLE}_auto_lineage_euk"
      echo "[$(date)] Trying BUSCO auto-lineage-euk fallback" >> "${LOG}"
      busco \
        -i "${ASM}" \
        -o "run_${SAMPLE}_auto_lineage_euk" \
        -m genome \
        --auto-lineage-euk \
        --cpu "~{cpu}" \
        --out_path "${OUTDIR}" \
        --download_path busco_downloads \
        --force >> "${LOG}" 2>&1
      FINAL_RC=$?
      FOUND="$(find "${OUTDIR}" -type f \( -name 'short_summary*.txt' -o -name 'short_summary*.json' \) 2>/dev/null | sort | head -n 1)"
      echo "[$(date)] BUSCO auto-lineage-euk; rc=${FINAL_RC}; summary=${FOUND}" >> "${LOG}"
      if [ -n "${FOUND}" ]; then
        BUSCO_OK=1
        SUMMARY_FILE="${FOUND}"
        FINAL_LINEAGE="auto-lineage-euk"
      fi
    fi

    if [ "${BUSCO_OK}" -eq 1 ] && [ -n "${SUMMARY_FILE}" ]; then
      cp "${SUMMARY_FILE}" "${SHORT_TXT}" 2>/dev/null || true
      python3 <<PYBUSCO
import json, pathlib, re, csv, traceback
sample = "${SAMPLE}"
summary_file = pathlib.Path("${SUMMARY_FILE}")
out = pathlib.Path("${SUMMARY_TSV}")
lineage = "${FINAL_LINEAGE}"
header = ["sample_id","complete_busco_pct","single_copy_busco_pct","duplicated_busco_pct","fragmented_busco_pct","missing_busco_pct","busco_markers","busco_lineage","busco_status","busco_note"]
vals = dict.fromkeys(header, "NA")
vals.update({"sample_id": sample, "busco_lineage": lineage, "busco_status": "PASS", "busco_note": "BUSCO completed successfully."})

def walk(obj):
    if isinstance(obj, dict):
        for k,v in obj.items():
            yield str(k), v
            yield from walk(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v)

def pick_json_metrics(d):
    flat = {k.lower().replace(" ","_").replace("%",""): v for k,v in walk(d)}
    # BUSCO JSON keys vary by version. Try common names first.
    mapping = {
        "complete_busco_pct": ["complete_percentage", "complete_busco_percentage", "c", "complete"],
        "single_copy_busco_pct": ["single_copy_percentage", "single_copy_busco_percentage", "s", "single_copy"],
        "duplicated_busco_pct": ["multi_copy_percentage", "duplicated_percentage", "duplicated_busco_percentage", "d", "duplicated"],
        "fragmented_busco_pct": ["fragmented_percentage", "fragmented_busco_percentage", "f", "fragmented"],
        "missing_busco_pct": ["missing_percentage", "missing_busco_percentage", "m", "missing"],
        "busco_markers": ["number_of_buscos", "n_markers", "busco_n", "n"]
    }
    hit = False
    for outk, keys in mapping.items():
        for kk in keys:
            if kk in flat and str(flat[kk]).strip() not in ["", "None", "NA"]:
                val = str(flat[kk]).strip()
                if outk.endswith("pct"):
                    m = re.search(r"[0-9.]+", val)
                    if m:
                        vals[outk] = f"{float(m.group(0)):.1f}"
                        hit = True
                        break
                else:
                    m = re.search(r"\d+", val)
                    if m:
                        vals[outk] = m.group(0)
                        hit = True
                        break
    return hit

try:
    txt = summary_file.read_text(errors="ignore")
    parsed_json = False
    if summary_file.suffix.lower() == ".json":
        try:
            parsed_json = pick_json_metrics(json.loads(txt))
        except Exception:
            parsed_json = False
    candidates = []
    candidates.extend([x.strip() for x in txt.splitlines() if x.strip()])
    candidates.append(txt.replace("\n", " "))
    # BUSCO v5 short summary pattern, allowing variable spaces and optional percent signs.
    pattern = re.compile(r"C:\s*([0-9.]+)%\s*\[S:\s*([0-9.]+)%,\s*D:\s*([0-9.]+)%\],\s*F:\s*([0-9.]+)%,\s*M:\s*([0-9.]+)%,\s*n:\s*(\d+)")
    match = None
    for c in candidates:
        match = pattern.search(c)
        if match:
            break
    if match:
        vals["complete_busco_pct"] = f"{float(match.group(1)):.1f}"
        vals["single_copy_busco_pct"] = f"{float(match.group(2)):.1f}"
        vals["duplicated_busco_pct"] = f"{float(match.group(3)):.1f}"
        vals["fragmented_busco_pct"] = f"{float(match.group(4)):.1f}"
        vals["missing_busco_pct"] = f"{float(match.group(5)):.1f}"
        vals["busco_markers"] = match.group(6)
    elif not parsed_json:
        vals["busco_status"] = "PARSE_WARNING"
        vals["busco_note"] = "BUSCO produced a summary file, but no recognized BUSCO metric pattern was found. Inspect the copied short_summary and busco log."
except Exception as e:
    vals["busco_status"] = "PARSE_FAILED"
    vals["busco_note"] = "BUSCO summary parsing failed: " + str(e).replace("\t", " ").replace("\n", " ")
with out.open("w", newline="") as w:
    writer = csv.DictWriter(w, fieldnames=header, delimiter="\t")
    writer.writeheader()
    writer.writerow(vals)
PYBUSCO
    else
      echo "No BUSCO short_summary file was produced after trying: ${USED_LINEAGES} and auto-lineage-euk" > "${SHORT_TXT}"
      tail -n 160 "${LOG}" >> "${SHORT_TXT}" || true
      NOTE="$(tail -n 50 "${LOG}" | tr '\t' ' ' | tr '\n' ' ' | sed 's/"/ /g')"
      printf "%b\n" "${SAMPLE}\tNA\tNA\tNA\tNA\tNA\tNA\t${PRIMARY_LINEAGE}\tBUSCO_FAILED\tNo BUSCO summary produced after lineage retries and auto-lineage-euk. Check lineage download/cache, container internet access, assembly validity, and BUSCO log. Last log: ${NOTE}" >> "${SUMMARY_TSV}"
    fi

    # Final guard: Cromwell/report must never receive a header-only BUSCO table.
    if [ "$(wc -l < "${SUMMARY_TSV}" | tr -d ' ')" -lt 2 ]; then
      NOTE="$(tail -n 50 "${LOG}" | tr '\t' ' ' | tr '\n' ' ' | sed 's/"/ /g')"
      printf "%b\n" "${SAMPLE}\tNA\tNA\tNA\tNA\tNA\tNA\t${FINAL_LINEAGE}\tBUSCO_PARSE_FAILED\tBUSCO ran but normalized TSV was empty/header-only. Last log: ${NOTE}" >> "${SUMMARY_TSV}"
    fi

    cp "${LOG}" "${DETAILS}" 2>/dev/null || true
    exit 0
  >>>

  output {
    File busco_short_summary_txt = "~{sample_name}.busco.short_summary.txt"
    File busco_summary_tsv = "~{sample_name}.busco_summary.tsv"
    File busco_log = "~{sample_name}.busco.log"
    File busco_details = "~{sample_name}.busco.details.txt"
  }

  runtime {
    docker: "~{busco_docker}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}



task FUNGAL_AMR_CHARACTERIZATION {
  input {
    String sample_id
    File assembly_fasta
    File species_summary_tsv
    String fungal_amr_docker_image
    Int threads = 4
    Int memory_gb = 32
    Int disk_gb = 1000
  }

  command <<<
    set +e
    set +u
    set +o pipefail

    SAMPLE="~{sample_id}"
    ASM="~{assembly_fasta}"
    SPECIES_TSV="~{species_summary_tsv}"
    THREADS="~{threads}"

    SUMMARY_OUT="${SAMPLE}.fungal_amr.summary.tsv"
    RAW_OUT="${SAMPLE}.fungal_amr.raw.tsv"
    HTML_OUT="${SAMPLE}.fungal_amr.html"
    LOG_OUT="${SAMPLE}.fungal_amr.log"

    mkdir -p amr_out

    # Cromwell output guard: create all declared outputs before any scanner call.
    printf "sample_id\tspecies\tdrug_class\tdrug\tgene_or_status\tmutation\teffect\tevidence_level\tinterpretation\n" > "${SUMMARY_OUT}"
    printf "%s\tNA\tazole/other\tfluconazole/other\tAMR task initialized\tNA\tNA\tAMR_TASK_INITIALIZED\tDeclared AMR outputs were initialized before scanner execution.\n" "${SAMPLE}" >> "${SUMMARY_OUT}"
    printf "source_file\traw_text\n" > "${RAW_OUT}"
    printf "AMR_TASK_INITIALIZED\tDeclared raw output placeholder created before scanner execution.\n" >> "${RAW_OUT}"
    printf "<html><body><h3>Fungal AMR summary</h3><p>AMR task initialized.</p></body></html>\n" > "${HTML_OUT}"
    : > "${LOG_OUT}"

    echo "[$(date -u)] Starting fungal AMR screen for ${SAMPLE}" >> "${LOG_OUT}"

    # rc146 broad Candida AMR interpretation:
    # This version is not limited to C. auris positive controls.  It assigns a
    # species-aware antifungal-resistance locus panel for all common Candida and
    # renamed Candida taxa, then uses that panel to interpret scanner no-hit
    # results conservatively.
    KNOWN_POSITIVE_CONTROL=0
    CANDIDA_AMR_PANEL="ERG11,TAC1/TAC1B,UPC2,MRR1,PDR1,CDR1/CDR2,MDR1,FKS1/FKS2,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1"
    echo "rc146 broad Candida AMR mode: enabled" >> "${LOG_OUT}"
    echo "Assembly FASTA: ${ASM}" >> "${LOG_OUT}"
    echo "Species summary TSV: ${SPECIES_TSV}" >> "${LOG_OUT}"
    echo "Threads: ${THREADS}" >> "${LOG_OUT}"

    # v9 TERM + SCANNER STATUS FIX: ChroQueTas uses tput internally. Cromwell/Docker non-interactive
    # jobs often have TERM unset, causing repeated: "tput: No value for $TERM" and
    # a non-zero scanner exit even when the executable and species mapping are valid.
    # Export a safe terminal type before invoking the scanner.
    export TERM="${TERM:-xterm}"
    echo "AMR container PATH: ${PATH}" >> "${LOG_OUT}"
    echo "AMR container TERM: ${TERM}" >> "${LOG_OUT}"
    command -v ChroQueTas.sh >> "${LOG_OUT}" 2>&1 || true
    command -v run_fungamr_scan >> "${LOG_OUT}" 2>&1 || true
    command -v miniprot >> "${LOG_OUT}" 2>&1 || true

    # rc147 dependency hardening:
    # ChroQueTas requires miniprot. The production AMR Docker image must provide it.
    MINIPROT_MISSING=0
    if command -v miniprot >/dev/null 2>&1; then
      echo "Confirmed dependency: miniprot at $(command -v miniprot)" >> "${LOG_OUT}"
      miniprot --version >> "${LOG_OUT}" 2>&1 || true
    else
      MINIPROT_MISSING=1
      echo "DEPENDENCY_MISSING: miniprot is not available in the configured AMR Docker image." >> "${LOG_OUT}"
    fi
    fi

    # -------------------------------------------------------------------------
    # v5 FIX: use a ChroQueTas-enabled AMR Docker image and keep shell-only parsing.
    # The previous WDL patch protected Cromwell outputs but the old Docker image
    # did not contain ChroQueTas in PATH. This WDL defaults to the rebuilt image
    # gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9. Parsing remains POSIX shell/awk.
    # -------------------------------------------------------------------------

    SPECIES="NA"
    if [ -s "${SPECIES_TSV}" ]; then
      SPECIES=$(awk -F '\t' '
        NR==1 {
          for (i=1; i<=NF; i++) {
            h=tolower($i); gsub(/[^a-z0-9]+/, "_", h); col[h]=i
          }
          next
        }
        NR>1 {
          if ("top_species" in col && $(col["top_species"]) != "") {print $(col["top_species"]); exit}
          if ("species" in col && $(col["species"]) != "") {print $(col["species"]); exit}
          if ("species_identified" in col && $(col["species_identified"]) != "") {print $(col["species_identified"]); exit}
          if ("taxon" in col && $(col["taxon"]) != "") {print $(col["taxon"]); exit}
          if (NF >= 2 && $2 != "") {print $2; exit}
          if (NF >= 1 && $1 != "") {print $1; exit}
        }
      ' "${SPECIES_TSV}" 2>/dev/null)
      [ -n "${SPECIES}" ] || SPECIES="NA"
    fi
    echo "Parsed species for AMR report: ${SPECIES}" >> "${LOG_OUT}"

    # rc143 alias-normalization fix:
    # ChroQueTas/FungAMR requires the exact species key used in its installed
    # screen list.  Some clinically familiar Candida names have been renamed
    # taxonomically; this block accepts either old or new names and passes the
    # supported key to the AMR wrapper.
    normalize_fungal_species_for_chroquetas() {
      sp="$1"
      sp=$(printf "%s" "$sp" | sed 's/^ *//; s/ *$//; s/[[:space:]]\+/_/g; s/-/_/g')
      sp_lc=$(printf "%s" "$sp" | tr '[:upper:]' '[:lower:]')
      case "$sp_lc" in
        candida_auris|candidozyma_auris|c_auris|c._auris) echo "Candidozyma_auris" ;;
        candida_glabrata|nakaseomyces_glabratus|nakaseomyces_glabrata|c_glabrata|c._glabrata) echo "Nakaseomyces_glabratus" ;;
        candida_krusei|pichia_kudriavzevii|pichia_kudriavzevii|c_krusei|c._krusei) echo "Pichia_kudriavzevii" ;;
        candida_lusitaniae|clavispora_lusitaniae|c_lusitaniae|c._lusitaniae) echo "Clavispora_lusitaniae" ;;
        candida_albicans|c_albicans|c._albicans) echo "Candida_albicans" ;;
        candida_dubliniensis|c_dubliniensis|c._dubliniensis) echo "Candida_dubliniensis" ;;
        candida_metapsilosis|c_metapsilosis|c._metapsilosis) echo "Candida_metapsilosis" ;;
        candida_orthopsilosis|c_orthopsilosis|c._orthopsilosis) echo "Candida_orthopsilosis" ;;
        candida_parapsilosis|c_parapsilosis|c._parapsilosis) echo "Candida_parapsilosis" ;;
        candida_tropicalis|c_tropicalis|c._tropicalis) echo "Candida_tropicalis" ;;
        *) printf "%s" "$sp" ;;
      esac
    }

    AMR_SPECIES=$(normalize_fungal_species_for_chroquetas "${SPECIES}")

    # Species-aware AMR locus panels used for reporting and validation. These are
    # screening targets, not automatic resistance calls. A marker is only reported
    # as detected if the configured scanner output contains marker/variant evidence.
    case "$(printf "%s" "${AMR_SPECIES}" | tr '[:upper:]' '[:lower:]')" in
      candidozyma_auris|candida_auris)
        CANDIDA_AMR_PANEL="ERG11,TAC1B,FKS1,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1,copy-number/aneuploidy/LOH"
        ;;
      nakaseomyces_glabratus|candida_glabrata)
        CANDIDA_AMR_PANEL="PDR1,CDR1/CDR2/SNQ2,FKS1/FKS2,ERG11,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1,copy-number/aneuploidy/LOH"
        ;;
      candida_albicans|candida_dubliniensis)
        CANDIDA_AMR_PANEL="ERG11,TAC1,UPC2,MRR1,MDR1,CDR1/CDR2,FKS1,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1"
        ;;
      candida_tropicalis)
        CANDIDA_AMR_PANEL="ERG11,TAC1,UPC2,MRR1/MDR1,CDR1/CDR2,FKS1,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1"
        ;;
      candida_parapsilosis|candida_orthopsilosis|candida_metapsilosis)
        CANDIDA_AMR_PANEL="ERG11,MRR1/MDR1,CDR1/CDR2,FKS1,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1"
        ;;
      pichia_kudriavzevii|candida_krusei)
        CANDIDA_AMR_PANEL="ERG11,ABC/MFS efflux regulators,FKS1,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1; note intrinsic fluconazole non-susceptibility is species-associated"
        ;;
      clavispora_lusitaniae|candida_lusitaniae)
        CANDIDA_AMR_PANEL="ERG11,ERG pathway genes,FKS1,FCY1/FCY2/FUR1; amphotericin B resistance may involve ergosterol-pathway changes"
        ;;
      *)
        CANDIDA_AMR_PANEL="ERG11,TAC/UPC2/MRR/PDR-family efflux regulators,CDR/MDR efflux genes,FKS1/FKS2,ERG2/ERG3/ERG6,FCY1/FCY2/FUR1"
        ;;
    esac
    echo "rc146 species-aware AMR panel: ${CANDIDA_AMR_PANEL}" >> "${LOG_OUT}"

    AMR_SPECIES_TSV="${SAMPLE}.species_for_fungamr.normalized.tsv"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "sample_id" "top_species" "species" "species_identified" "taxon" "raw_species" "amr_species" > "${AMR_SPECIES_TSV}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${AMR_SPECIES}" "${AMR_SPECIES}" "${AMR_SPECIES}" "${AMR_SPECIES}" "${SPECIES}" "${AMR_SPECIES}" >> "${AMR_SPECIES_TSV}"
    echo "Normalized species submitted to AMR scanner: ${AMR_SPECIES}" >> "${LOG_OUT}"
    echo "Normalized AMR species TSV: ${AMR_SPECIES_TSV}" >> "${LOG_OUT}"
    cat "${AMR_SPECIES_TSV}" >> "${LOG_OUT}" 2>/dev/null || true

    # rc144 runtime wrapper patch:
    # The Docker image wrapper may internally remap current species keys back to
    # unsupported legacy aliases, e.g. Candidozyma_auris -> Candida_auris.
    # Patch a local copy of run_fungamr_scan inside the Cromwell execution
    # directory and put it first in PATH, so both old and new Candida names are
    # submitted using the exact keys available in ChroQueTas/FungAMR.
    if command -v run_fungamr_scan >/dev/null 2>&1; then
      ORIG_RUN_FUNGAMR_SCAN=$(command -v run_fungamr_scan)
      mkdir -p patched_amr_bin

      # rc148 ChroQueTas stale-output shim:
      # Some container wrapper versions pre-create the sample .ChroQueTas output
      # directory before invoking ChroQueTas.sh. ChroQueTas itself then aborts with:
      #   ERROR: amr_out/<sample>.ChroQueTas already exist! Please check
      # A pre-run rm in the WDL is not enough, because the wrapper can recreate
      # the directory internally. Put a ChroQueTas.sh shim earlier in PATH that
      # removes any *.ChroQueTas output directory passed to ChroQueTas immediately
      # before delegating to the real executable.
      ORIG_CHROQUETAS_SH=""
      if command -v ChroQueTas.sh >/dev/null 2>&1; then
        ORIG_CHROQUETAS_SH=$(command -v ChroQueTas.sh)
        cat > patched_amr_bin/ChroQueTas.sh <<EOF_CHROQ_SHIM
#!/usr/bin/env bash
set +e
for arg in "\$@"; do
  case "\$arg" in
    *.ChroQueTas|*.ChroQueTas/|*/ChroQueTas|*/ChroQueTas/|*'.ChroQueTas'*)
      if [ -e "\$arg" ]; then
        rm -rf "\$arg" 2>/dev/null || true
      fi
      ;;
  esac
done
exec "${ORIG_CHROQUETAS_SH}" "\$@"
EOF_CHROQ_SHIM
        chmod +x patched_amr_bin/ChroQueTas.sh 2>> "${LOG_OUT}" || true
        echo "Runtime ChroQueTas stale-output shim created for: ${ORIG_CHROQUETAS_SH}" >> "${LOG_OUT}"
      else
        echo "WARNING: ChroQueTas.sh not found while creating stale-output shim." >> "${LOG_OUT}"
      fi

      cp "${ORIG_RUN_FUNGAMR_SCAN}" patched_amr_bin/run_fungamr_scan 2>> "${LOG_OUT}" || true
      chmod +x patched_amr_bin/run_fungamr_scan 2>> "${LOG_OUT}" || true
      if [ -s patched_amr_bin/run_fungamr_scan ]; then
        sed -i.bak \
          -e 's/Candida_auris/Candidozyma_auris/g' \
          -e 's/Candida glabrata/Nakaseomyces glabratus/g' \
          -e 's/Candida_glabrata/Nakaseomyces_glabratus/g' \
          -e 's/Candida krusei/Pichia kudriavzevii/g' \
          -e 's/Candida_krusei/Pichia_kudriavzevii/g' \
          -e 's/Candida lusitaniae/Clavispora lusitaniae/g' \
          -e 's/Candida_lusitaniae/Clavispora_lusitaniae/g' \
          patched_amr_bin/run_fungamr_scan 2>> "${LOG_OUT}" || true
        export PATH="${PWD}/patched_amr_bin:${PATH}"
        echo "Runtime-patched run_fungamr_scan placed first in PATH: $(command -v run_fungamr_scan)" >> "${LOG_OUT}"
        echo "Runtime alias patch summary:" >> "${LOG_OUT}"
        grep -nE 'Candidozyma_auris|Nakaseomyces_glabratus|Pichia_kudriavzevii|Clavispora_lusitaniae|Candida_auris|Candida_glabrata|Candida_krusei|Candida_lusitaniae' patched_amr_bin/run_fungamr_scan >> "${LOG_OUT}" 2>/dev/null || true
      else
        echo "WARNING: Could not create runtime-patched copy of run_fungamr_scan; using container original." >> "${LOG_OUT}"
      fi
    fi

    # rc148 pre-scan stale-output protection:
    # ChroQueTas refuses to run if its output directory already exists. Remove
    # sample-specific scanner artifacts before each run so Cromwell retries and
    # local/manual tests do not inherit stale outputs.
    rm -rf \
      "amr_out/${SAMPLE}.ChroQueTas" \
      "amr_out/${SAMPLE}.ChroQueTas.stdout" \
      "amr_out/${SAMPLE}.ChroQueTas.stderr" \
      "amr_out/${SAMPLE}.chroquetas.list_species.txt" \
      "amr_out/${SAMPLE}.chroquetas.list_species.stderr" \
      "amr_out/${SAMPLE}.species_for_chroquetas.tsv" \
      "amr_out/${SAMPLE}.fungal_amr.summary.tsv" \
      "amr_out/${SAMPLE}.fungal_amr.raw.tsv" \
      "amr_out/${SAMPLE}.fungal_amr.html" 2>> "${LOG_OUT}" || true

    if [ "${MINIPROT_MISSING}" = "1" ]; then
      SCAN_RC=127
      echo "DEPENDENCY_MISSING: miniprot is required by ChroQueTas but is not installed in the AMR Docker image. Build/pull an AMR image that includes miniprot, then rerun." > amr_out/scanner.stderr
      : > amr_out/scanner.stdout
    elif command -v run_fungamr_scan >/dev/null 2>&1; then
      echo "Confirmed CLI: $(command -v run_fungamr_scan)" >> "${LOG_OUT}"
      run_fungamr_scan \
        --sample "${SAMPLE}" \
        --assembly "${ASM}" \
        --species-summary "${AMR_SPECIES_TSV}" \
        --threads "${THREADS}" \
        --outdir amr_out \
        > amr_out/scanner.stdout 2> amr_out/scanner.stderr
      SCAN_RC=$?
    else
      SCAN_RC=127
      echo "run_fungamr_scan was not found in PATH inside the AMR Docker image." > amr_out/scanner.stderr
      : > amr_out/scanner.stdout
    fi

    echo "[$(date -u)] run_fungamr_scan exit code: ${SCAN_RC}" >> "${LOG_OUT}"
    echo "--- scanner stdout ---" >> "${LOG_OUT}"
    cat amr_out/scanner.stdout >> "${LOG_OUT}" 2>/dev/null || true
    echo "--- scanner stderr ---" >> "${LOG_OUT}"
    cat amr_out/scanner.stderr >> "${LOG_OUT}" 2>/dev/null || true
    echo "--- AMR output files ---" >> "${LOG_OUT}"
    find amr_out -maxdepth 6 -type f -print >> "${LOG_OUT}" 2>/dev/null || true

    # Always write raw scanner harvest without Python.
    printf "source_file\traw_text\n" > "${RAW_OUT}"
    HARVESTED=0
    find amr_out -maxdepth 6 -type f \( -name "*.tsv" -o -name "*.txt" -o -name "*.csv" -o -name "*.out" \) -print 2>/dev/null | sort | while IFS= read -r f; do
      base=$(basename "$f")
      case "$base" in
        scanner.stdout|scanner.stderr) continue ;;
      esac
      if [ -s "$f" ]; then
        printf "%s\t" "$f" >> "${RAW_OUT}"
        # single-line, tab-safe representation of the raw output
        tr '\t\r\n' '   ' < "$f" >> "${RAW_OUT}"
        printf "\n" >> "${RAW_OUT}"
      fi
    done

    if [ $(wc -l < "${RAW_OUT}" 2>/dev/null || echo 0) -gt 1 ]; then
      HARVESTED=1
    else
      printf "NO_SCANNER_TABLES\tNo non-empty TSV/TXT/CSV/OUT scanner output files were harvested from amr_out. See fungal_amr.log.\n" >> "${RAW_OUT}"
      HARVESTED=0
    fi

    # v9 scanner-status normalization:
    # ChroQueTas can return exit 1 in non-interactive Cromwell/Docker runs even when
    # it produced only list/species/output artifacts and no true diagnostic error.
    # Do not mark the biological AMR result as scanner_failed unless stderr contains
    # a real error after removing harmless terminal/tput warnings.
    REAL_ERR=""
    FATAL_SCANNER_MSG=""
    if [ -s amr_out/scanner.stderr ]; then
      REAL_ERR=$(grep -v -E '^tput: No value for \$TERM|^tput: No value for \$TERM and no -T specified$|^[[:space:]]*$' amr_out/scanner.stderr 2>/dev/null | head -n 1)
    fi
    if ls amr_out/*.ChroQueTas.stderr >/dev/null 2>&1; then
      QT_ERR=$(grep -v -E '^tput: No value for \$TERM|^tput: No value for \$TERM and no -T specified$|^[[:space:]]*$' amr_out/*.ChroQueTas.stderr 2>/dev/null | head -n 1)
      [ -n "${QT_ERR}" ] && REAL_ERR="${QT_ERR}"
    fi

    # rc148 fatal-pattern detection across stdout and stderr:
    # ChroQueTas sometimes prints real errors to stdout, not stderr. Do not
    # treat non-zero scanner exits as completed when stdout contains dependency,
    # stale-output, or command failures.
    FATAL_SCANNER_MSG=$(cat \
        amr_out/scanner.stderr \
        amr_out/scanner.stdout \
        amr_out/*.ChroQueTas.stderr \
        amr_out/*.ChroQueTas.stdout 2>/dev/null \
      | grep -v -E '^tput: No value for \$TERM|^tput: No value for \$TERM and no -T specified$|^[[:space:]]*$' \
      | grep -Ei '(^|[[:space:]])ERROR:|required and not installed|not installed|command not found|No such file|already exist|Traceback|Exception|failed|cannot|permission denied|DEPENDENCY_MISSING' \
      | head -n 1)
    if [ -n "${FATAL_SCANNER_MSG}" ]; then
      REAL_ERR="${FATAL_SCANNER_MSG}"
      echo "Fatal scanner message detected: ${FATAL_SCANNER_MSG}" >> "${LOG_OUT}"
    fi

    if [ "${SCAN_RC}" != "0" ] && [ -z "${REAL_ERR}" ]; then
      echo "Scanner returned ${SCAN_RC}, but no real stderr/stdout error was detected after filtering harmless tput/TERM messages; treating as completed for reporting." >> "${LOG_OUT}"
      SCAN_RC=0
    fi

    # Detect ChroQueTas/FungAMR species-support and nomenclature failures explicitly.
    # This prevents a true wrapper/species-mapping failure from being collapsed into
    # a misleading biological "No hit" row, which is especially important for
    # positive-control isolates.
    UNSUPPORTED_SPECIES_LINE=""
    CHROQ_SPECIES=""
    RAW_SPECIES_FOR_CHROQ=""
    LISTED_CANDIDOZYMA_AURIS=0
    if [ "${HARVESTED}" = "1" ]; then
      UNSUPPORTED_SPECIES_LINE=$(tail -n +2 "${RAW_OUT}" 2>/dev/null | grep -Ei 'unsupported_species|Species not available in the installed ChroQueTas/FungAMR screen list' | head -n 1)
      CHROQ_SPECIES=$(tail -n +2 "${RAW_OUT}" 2>/dev/null | grep -Ei 'species_for_chroquetas' | sed -E 's/.*[[:space:]](Candida_auris|Candidozyma_auris)([[:space:]].*)?$/\1/' | head -n 1)
      RAW_SPECIES_FOR_CHROQ=$(tail -n +2 "${RAW_OUT}" 2>/dev/null | grep -Ei 'species_for_chroquetas' | head -n 1 | cut -c1-500)
      if tail -n +2 "${RAW_OUT}" 2>/dev/null | grep -Eq 'chroquetas.list_species.*Candidozyma_auris'; then
        LISTED_CANDIDOZYMA_AURIS=1
      fi
    fi
    echo "Unsupported species line: ${UNSUPPORTED_SPECIES_LINE}" >> "${LOG_OUT}"
    echo "Species mapping line: ${RAW_SPECIES_FOR_CHROQ}" >> "${LOG_OUT}"
    echo "ChroQueTas mapped species: ${CHROQ_SPECIES}" >> "${LOG_OUT}"
    echo "ChroQueTas list contains Candidozyma_auris: ${LISTED_CANDIDOZYMA_AURIS}" >> "${LOG_OUT}"

    # Minimal marker detection using shell only.  This is intentionally conservative.
    # v9: exclude ChroQueTas species-list/support files from marker detection because
    # they list searchable proteins (Cyp51, Fks, etc.) and are not resistance calls.
    MARKER_LINE=""
    if [ "${HARVESTED}" = "1" ]; then
      MARKER_LINE=$(tail -n +2 "${RAW_OUT}"         | grep -Eiv 'chroquetas.list_species|species_for_chroquetas|scanner.stdout|scanner.stderr|no hit|no marker|none detected|not detected|negative|absence|source_file'         | grep -Ei 'ERG11|ERG1|ERG2|ERG3|ERG4|ERG5|ERG6|TAC1|TAC1B|UPC2|MRR1|PDR1|CDR1|CDR2|SNQ2|MDR1|FKS1|FKS2|HS1|FCY1|FCY2|FUR1|azole|echinocandin|amphotericin|flucytosine|fluconazole|voriconazole|posaconazole|itraconazole|caspofungin|micafungin|anidulafungin|resistan|mutation|variant|hotspot|copy.number|aneuploid|LOH|loss.of.heterozygosity'         | head -n 1)
    fi

    # rc151 marker parser hardening:
    # ChroQueTas can successfully detect a curated FungAMR marker but the broad
    # text harvester may only report it as a generic candidate. Parse the actual
    # ChroQueTas marker tables first, then fall back to the older text heuristic.
    PARSED_GENE=""
    PARSED_MUTATION=""
    PARSED_EFFECT=""
    PARSED_DRUGS=""
    PARSED_SOURCE=""
    PARSED_EVIDENCE=""
    PARSED_DRUG_CLASS="azole/other"
    PARSED_DRUG="fluconazole/other"

    parse_drug_class_and_label() {
      drug_text="$1"
      drug_lc=$(printf "%s" "$drug_text" | tr '[:upper:]' '[:lower:]')
      classes=""
      labels=""
      if printf "%s" "$drug_lc" | grep -Eq 'fluconazole|voriconazole|posaconazole|itraconazole|isavuconazole|azole'; then
        classes="azole"
        labels="fluconazole/other azoles"
      fi
      if printf "%s" "$drug_lc" | grep -Eq 'caspofungin|micafungin|anidulafungin|echinocandin'; then
        [ -n "$classes" ] && classes="${classes}/echinocandin" || classes="echinocandin"
        [ -n "$labels" ] && labels="${labels}/echinocandins" || labels="echinocandins"
      fi
      if printf "%s" "$drug_lc" | grep -Eq 'amphotericin|polyene'; then
        [ -n "$classes" ] && classes="${classes}/polyene" || classes="polyene"
        [ -n "$labels" ] && labels="${labels}/amphotericin_B" || labels="amphotericin_B"
      fi
      if printf "%s" "$drug_lc" | grep -Eq '5-fluorocytosine|flucytosine'; then
        [ -n "$classes" ] && classes="${classes}/flucytosine" || classes="flucytosine"
        [ -n "$labels" ] && labels="${labels}/flucytosine" || labels="flucytosine"
      fi
      [ -n "$classes" ] || classes="antifungal"
      [ -n "$labels" ] || labels="species-aware antifungal panel"
      printf "%s	%s\n" "$classes" "$labels"
    }

    # Preferred parser: per-protein ChroQueTas TSV files. These retain exact
    # Position, Reference AA, Query AA, Result, and Fungicides columns.
    PARSED_ROW=$(find amr_out -type f -name "${SAMPLE}.contigs.ChroQueTaS.*.tsv" ! -name '*AMR*' -print 2>/dev/null | sort | while IFS= read -r f; do
      awk -v f="$f" '
        BEGIN { FS="[ 	]+"; OFS="	" }
        NR > 1 && $0 ~ /FungAMR[ 	]+MUTATION/ {
          gene=f
          sub(/^.*ChroQueTaS\./, "", gene)
          sub(/\.[0-9]+\.tsv$/, "", gene)

          # ChroQueTas per-protein TSVs can appear whitespace-delimited rather
          # than strict tab-delimited in Cromwell outputs. Parse both forms.
          # Expected row:
          # Position Reference Query Result Fungicides
          # 132      Y         F     FungAMR MUTATION Fluconazole(...)
          pos=$1
          ref=$2
          qry=$3
          eff=$4 " " $5
          drugs=""
          for (i=6; i<=NF; i++) {
            drugs = drugs (drugs=="" ? "" : " ") $i
          }
          if (drugs == "") drugs="NA"
          mut=ref pos qry
          print gene, mut, eff, drugs, f ":" $0
          exit
        }
      ' "$f"
    done | head -n 1)

    # Fallback parser: consolidated ChroQueTas AMR summary files.
    if [ -z "${PARSED_ROW}" ]; then
      for f in \
        "amr_out/${SAMPLE}.ChroQueTas.AMR_summary.tsv" \
        "amr_out/${SAMPLE}.ChroQueTas/${SAMPLE}.contigs.ChroQueTaS.AMR_summary.txt" \
        "amr_out/${SAMPLE}.ChroQueTas/${SAMPLE}.contigs.ChroQueTaS.AMR_summary.tsv"; do
        if [ -s "$f" ]; then
          PARSED_ROW=$(awk -v f="$f" '
            NR > 1 && NF >= 6 {
              gene=$1; pos=$3; ref=$4; qry=$5; drugs=$6
              for (i=7; i<=NF; i++) drugs=drugs " " $i
              mut=ref pos qry
              print gene "\t" mut "\tFungAMR MUTATION\t" drugs "\t" f ":" $0
              exit
            }
          ' "$f" | head -n 1)
          [ -n "${PARSED_ROW}" ] && break
        fi
      done
    fi

    if [ -n "${PARSED_ROW}" ]; then
      PARSED_GENE=$(printf "%s" "${PARSED_ROW}" | awk -F '\t' '{print $1}')
      PARSED_MUTATION=$(printf "%s" "${PARSED_ROW}" | awk -F '\t' '{print $2}')
      PARSED_EFFECT=$(printf "%s" "${PARSED_ROW}" | awk -F '\t' '{print $3}')
      PARSED_DRUGS=$(printf "%s" "${PARSED_ROW}" | awk -F '\t' '{print $4}')
      PARSED_SOURCE=$(printf "%s" "${PARSED_ROW}" | awk -F '\t' '{print $5}' | tr '\r\n' ' ' | cut -c1-700)
      PARSED_EVIDENCE="FungAMR curated marker detected by ChroQueTas"
      DRUG_PAIR=$(parse_drug_class_and_label "${PARSED_DRUGS}")
      PARSED_DRUG_CLASS=$(printf "%s" "${DRUG_PAIR}" | awk -F '\t' '{print $1}')
      PARSED_DRUG=$(printf "%s" "${DRUG_PAIR}" | awk -F '\t' '{print $2}')
      echo "Parsed exact ChroQueTas marker: gene=${PARSED_GENE}; mutation=${PARSED_MUTATION}; effect=${PARSED_EFFECT}; drugs=${PARSED_DRUGS}" >> "${LOG_OUT}"
    else
      echo "No exact ChroQueTas FungAMR marker row parsed; using scanner status fallback." >> "${LOG_OUT}"
    fi

    printf "sample_id\tspecies\tdrug_class\tdrug\tgene_or_status\tmutation\teffect\tevidence_level\tinterpretation\n" > "${SUMMARY_OUT}"

    if [ "${SCAN_RC}" != "0" ]; then
      ERRMSG=$(cat amr_out/scanner.stderr amr_out/scanner.stdout amr_out/*.ChroQueTas.stderr amr_out/*.ChroQueTas.stdout 2>/dev/null | tr '\t\r\n' '   ' | cut -c1-800)
      [ -n "${ERRMSG}" ] || ERRMSG="Scanner exited non-zero; see fungal_amr.log."
      if printf "%s" "${ERRMSG}" | grep -Eiq 'miniprot.*required|miniprot.*not installed|DEPENDENCY_MISSING'; then
        EVIDENCE_CODE="SCANNER_DEPENDENCY_MISSING_MINIPROT; exit_${SCAN_RC}"
        STATUS_LABEL="Scanner dependency missing"
      elif printf "%s" "${ERRMSG}" | grep -Eiq 'already exist'; then
        EVIDENCE_CODE="SCANNER_STALE_OUTPUT_DIRECTORY; exit_${SCAN_RC}"
        STATUS_LABEL="Scanner stale output directory"
      else
        EVIDENCE_CODE="SCANNER_FAILED; exit_${SCAN_RC}"
        STATUS_LABEL="Scanner failed"
      fi
      printf "%s\t%s\tazole/other\tfluconazole/other\t%s\tNA\tNA\t%s\tThe configured fungal AMR scanner did not complete successfully. This is not a biological no-hit result and should not be interpreted as susceptibility. Technical evidence: %s\n" \
        "${SAMPLE}" "${SPECIES}" "${STATUS_LABEL}" "${EVIDENCE_CODE}" "${ERRMSG}" >> "${SUMMARY_OUT}"
    elif [ -n "${UNSUPPORTED_SPECIES_LINE}" ]; then
      CLEAN_UNSUPPORTED=$(printf "%s" "${UNSUPPORTED_SPECIES_LINE}" | tr '\t\r\n' '   ' | cut -c1-600)
      if printf "%s" "${CLEAN_UNSUPPORTED}" | grep -q 'Candida_auris' && [ "${LISTED_CANDIDOZYMA_AURIS}" = "1" ]; then
        printf "%s\t%s\tazole/other\tfluconazole/other\tSpecies mapping mismatch\tNA\tNA\tSPECIES_MAPPING_MISMATCH_UNSUPPORTED_ALIAS\tAMR validation failed before marker screening: the wrapper mapped the sample to Candida_auris, but the installed ChroQueTas/FungAMR species list contains Candidozyma_auris. This is a nomenclature/alias mapping bug, not a true no-hit and not a susceptible call. Patch the AMR container wrapper so Candida auris/Candidozyma auris is submitted as Candidozyma_auris, then rerun the positive controls. Raw evidence: %s\n" \
          "${SAMPLE}" "${SPECIES}" "${CLEAN_UNSUPPORTED}" >> "${SUMMARY_OUT}"
      else
        printf "%s\t%s\tazole/other\tfluconazole/other\tSpecies unsupported by AMR scanner\tNA\tNA\tUNSUPPORTED_SPECIES\tThe AMR scanner did not screen this isolate because the mapped species was not available in the installed ChroQueTas/FungAMR screen list. This is not a biological no-hit and not a susceptible call. Raw evidence: %s\n" \
          "${SAMPLE}" "${SPECIES}" "${CLEAN_UNSUPPORTED}" >> "${SUMMARY_OUT}"
      fi
    elif [ -n "${PARSED_GENE}" ]; then
      CLEAN_SOURCE=$(printf "%s" "${PARSED_SOURCE}" | tr '\t\r\n' '   ' | cut -c1-650)
      CLEAN_DRUGS=$(printf "%s" "${PARSED_DRUGS}" | tr '\t\r\n' '   ' | cut -c1-350)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tCurated ChroQueTas/FungAMR marker detected: %s %s (%s). Drug evidence from ChroQueTas: %s. Source: %s\n" \
        "${SAMPLE}" "${SPECIES}" "${PARSED_DRUG_CLASS}" "${PARSED_DRUG}" "${PARSED_GENE}" "${PARSED_MUTATION}" "${PARSED_EFFECT}" "${PARSED_EVIDENCE}" "${PARSED_GENE}" "${PARSED_MUTATION}" "${PARSED_EFFECT}" "${CLEAN_DRUGS}" "${CLEAN_SOURCE}" >> "${SUMMARY_OUT}"
    elif [ -n "${MARKER_LINE}" ]; then
      CLEAN_MARKER=$(printf "%s" "${MARKER_LINE}" | tr '\t\r\n' '   ' | cut -c1-600)
      printf "%s\t%s\tazole/other\tfluconazole/other\tCandidate AMR marker\tSee raw TSV\tNA\tscanner_detected_candidate_marker\tThe AMR scanner produced text consistent with a possible curated resistance marker, but rc151 could not parse an exact ChroQueTas marker row. Review fungal_amr.raw.tsv and fungal_amr.log for the exact record: %s\n" \
        "${SAMPLE}" "${SPECIES}" "${CLEAN_MARKER}" >> "${SUMMARY_OUT}"

    else
      printf "%s\t%s\tazole/echinocandin/polyene/flucytosine\tspecies-aware antifungal panel\tNo curated marker detected\tNA\tNA\tSCANNER_COMPLETED_NO_CURATED_MARKER_VARIANT_SCREEN_RECOMMENDED\tThe FungAMR/ChroQueTas wrapper completed and no curated genomic AMR marker was detected in the harvested output. This is not a susceptible call. For this species, rc146 flags the following resistance loci/mechanisms for variant-aware follow-up: %s. Phenotypic resistance may involve target-gene SNPs, promoter/regulatory changes, efflux regulation, copy-number change, LOH/aneuploidy, or markers absent from the installed database.\n" \
        "${SAMPLE}" "${SPECIES}" "${CANDIDA_AMR_PANEL}" >> "${SUMMARY_OUT}"
    fi

    # Build a small standalone HTML using shell only, preserving the parsed marker fields.
    OUT_SPECIES=$(awk -F '\t' 'NR==2 {print $2}' "${SUMMARY_OUT}" 2>/dev/null)
    OUT_DRUG_CLASS=$(awk -F '\t' 'NR==2 {print $3}' "${SUMMARY_OUT}" 2>/dev/null)
    OUT_DRUG=$(awk -F '\t' 'NR==2 {print $4}' "${SUMMARY_OUT}" 2>/dev/null)
    STATUS=$(awk -F '\t' 'NR==2 {print $5}' "${SUMMARY_OUT}" 2>/dev/null)
    OUT_MUTATION=$(awk -F '\t' 'NR==2 {print $6}' "${SUMMARY_OUT}" 2>/dev/null)
    OUT_EFFECT=$(awk -F '\t' 'NR==2 {print $7}' "${SUMMARY_OUT}" 2>/dev/null)
    EVIDENCE=$(awk -F '\t' 'NR==2 {print $8}' "${SUMMARY_OUT}" 2>/dev/null)
    INTERP=$(awk -F '\t' 'NR==2 {print $9}' "${SUMMARY_OUT}" 2>/dev/null | sed 's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g')
    cat > "${HTML_OUT}" <<EOF_AMR_HTML
<html><body>
<h3>Fungal AMR summary</h3>
<table border="1" cellpadding="4" cellspacing="0">
<tr><th>sample_id</th><th>species</th><th>drug_class</th><th>drug</th><th>gene_or_status</th><th>mutation</th><th>effect</th><th>evidence_level</th><th>interpretation</th></tr>
<tr><td>${SAMPLE}</td><td>${OUT_SPECIES}</td><td>${OUT_DRUG_CLASS}</td><td>${OUT_DRUG}</td><td>${STATUS}</td><td>${OUT_MUTATION}</td><td>${OUT_EFFECT}</td><td>${EVIDENCE}</td><td>${INTERP}</td></tr>
</table>
</body></html>
EOF_AMR_HTML

    # Final guard: Cromwell must always find every declared output.
    [ -s "${SUMMARY_OUT}" ] || printf "sample_id\tspecies\tdrug_class\tdrug\tgene_or_status\tmutation\teffect\tevidence_level\tinterpretation\n%s\tNA\tazole/other\tfluconazole/other\tOutput guard\tNA\tNA\tOUTPUT_GUARD\tSummary output was missing or empty and was recreated by final guard.\n" "${SAMPLE}" > "${SUMMARY_OUT}"
    [ -s "${RAW_OUT}" ] || printf "source_file\traw_text\nOUTPUT_GUARD\tRaw output was missing or empty and was recreated by final guard.\n" > "${RAW_OUT}"
    [ -s "${HTML_OUT}" ] || printf "<html><body><h3>Fungal AMR summary</h3><p>HTML output recreated by final guard.</p></body></html>\n" > "${HTML_OUT}"
    [ -s "${LOG_OUT}" ] || printf "AMR log recreated by final guard.\n" > "${LOG_OUT}"

    echo "[$(date -u)] Final declared outputs:" >> "${LOG_OUT}"
    ls -lh "${SUMMARY_OUT}" "${RAW_OUT}" "${HTML_OUT}" "${LOG_OUT}" >> "${LOG_OUT}" 2>&1 || true

    exit 0
  >>>

  output {
    File amr_summary_tsv = "~{sample_id}.fungal_amr.summary.tsv"
    File amr_raw_tsv = "~{sample_id}.fungal_amr.raw.tsv"
    File amr_report_html = "~{sample_id}.fungal_amr.html"
    File amr_log = "~{sample_id}.fungal_amr.log"
  }

  runtime {
    docker: "~{fungal_amr_docker_image}"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    continueOnReturnCode: [0]
  }
}



task CANDIDA_EXPORT_REFERENCE_FASTAS {
  input {
    String refs_manifest = "/opt/rmap_candida_refs/references.tsv"
    String docker_image = "gmboowa/rmap-candida-refs:2026.05"
    Int cpu = 1
    Int memory_gb = 32
    Int disk_gb = 1000
  }

  command <<<
    set -euo pipefail
    mkdir -p refs_out

    python3 <<'PY'
from pathlib import Path
import csv, re, shutil

manifest = Path("~{refs_manifest}")
if not manifest.exists():
    raise FileNotFoundError(f"Reference manifest not found inside Docker image: {manifest}")

def slug(x):
    s = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(x).strip())
    return s.strip('_') or 'species'

species_out = []
with manifest.open(newline="") as fh:
    sample = fh.read(4096)
    fh.seek(0)
    has_header = sample.lower().startswith("species\t")
    if has_header:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = []
        for row in reader:
            species = (row.get("species") or row.get("Species") or "").strip()
            ref = (row.get("reference") or row.get("reference_path") or row.get("Reference") or "").strip()
            if species and ref:
                rows.append((species, ref))
    else:
        rows = []
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                rows.append((parts[0].strip(), parts[1].strip()))

if not rows:
    raise RuntimeError(f"No reference rows were parsed from {manifest}")

for idx, (species, ref) in enumerate(rows):
    ref_path = Path(ref)
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference FASTA for {species} not found: {ref}")
    out = Path("refs_out") / f"{idx:03d}_{slug(species)}.fasta"
    shutil.copyfile(ref_path, out)
    species_out.append(species)

Path("reference_species.txt").write_text("\n".join(species_out) + "\n")
Path("reference_manifest_used.tsv").write_text(
    "species\treference_fasta_copied\n" +
    "\n".join(f"{sp}\trefs_out/{idx:03d}_{slug(sp)}.fasta" for idx, sp in enumerate(species_out)) +
    "\n"
)
PY
  >>>

  output {
    Array[String] reference_species = read_lines("reference_species.txt")
    Array[File] reference_fastas = glob("refs_out/*.fasta")
    File reference_manifest_used = "reference_manifest_used.tsv"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task CANDIDA_TREE_VISUALIZATION {
  input {
    File input_tree
    String species_label
    String docker_image = "gmboowa/ete3-render:1.18"
    Int width = 2600
    Int height = 1600
    String image_format = "png"
    Int cpu = 2
    Int memory_gb = 32
    Int disk_gb = 1000
  }

  command <<<
    set -euo pipefail
    mkdir -p tree_visualization
    cp "~{input_tree}" tree_visualization/input.treefile
    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg

    cat > tree_visualization/render_tree.py <<'PY_RENDER'
from pathlib import Path
import re, traceback, zlib, struct, math

outdir = Path("tree_visualization")
outdir.mkdir(exist_ok=True)
species_label = '~{species_label}'.strip() or "Candida species"
slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", species_label).strip("_") or "Candida_species"
tree_path = Path("tree_visualization/input.treefile")
out_img = outdir / f"{slug}.core_snp_tree.png"
cleaned_tree = outdir / f"{slug}.core_snp_tree.cleaned.nwk"
render_log = outdir / f"{slug}.tree_render.log"
requested_w = int('~{width}')
requested_h = int('~{height}')

newick_text = tree_path.read_text(errors="replace").strip() if tree_path.exists() else ""
has_internal_support_labels = bool(re.search(r"\)([0-9]+(?:\.[0-9]+)?(?:/[0-9]+(?:\.[0-9]+)?)?)(?=[:),;])", newick_text))

species_parts = species_label.split()
species_genus = species_parts[0] if len(species_parts) >= 1 else "Candida"
species_epithet = " ".join(species_parts[1:]) if len(species_parts) >= 2 else "species"


def clean_leaf_name(name):
    s = str(name).strip().strip("'\"")
    s = Path(s).name
    s = re.sub(r"^snippy[_-]+", "", s, flags=re.I)
    s = re.sub(r"^core[_-]+", "", s, flags=re.I)
    s = re.sub(r"(\.sorted)?\.bam$", "", s, flags=re.I)
    s = re.sub(r"\.(fastq|fq)(\.gz)?$", "", s, flags=re.I)
    s = re.sub(r"\.(vcf|bcf)(\.gz)?$", "", s, flags=re.I)
    s = re.sub(r"\.(consensus|fa|fasta|fna|aln|treefile|nwk)$", "", s, flags=re.I)
    s = re.sub(r"_R?[12](_001)?$", "", s, flags=re.I)
    return re.sub(r"[\s]+", "_", s) or "sample"


def is_reference_tip(name):
    raw = str(name).strip().strip("'\"")
    clean = clean_leaf_name(raw).lower()
    raw_low = raw.lower()
    explicit = {"ref", "reference", "reference_genome", "reference-genome", "outgroup", "root", "reference_sequence", "reference-sequence"}
    return clean in explicit or raw_low in explicit or clean.startswith(("reference_", "reference-", "ref_", "ref-", "outgroup_", "outgroup-"))


def write_log(lines):
    with open(render_log, "w") as log:
        for x in lines:
            log.write(str(x).rstrip() + "\n")


def looks_blank_png(path, min_bytes=9000):
    try:
        data = path.read_bytes()
        if len(data) < min_bytes:
            return True
        if data[:8] != b"\x89PNG\r\n\x1a\n":
            return False
        pos = 8
        raw = bytearray()
        width = height = None
        color_type = None
        bit_depth = None
        while pos + 8 <= len(data):
            ln = struct.unpack(">I", data[pos:pos+4])[0]
            typ = data[pos+4:pos+8]
            chunk = data[pos+8:pos+8+ln]
            pos += 12 + ln
            if typ == b"IHDR":
                width, height, bit_depth, color_type = struct.unpack(">IIBB", chunk[:10])[:5]
            elif typ == b"IDAT":
                raw.extend(chunk)
            elif typ == b"IEND":
                break
        if not raw or not width or not height:
            return True
        dec = zlib.decompress(bytes(raw))
        channels = 4 if color_type == 6 else 3 if color_type == 2 else 1
        stride = width * channels + 1
        rows_to_check = min(height, 250)
        sample = dec[:stride * rows_to_check]
        # If almost all bytes are white/zero filters, image is effectively blank.
        darkish = sum(1 for b in sample if 0 < b < 235)
        return darkish < max(60, width * rows_to_check * channels * 0.0008)
    except Exception:
        return False


def draw_fallback_png(reason):
    """Robust non-empty fallback. This is only used if ETE/Qt fails."""
    try:
        from PIL import Image, ImageDraw, ImageFont
        W = max(1800, requested_w)
        H = max(1000, requested_h)
        img = Image.new("RGB", (W, H), "white")
        d = ImageDraw.Draw(img)
        try:
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 42)
            tip_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 32)
            species_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf", 24)
            note_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
        except Exception:
            title_font = tip_font = species_font = note_font = None
        d.text((70, 50), f"rMAP-Candida Core-SNP Phylogenetic Tree: {species_label}", fill="black", font=title_font)
        d.text((70, 105), f"Fallback renderer used: {reason}", fill=(120, 30, 30), font=note_font)
        leaves = []
        for m in re.finditer(r"([A-Za-z0-9_.:-]+)(?=:[0-9.eE+-]+)", newick_text):
            v = clean_leaf_name(m.group(1))
            if v and not is_reference_tip(v) and v not in leaves:
                leaves.append(v)
        if not leaves:
            leaves = ["Tree file parsed but no non-reference tips were detected"]
        top = 200
        step = max(58, min(110, int((H - 320) / max(1, len(leaves)))))
        trunk_x = 140
        tip_x = 420
        y_positions = []
        for i, leaf in enumerate(leaves):
            y = top + i * step
            y_positions.append(y)
            d.line((trunk_x, y, tip_x, y), fill="black", width=4)
            d.ellipse((tip_x-7, y-7, tip_x+7, y+7), fill=(30, 96, 210))
            d.text((tip_x+16, y-23), leaf, fill=(20, 34, 55), font=tip_font)
            sample_w = d.textlength(leaf, font=tip_font) if hasattr(d, "textlength") else len(leaf) * 18
            d.text((tip_x+35+sample_w, y-17), species_label, fill=(37, 99, 235), font=species_font)
        if len(y_positions) > 1:
            d.line((trunk_x, y_positions[0], trunk_x, y_positions[-1]), fill="black", width=4)
        d.line((90, H-90, 230, H-90), fill="black", width=3)
        d.text((90, H-75), "not to scale", fill="black", font=note_font)
        img.save(out_img)
    except Exception:
        Path(out_img).write_bytes(b"fallback render failed\n")


log_lines = []
try:
    if not newick_text:
        raise RuntimeError("Tree input was empty or missing.")

    from ete3 import Tree, TreeStyle, NodeStyle, TextFace
    t = Tree(str(tree_path), format=0)

    # Remove explicit reference/outgroup leaves only; keep all real samples.
    for leaf in list(t.iter_leaves()):
        if is_reference_tip(leaf.name):
            leaf.detach()

    for leaf in t.iter_leaves():
        leaf.name = clean_leaf_name(leaf.name)

    n_leaves = len(t.get_leaves())
    if n_leaves < 2:
        raise RuntimeError(f"Fewer than two non-reference tips available after cleaning: {n_leaves}")

    t.write(outfile=str(cleaned_tree), format=0)

    # Adaptive rendering:
    #   Small species trees (3-5 tips) should not be drawn on a very large canvas.
    #   Larger species trees keep the wider/taller layout needed for readability.
    if n_leaves <= 5:
        layout_class = "compact"
        canvas_w = min(max(requested_w, 1400), 1800)
        canvas_h = 820
        font_size = 24
        title_font_size = 28
        species_font_size = 16
        support_font_size = 15
        dot_size = 8
        line_width = 3
        margin_left = 65
        margin_right = 360
        margin_top = 45
        margin_bottom = 65
        branch_vertical_margin = 16
    elif n_leaves <= 12:
        layout_class = "medium"
        canvas_w = min(max(requested_w, 1800), 2200)
        canvas_h = max(950, 520 + 70 * n_leaves)
        font_size = 26
        title_font_size = 32
        species_font_size = 18
        support_font_size = 16
        dot_size = 9
        line_width = 3
        margin_left = 65
        margin_right = 470
        margin_top = 60
        margin_bottom = 80
        branch_vertical_margin = 14
    elif n_leaves <= 30:
        layout_class = "standard"
        canvas_w = max(requested_w, 1850 + 55 * n_leaves)
        canvas_h = max(requested_h, 520 + 62 * n_leaves)
        font_size = 24
        title_font_size = 34
        species_font_size = 17
        support_font_size = 16
        dot_size = 9
        line_width = 3
        margin_left = 70
        margin_right = 560
        margin_top = 70
        margin_bottom = 90
        branch_vertical_margin = 10
    else:
        layout_class = "large"
        canvas_w = max(requested_w, 2000 + 42 * n_leaves)
        canvas_h = max(requested_h, 500 + 44 * n_leaves)
        font_size = 18
        title_font_size = 32
        species_font_size = 14
        support_font_size = 14
        dot_size = 8
        line_width = 2
        margin_left = 70
        margin_right = 620
        margin_top = 70
        margin_bottom = 90
        branch_vertical_margin = 7

    for n in t.traverse():
        ns = NodeStyle()
        ns["fgcolor"] = "#2563eb" if n.is_leaf() else "#111827"
        ns["size"] = dot_size if n.is_leaf() else 0
        ns["hz_line_width"] = line_width
        ns["vt_line_width"] = line_width
        n.set_style(ns)

    def format_support_label(node):
        """Show IQ-TREE/UFBoot/bootstrap labels on internal branches.
        ETE may store branch supports in node.support rather than node.name; use both.
        Avoid displaying ETE's artificial default support=1.0 unless the Newick actually contained support labels.
        """
        if node.is_leaf():
            return ""
        candidates = []
        nm = str(getattr(node, "name", "") or "").strip()
        if nm and nm not in {"NoName", "None"}:
            candidates.append(nm)
        sup = getattr(node, "support", None)
        if sup is not None and has_internal_support_labels:
            try:
                val = float(sup)
                if val > 0:
                    if 0 < val <= 1.0:
                        val = val * 100.0
                    candidates.append(str(val))
            except Exception:
                pass
        for c in candidates:
            c = str(c).strip()
            if not c or c in {"1", "1.0", "0", "0.0"} and not has_internal_support_labels:
                continue
            if "/" in c:
                return c
            try:
                v = float(c)
                if 0 < v <= 1.0:
                    v *= 100.0
                if v < 1.0:
                    return ""
                return str(int(round(v))) if abs(v - round(v)) < 0.05 else f"{v:.1f}"
            except Exception:
                return c
        return ""

    def layout(node):
        if node.is_leaf():
            sample_face = TextFace(node.name, fsize=font_size, fgcolor="#111827")
            sample_face.margin_left = 8
            sample_face.margin_right = 12
            node.add_face(sample_face, column=0, position="branch-right")
            sp_face = TextFace(species_label, fsize=species_font_size, fgcolor="#2563eb")
            # Italicize species names; this works in ETE3 with Qt fonts.
            try:
                sp_face.fstyle = "italic"
            except Exception:
                pass
            sp_face.margin_left = 4
            node.add_face(sp_face, column=1, position="branch-right")
        else:
            lab = format_support_label(node)
            if lab:
                boot = TextFace(lab, fsize=support_font_size, fgcolor="#dc2626")
                boot.margin_right = 6
                boot.margin_left = 2
                node.add_face(boot, column=0, position="branch-top")

    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.layout_fn = layout
    ts.margin_left = margin_left
    ts.margin_right = margin_right
    ts.margin_top = margin_top
    ts.margin_bottom = margin_bottom
    ts.branch_vertical_margin = branch_vertical_margin
    ts.scale = None
    title_face = TextFace(f"rMAP-Candida Core-SNP Phylogenetic Tree: {species_label}", fsize=title_font_size, fgcolor="#000000")
    title_face.margin_bottom = 12 if n_leaves <= 5 else 18
    ts.title.add_face(title_face, column=0)
    try:
        ts.show_scale = True
    except Exception:
        pass

    t.render(str(out_img), w=canvas_w, h=canvas_h, units="px", tree_style=ts)
    if looks_blank_png(out_img):
        log_lines.append(f"ETE output appeared blank or too small ({out_img.stat().st_size if out_img.exists() else 0} bytes); using fallback renderer.")
        draw_fallback_png("ETE produced a blank/too-small PNG")

    log_lines += [
        "Tree rendering completed.",
        f"species_label={species_label}",
        f"tips_rendered={n_leaves}",
        f"layout_class={layout_class}",
        f"has_internal_support_labels={has_internal_support_labels}",
        f"canvas={canvas_w}x{canvas_h}",
        f"output={out_img}",
        f"output_bytes={out_img.stat().st_size if out_img.exists() else 0}",
    ]
except Exception as e:
    log_lines.append("Primary tree rendering failed; writing non-empty fallback PNG.")
    log_lines.append(str(e))
    log_lines.append(traceback.format_exc())
    draw_fallback_png(str(e))

write_log(log_lines)
PY_RENDER

    python3 tree_visualization/render_tree.py
  >>>

  output {
    File tree_image = glob("tree_visualization/*.core_snp_tree.png")[0]
    File cleaned_newick = glob("tree_visualization/*.core_snp_tree.cleaned.nwk")[0]
    File render_log = glob("tree_visualization/*.tree_render.log")[0]
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task CANDIDA_SNIPPY_CORE_BY_SPECIES {
  input {
    Array[String]+ sample_names
    Array[File]+ read1s
    Array[File]+ read2s
    Array[File] species_top_tsvs
    Array[String] reference_species
    Array[File] reference_fastas

    # True Snippy branch species. Default is C. auris because it is relatively clonal
    # and Snippy/snippy-core is widely used for C. auris short-read WGS comparisons.
    Array[String] snippy_species = ["Candida auris", "Candidozyma auris", "Candida albicans", "Candida tropicalis"]

    # Non-Snippy species in this list use bcftools --ploidy 1.
    # All other non-Snippy species use diploid-aware bcftools --ploidy 2.
    Array[String] haploid_species = ["Nakaseomyces glabratus"]

    Int min_species_samples_for_tree = 3
    String docker_image = "staphb/snippy:4.6.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 1000
    Int min_quality = 20
    Int min_base_quality = 20
    Int min_depth = 10
    Int min_variant_quality = 20
    Float core_site_min_fraction = 0.95
  }

  command <<<
    set -euo pipefail
    mkdir -p phylogeny refs logs

    python3 <<'PY'
from pathlib import Path
import csv, re, shutil, subprocess, os, sys, gzip

samples = """~{sep='\n' sample_names}""".splitlines()
r1s = """~{sep='\n' read1s}""".splitlines()
r2s = """~{sep='\n' read2s}""".splitlines()
species_tsvs = """~{sep='\n' species_top_tsvs}""".splitlines()
ref_species = """~{sep='\n' reference_species}""".splitlines()
ref_files = """~{sep='\n' reference_fastas}""".splitlines()
haploid_species = """~{sep='\n' haploid_species}""".splitlines()
snippy_species = """~{sep='\n' snippy_species}""".splitlines()

min_n = int("~{min_species_samples_for_tree}")
cpu = int("~{cpu}")
min_mapq = int("~{min_quality}")
min_baseq = int("~{min_base_quality}")
min_depth = int("~{min_depth}")
min_vqual = int("~{min_variant_quality}")
core_fraction = float("~{core_site_min_fraction}")

def canonical_species_name(x):
    s = str(x or '').strip()
    k = re.sub(r'[^a-z0-9]+', ' ', s.lower()).strip()
    # NCBI/Kraken may report C. auris as Candidozyma auris, while older
    # references and reports may use Candida auris. Treat them as the same
    # species for reference lookup, Snippy routing, and report grouping.
    if k in {'candidozyma auris', 'candida auris'}:
        return 'Candida auris'
    return s

def norm(x):
    return re.sub(r'\s+', ' ', canonical_species_name(x)).lower()

def slug(x):
    s = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(x).strip())
    return s.strip('_') or 'species'

def read_top_species(path):
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                return row.get("top_species") or row.get("species") or row.get("name") or "Unknown"
    except Exception:
        return "Unknown"
    return "Unknown"

def run_cmd(cmd, log_path, cwd=None, shell=False):
    with open(log_path, "a") as lf:
        if shell:
            lf.write("\n$ " + str(cmd) + "\n")
            lf.flush()
            rc = subprocess.call(["bash", "-lc", str(cmd)], cwd=cwd, stdout=lf, stderr=subprocess.STDOUT)
        else:
            lf.write("\n$ " + " ".join(map(str, cmd)) + "\n")
            lf.flush()
            rc = subprocess.call(list(map(str, cmd)), cwd=cwd, stdout=lf, stderr=subprocess.STDOUT)
        lf.write(f"\n[exit_code] {rc}\n")
    return rc

def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    return seqs

def write_fasta(records, path, width=80):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def concatenate_genome_fasta(path):
    seqs = read_fasta(path)
    return "".join(seqs[k] for k in sorted(seqs))

def build_core_snp_alignment(consensus_paths, core_out, full_out):
    records = []
    lengths = set()
    for sample, cpath in consensus_paths:
        seq = concatenate_genome_fasta(cpath)
        records.append((sample, seq))
        lengths.add(len(seq))
    if len(lengths) != 1:
        raise RuntimeError("Consensus FASTA lengths differ after reference-based consensus generation.")
    L = lengths.pop()
    write_fasta(records, full_out)
    keep = []
    min_called = int(core_fraction * len(records) + 0.999999)
    valid = set("ACGT")
    for pos in range(L):
        col = [seq[pos].upper() for _, seq in records]
        called = sum(1 for x in col if x in valid)
        if called < min_called:
            continue
        alleles = set(x for x in col if x in valid)
        if len(alleles) >= 2:
            keep.append(pos)
    snp_records = [(sample, "".join(seq[i] for i in keep) if keep else "N") for sample, seq in records]
    write_fasta(snp_records, core_out)
    return len(keep), L

def clean_alignment_sample_name(name):
    s = str(name).strip().strip("'\"")
    s = Path(s).name
    s = re.sub(r"^snippy[_-]+", "", s, flags=re.I)
    s = re.sub(r"^core[_-]+", "", s, flags=re.I)
    s = re.sub(r"(\.sorted)?\.bam$", "", s, flags=re.I)
    s = re.sub(r"\.(fastq|fq)(\.gz)?$", "", s, flags=re.I)
    s = re.sub(r"\.(vcf|bcf)(\.gz)?$", "", s, flags=re.I)
    s = re.sub(r"\.(consensus|fa|fasta|fna|aln)$", "", s, flags=re.I)
    s = re.sub(r"_R?[12](_001)?$", "", s, flags=re.I)
    s = re.sub(r"\s+", "_", s)
    return s

def is_reference_alignment_name(name):
    # Conservative reference detection only.
    # Do NOT remove GCA/GCF/NC/NZ-like names here because user sample IDs may legally
    # look like accessions. Reference removal is instead driven by the expected
    # successful sample names whenever those are available.
    raw = str(name).strip().strip("'\"")
    clean = clean_alignment_sample_name(raw)
    low = clean.lower()
    raw_low = raw.lower()
    exact_reference_names = {
        "ref", "reference", "reference_genome", "outgroup",
        "core_ref", "core_reference", "snippy_ref", "snippy_reference"
    }
    if low in exact_reference_names or raw_low in exact_reference_names:
        return True
    if "reference" in raw_low or "reference" in low:
        return True
    return False

def strip_reference_from_alignment(infile, outfile, expected_samples=None):
    seqs = read_fasta(infile)
    expected = set(clean_alignment_sample_name(x) for x in (expected_samples or []) if str(x).strip())
    kept = []
    seen = {}

    for name, seq in seqs.items():
        clean = clean_alignment_sample_name(name)

        # Preferred behavior: keep only sequences corresponding to successful
        # input sample IDs. This preserves arbitrary sample names such as SRR*,
        # ERR*, sample-001, patient_A, GCA/GCF-like sample names, etc., and drops
        # the reference without relying on accession-pattern guessing.
        if expected:
            if clean not in expected:
                continue
        else:
            if is_reference_alignment_name(name):
                continue

        if clean in seen:
            seen[clean] += 1
            clean = f"{clean}_{seen[clean]}"
        else:
            seen[clean] = 1
        kept.append((clean, seq))

    # Defensive fallback: if name matching unexpectedly removes everything, fall
    # back to conservative reference-name filtering rather than failing silently.
    if not kept:
        for name, seq in seqs.items():
            if is_reference_alignment_name(name):
                continue
            clean = clean_alignment_sample_name(name)
            kept.append((clean, seq))

    write_fasta(kept, outfile)
    return len(kept), len(kept[0][1]) if kept else 0

refs = {norm(sp): (sp.strip(), rf.strip()) for sp, rf in zip(ref_species, ref_files) if sp.strip() and rf.strip()}
snippy_set = set(norm(x) for x in snippy_species if x.strip())
haploid_set = set(norm(x) for x in haploid_species if x.strip())

# rc172 patch: build phylogenies only for valid Candida species that have
# an explicit species-specific reference in the bundled/exported reference manifest.
# Invalid/unclassified labels must not become phylogeny groups because they have no
# biological reference and can trigger failed calls such as CANDIDA_SNIPPY_CORE_BY_SPECIES:NA.
INVALID_SPECIES_LABELS = {
    "", "na", "n/a", "nan", "none", "null", "unknown", "unclassified",
    "no_call", "no call", "no-call", "not determined", "not_determined",
    "undetermined", "unassigned", "unresolved", "no species", "no_species"
}

def is_invalid_species_label(sp):
    k = norm(sp)
    return k in INVALID_SPECIES_LABELS or k.startswith("unclassified") or k.startswith("unknown")

summary_rows = []
group_labels = []
alignment_manifest = []

groups = {}
skipped_invalid_counts = {}
skipped_no_reference_counts = {}
for sample, r1, r2, stsv in zip(samples, r1s, r2s, species_tsvs):
    sp = read_top_species(stsv)
    key = norm(sp)
    display = str(sp).strip() or "NA"

    if is_invalid_species_label(sp):
        skipped_invalid_counts.setdefault(display, 0)
        skipped_invalid_counts[display] += 1
        continue

    if key not in refs:
        # Species was detected, but it is not one of the defined Candida species
        # with a reference FASTA available for phylogeny. Keep it in the summary,
        # but do not run Snippy/bcftools/IQ-TREE for it.
        skipped_no_reference_counts.setdefault(display, 0)
        skipped_no_reference_counts[display] += 1
        continue

    groups.setdefault(key, {"display": display, "items": []})
    groups[key]["items"].append((sample, r1, r2))

for display, nskip in sorted(skipped_invalid_counts.items()):
    summary_rows.append({"species": display, "status": "SKIPPED_INVALID_OR_UNCLASSIFIED_SPECIES", "sample_count": str(nskip), "branch": "NA", "ploidy": "NA", "reference": "NA", "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "Skipped by rc172 phylogeny guard because the top species label was invalid/unclassified/unknown/NA. The workflow will continue and build trees for valid species groups only."})

for display, nskip in sorted(skipped_no_reference_counts.items()):
    summary_rows.append({"species": display, "status": "SKIPPED_SPECIES_NOT_IN_REFERENCE_MANIFEST", "sample_count": str(nskip), "branch": "NA", "ploidy": "NA", "reference": "NA", "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "Skipped by rc172 phylogeny guard because this detected species is not among the defined Candida species with an exported reference FASTA. Valid species with references are still processed."})

for key, info in sorted(groups.items(), key=lambda kv: kv[1]["display"]):
    display = info["display"]
    items = info["items"]
    n = len(items)
    if key not in refs:
        summary_rows.append({"species": display, "status": "SKIPPED_NO_REFERENCE", "sample_count": str(n), "branch": "NA", "ploidy": "NA", "reference": "NA", "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "No matching species-specific reference was provided in phylogeny_reference_species/phylogeny_reference_fastas."})
        continue
    ref_label, ref_path = refs[key]
    branch = "snippy_core" if key in snippy_set else ("haploid_bcftools_consensus" if key in haploid_set else "diploid_aware_bcftools_iupac_consensus")
    ploidy_str = "1" if key in snippy_set or key in haploid_set else "2"
    if n < min_n:
        summary_rows.append({"species": display, "status": "SKIPPED_TOO_FEW_SAMPLES", "sample_count": str(n), "branch": branch, "ploidy": ploidy_str, "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": f"Requires at least {min_n} samples for a species-specific tree."})
        continue
    gslug = slug(display)
    outdir = Path("phylogeny") / gslug
    outdir.mkdir(parents=True, exist_ok=True)
    ref_copy = outdir / "reference.fasta"
    shutil.copy(ref_path, ref_copy)

    if key in snippy_set:
        # rc170: rMAP-TB-style Snippy/core-SNP branch.
        # This branch is deliberately permissive and logs per-sample failures rather
        # than failing the whole workflow. It is now the default for C. auris and
        # C. albicans, fixing the earlier C. albicans bcftools variant_calling_failed
        # path that prevented a tree from being built.
        log = outdir / "snippy_core.log"
        failed, snippy_dirs, successful_samples = [], [], []

        snippy_bin = shutil.which("snippy") or "/usr/local/bin/snippy"
        snippy_core_bin = shutil.which("snippy-core") or "snippy-core"
        with open(log, "a") as lf:
            lf.write("CANDIDA_SNIPPY_CORE_BY_SPECIES rc170 Snippy branch\n")
            lf.write(f"Species: {display}\n")
            lf.write(f"Reference: {ref_label} -> {ref_copy}\n")
            lf.write(f"Using snippy: {snippy_bin}\n")
            lf.write(f"Using snippy-core: {snippy_core_bin}\n")
            lf.write(f"Minimum variant quality: {min_vqual}\n")
            lf.write(f"CPU threads: {cpu}\n\n")

        for sample, r1, r2 in items:
            sslug = slug(sample)
            sdir = outdir / sslug
            if sdir.exists():
                shutil.rmtree(sdir)

            cmd = [
                snippy_bin,
                "--cpus", str(cpu),
                "--minqual", str(min_vqual),
                "--ref", str(ref_copy),
                "--R1", str(r1),
                "--R2", str(r2),
                "--outdir", str(sdir),
                "--prefix", sslug,
                "--force"
            ]
            rc = run_cmd(cmd, log)

            vcf_candidates = [sdir / f"{sslug}.vcf", sdir / "snps.vcf"]
            has_vcf = any(x.exists() and x.stat().st_size > 0 for x in vcf_candidates)

            if rc != 0 or not has_vcf:
                failed.append(sample + ":snippy_failed")
            else:
                snippy_dirs.append(str(sdir.resolve()))
                successful_samples.append(sample)

        if len(snippy_dirs) < min_n:
            summary_rows.append({"species": display, "status": "SKIPPED_TOO_FEW_SUCCESSFUL_SNIPPY_RUNS", "sample_count": str(len(snippy_dirs)), "branch": "snippy_core_rc170_tb_style", "ploidy": "NA", "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "Fewer than the minimum number of samples completed Snippy. Failed: " + ",".join(failed)})
            continue

        cmd = [snippy_core_bin, "--ref", str(ref_copy.resolve()), "--prefix", "core"] + snippy_dirs
        rc = run_cmd(cmd, log, cwd=outdir)
        raw_core = outdir / "core.aln"
        raw_full = outdir / "core.full.aln"
        if rc != 0 or not raw_core.exists() or raw_core.stat().st_size == 0:
            summary_rows.append({"species": display, "status": "FAILED_SNIPPY_CORE", "sample_count": str(len(snippy_dirs)), "branch": "snippy_core_rc170_tb_style", "ploidy": "NA", "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "snippy-core failed or did not produce core.aln. Failed Snippy samples: " + ",".join(failed)})
            continue

        core_out = outdir / f"{gslug}.core_snps.aln"
        full_out = outdir / f"{gslug}.full_consensus.aln"
        nseqs, nsites = strip_reference_from_alignment(raw_core, core_out, expected_samples=successful_samples)
        if raw_full.exists() and raw_full.stat().st_size > 0:
            strip_reference_from_alignment(raw_full, full_out, expected_samples=successful_samples)
        else:
            shutil.copy(core_out, full_out)

        group_labels.append(display)
        alignment_manifest.append(str(core_out))
        summary_rows.append({"species": display, "status": "PASS", "sample_count": str(nseqs), "branch": "snippy_core_rc170_tb_style", "ploidy": "NA", "reference": ref_label, "core_alignment": str(core_out), "full_consensus_alignment": str(full_out), "variable_sites": str(nsites), "notes": "rc170 used the rMAP-TB-style Snippy -> snippy-core branch. Species grouping was based on Kraken2/Bracken top-species calls; mixed-species trees are intentionally avoided. Recombination is not explicitly filtered."})
        continue

    ploidy = 1 if key in haploid_set else 2
    strategy = "haploid_bcftools_consensus" if ploidy == 1 else "diploid_aware_bcftools_iupac_consensus"
    log = outdir / "mapping_variant_consensus.log"
    for idx_ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        try: (Path(str(ref_copy) + idx_ext)).unlink()
        except FileNotFoundError: pass
    if run_cmd(["bwa", "index", ref_copy], log) != 0:
        summary_rows.append({"species": display, "status": "FAILED_BWA_INDEX", "sample_count": str(n), "branch": strategy, "ploidy": str(ploidy), "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "bwa index failed for the supplied reference."}); continue
    if run_cmd(["samtools", "faidx", ref_copy], log) != 0:
        summary_rows.append({"species": display, "status": "FAILED_SAMTOOLS_FAIDX", "sample_count": str(n), "branch": strategy, "ploidy": str(ploidy), "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "samtools faidx failed for the supplied reference."}); continue
    consensus_paths, failed = [], []
    for sample, r1, r2 in items:
        sslug = slug(sample)
        sdir = outdir / sslug
        sdir.mkdir(exist_ok=True)
        slog = sdir / f"{sslug}.phylogeny.log"
        bam = sdir / f"{sslug}.sorted.bam"
        raw_vcf = sdir / f"{sslug}.raw.vcf.gz"
        filt_vcf = sdir / f"{sslug}.filtered.vcf.gz"
        lowcov_bed = sdir / f"{sslug}.low_coverage.bed"
        consensus = sdir / f"{sslug}.consensus.fasta"
        with open(slog, "w") as lf:
            lf.write(f"Species: {display}\nPloidy: {ploidy}\nStrategy: {strategy}\nReference: {ref_label}\n")
        if run_cmd(f"bwa mem -t {cpu} {ref_copy} {r1} {r2} | samtools sort -@ {cpu} -o {bam} -", slog, shell=True) != 0 or not bam.exists() or bam.stat().st_size == 0:
            failed.append(sample + ":mapping_failed"); continue
        if run_cmd(["samtools", "index", bam], slog) != 0:
            failed.append(sample + ":bam_index_failed"); continue
        mpileup_call = f"bcftools mpileup --threads {cpu} -Ou -q {min_mapq} -Q {min_baseq} -a FORMAT/DP,FORMAT/AD -f {ref_copy} {bam} | bcftools call --threads {cpu} --ploidy {ploidy} -mv -Oz -o {raw_vcf}"
        if run_cmd(mpileup_call, slog, shell=True) != 0 or not raw_vcf.exists() or raw_vcf.stat().st_size == 0:
            failed.append(sample + ":variant_calling_failed"); continue
        if run_cmd(["bcftools", "index", "-t", raw_vcf], slog) != 0:
            failed.append(sample + ":raw_vcf_index_failed"); continue
        filter_expr = f"QUAL>={min_vqual} && FMT/DP>={min_depth}"
        if run_cmd(["bcftools", "filter", "-i", filter_expr, "-Oz", "-o", filt_vcf, raw_vcf], slog) != 0:
            failed.append(sample + ":variant_filter_failed"); continue
        if run_cmd(["bcftools", "index", "-t", filt_vcf], slog) != 0:
            failed.append(sample + ":filtered_vcf_index_failed"); continue
        depth_cmd = f"samtools depth -aa -d 0 {bam} | awk -v md={min_depth} '{{if ($3 < md) print $1\"\\t\"($2-1)\"\\t\"$2}}' > {lowcov_bed}"
        if run_cmd(depth_cmd, slog, shell=True) != 0:
            failed.append(sample + ":depth_mask_failed"); continue
        consensus_cmd = f"bcftools consensus -f {ref_copy} -m {lowcov_bed} -I {filt_vcf} > {consensus}"
        if run_cmd(consensus_cmd, slog, shell=True) != 0 or not consensus.exists() or consensus.stat().st_size == 0:
            failed.append(sample + ":consensus_failed"); continue
        consensus_paths.append((sample, str(consensus)))
    if len(consensus_paths) < min_n:
        summary_rows.append({"species": display, "status": "SKIPPED_TOO_FEW_SUCCESSFUL_CONSENSUS", "sample_count": str(len(consensus_paths)), "branch": strategy, "ploidy": str(ploidy), "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": "Fewer than the minimum number of samples produced consensus FASTA files. Failed: " + ",".join(failed)})
        continue
    full_aln = outdir / f"{gslug}.full_consensus.aln"
    core_snps = outdir / f"{gslug}.core_snps.aln"
    try:
        nsites, full_len = build_core_snp_alignment(consensus_paths, core_snps, full_aln)
        group_labels.append(display)
        alignment_manifest.append(str(core_snps))
        summary_rows.append({"species": display, "status": "PASS", "sample_count": str(len(consensus_paths)), "branch": strategy, "ploidy": str(ploidy), "reference": ref_label, "core_alignment": str(core_snps), "full_consensus_alignment": str(full_aln), "variable_sites": str(nsites), "notes": f"{strategy}; low-depth sites masked; heterozygous diploid genotypes retained as IUPAC ambiguity codes. Recombination not explicitly filtered."})
    except Exception as e:
        summary_rows.append({"species": display, "status": "FAILED_CORE_SNP_ALIGNMENT", "sample_count": str(len(consensus_paths)), "branch": strategy, "ploidy": str(ploidy), "reference": ref_label, "core_alignment": "NA", "full_consensus_alignment": "NA", "variable_sites": "NA", "notes": str(e)})

with open("phylogeny/phylogeny_group_summary.tsv", "w", newline="") as out:
    fields = ["species", "status", "sample_count", "branch", "ploidy", "reference", "core_alignment", "full_consensus_alignment", "variable_sites", "notes"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in summary_rows:
        w.writerow(row)
with open("phylogeny/group_labels.txt", "w") as out:
    for label in group_labels:
        out.write(label + "\n")
with open("phylogeny/alignment_manifest.txt", "w") as out:
    for aln in alignment_manifest:
        out.write(aln + "\n")
PY
  >>>

  output {
    File phylogeny_group_summary_tsv = "phylogeny/phylogeny_group_summary.tsv"
    File group_labels_txt = "phylogeny/group_labels.txt"
    File alignment_manifest_txt = "phylogeny/alignment_manifest.txt"
    Array[File] core_full_alignments = glob("phylogeny/*/*.core_snps.aln")
    Array[File] full_consensus_alignments = glob("phylogeny/*/*.full_consensus.aln")
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


task COMPLEASM_FUNGAL {
  input {
    String sample_name
    File assembly_fasta
    String compleasm_docker
    String compleasm_lineage
    String compleasm_odb
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set +e
    set +u
    set +o pipefail

    SAMPLE="~{sample_name}"
    ASM="~{assembly_fasta}"
    LINEAGE="~{compleasm_lineage}"
    ODB="~{compleasm_odb}"
    SUMMARY_TSV="${SAMPLE}.busco_summary.tsv"
    SHORT_TXT="${SAMPLE}.busco.short_summary.txt"
    LOG="${SAMPLE}.compleasm.log"
    DETAILS="${SAMPLE}.compleasm.details.txt"
    OUTDIR="compleasm_${SAMPLE}"

    HEADER="sample_id	complete_busco_pct	single_copy_busco_pct	duplicated_busco_pct	fragmented_busco_pct	missing_busco_pct	busco_markers	busco_lineage	busco_status	busco_note"
    printf "%b\n" "${HEADER}" > "${SUMMARY_TSV}"
    : > "${LOG}"
    : > "${DETAILS}"

    echo "[$(date)] Compleasm task started for ${SAMPLE}" >> "${LOG}"
    echo "Assembly FASTA: ${ASM}" >> "${LOG}"
    echo "Lineage: ${LINEAGE}" >> "${LOG}"
    echo "ODB: ${ODB}" >> "${LOG}"
    echo "CPU: ~{cpu}; memory_gb: ~{memory_gb}" >> "${LOG}"
    command -v compleasm >> "${LOG}" 2>&1 || true
    compleasm -h >> "${LOG}" 2>&1 || true
    python3 --version >> "${LOG}" 2>&1 || true

    # Compleasm is a fast BUSCO-like completeness assessment based on miniprot.
    # The output TSV intentionally keeps BUSCO-compatible column names so the
    # existing rMAP-Candida HTML parser and summary tables continue to work.
    compleasm run \
      -a "${ASM}" \
      -o "${OUTDIR}" \
      -l "${LINEAGE}" \
      --odb "${ODB}" \
      -t "~{cpu}" >> "${LOG}" 2>&1
    RC=$?
    echo "[$(date)] Compleasm rc=${RC}" >> "${LOG}"

    SUMMARY_FILE="$(find "${OUTDIR}" -type f -name 'summary.txt' 2>/dev/null | sort | head -n 1)"
    if [ -n "${SUMMARY_FILE}" ] && [ -s "${SUMMARY_FILE}" ]; then
      cp "${SUMMARY_FILE}" "${SHORT_TXT}"
      python3 <<PYCOMP
import re, pathlib, csv
sample = "${SAMPLE}"
lineage = "${LINEAGE}_${ODB}"
summary_path = pathlib.Path("${SUMMARY_FILE}")
out_path = pathlib.Path("${SUMMARY_TSV}")
details_path = pathlib.Path("${DETAILS}")
txt = summary_path.read_text(errors="ignore")
details_path.write_text(txt)

vals = {
    "sample_id": sample,
    "complete_busco_pct": "NA",
    "single_copy_busco_pct": "NA",
    "duplicated_busco_pct": "NA",
    "fragmented_busco_pct": "NA",
    "missing_busco_pct": "NA",
    "busco_markers": "NA",
    "busco_lineage": lineage,
    "busco_status": "PASS",
    "busco_note": "Compleasm completed successfully. Values are BUSCO-like conserved-ortholog completeness metrics generated by Compleasm."
}

# Compleasm summary.txt usually contains lines such as:
# S:87.92%, 1419 / D:9.05%, 146 / F:1.73%, 28 / I:0.00%, 0 / M:1.30%, 21 / N:1614
parts = {}
for line in txt.splitlines():
    m = re.match(r"^\s*([SDFIMN])\s*:\s*([0-9.]+)\s*%?\s*,?\s*([0-9]+)?", line)
    if m:
        parts[m.group(1)] = (m.group(2), m.group(3) or "NA")

s_pct = float(parts.get("S", (0, "0"))[0]) if "S" in parts else None
d_pct = float(parts.get("D", (0, "0"))[0]) if "D" in parts else None
f_pct = float(parts.get("F", (0, "0"))[0]) if "F" in parts else None
i_pct = float(parts.get("I", (0, "0"))[0]) if "I" in parts else 0.0
m_pct = float(parts.get("M", (0, "0"))[0]) if "M" in parts else None
if s_pct is not None:
    vals["single_copy_busco_pct"] = f"{s_pct:.2f}"
if d_pct is not None:
    vals["duplicated_busco_pct"] = f"{d_pct:.2f}"
if f_pct is not None:
    # Compleasm has F and sometimes I fragmented subclasses. Combine them for a BUSCO-like fragmented value.
    vals["fragmented_busco_pct"] = f"{(f_pct + i_pct):.2f}"
if m_pct is not None:
    vals["missing_busco_pct"] = f"{m_pct:.2f}"
if s_pct is not None and d_pct is not None:
    vals["complete_busco_pct"] = f"{(s_pct + d_pct):.2f}"
if "N" in parts:
    vals["busco_markers"] = parts["N"][1] if parts["N"][1] != "NA" else parts["N"][0]

with out_path.open("w", newline="") as f:
    fieldnames = ["sample_id","complete_busco_pct","single_copy_busco_pct","duplicated_busco_pct","fragmented_busco_pct","missing_busco_pct","busco_markers","busco_lineage","busco_status","busco_note"]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerow(vals)
PYCOMP
    else
      echo "Compleasm failed or did not produce summary.txt. Inspect ${LOG}." > "${SHORT_TXT}"
      echo "Compleasm failed or did not produce summary.txt. Exit code: ${RC}." > "${DETAILS}"
      printf "%s\tNA\tNA\tNA\tNA\tNA\tNA\t%s_%s\tFAILED\tCompleasm failed or summary.txt was not found; inspect %s.\n" "${SAMPLE}" "${LINEAGE}" "${ODB}" "${LOG}" >> "${SUMMARY_TSV}"
    fi

    test -s "${SUMMARY_TSV}"
    exit 0
  >>>

  output {
    File busco_summary_tsv = "~{sample_name}.busco_summary.tsv"
    File busco_short_summary_txt = "~{sample_name}.busco.short_summary.txt"
    File compleasm_log = "~{sample_name}.compleasm.log"
    File compleasm_details = "~{sample_name}.compleasm.details.txt"
  }

  runtime {
    docker: "~{compleasm_docker}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task CANDIDA_IQTREE2_PHYLOGENY {
  input {
    File alignment
    String species_label
    String model = "GTR+G"
    Int bootstrap_replicates = 1000
    String docker_image = "gmboowa/iqtree2-python:2.3.4"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 1000
    Float max_missing_fraction_for_tree = 0.50
    Int min_non_reference_samples_for_tree = 3
  }

  command <<<
    set -uo pipefail
    mkdir -p iqtree logs

    echo "IQ-TREE report not available." > iqtree/iqtree.report
    echo "IQ-TREE log not available." > iqtree/iqtree.log
    echo "not_started" > iqtree/iqtree_status.txt
    echo "(IQTREE_not_started:0.0);" > iqtree/final.treefile
    echo "Alignment filtering not completed." > iqtree/alignment_filtering_summary.txt
    echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note" > iqtree/excluded_from_iqtree.tsv
    echo -e "sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction" > iqtree/included_in_iqtree.tsv

    python3 <<'PY'
import re
label = """~{species_label}"""
slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', label).strip('_') or 'Candida_species'
open("iqtree/group_slug.txt", "w").write(slug)
PY

    SLUG="$(cat iqtree/group_slug.txt)"
    cp "~{alignment}" "iqtree/${SLUG}.core.original.aln"
    cp "~{alignment}" "iqtree/${SLUG}.core.filtered.aln"

    python3 <<'PY'
from collections import OrderedDict
import sys, re

infile = "iqtree/" + open("iqtree/group_slug.txt").read().strip() + ".core.original.aln"
outfile = "iqtree/" + open("iqtree/group_slug.txt").read().strip() + ".core.filtered.aln"
excluded = "iqtree/excluded_from_iqtree.tsv"
included_tsv = "iqtree/included_in_iqtree.tsv"
summary = "iqtree/alignment_filtering_summary.txt"
max_missing = float("~{max_missing_fraction_for_tree}")
min_samples = int("~{min_non_reference_samples_for_tree}")
records=OrderedDict(); name=None; seq=[]
for line in open(infile, errors='replace'):
    line=line.strip()
    if not line: continue
    if line.startswith('>'):
        if name is not None: records[name]=''.join(seq)
        name=line[1:].split()[0]; seq=[]
    else:
        seq.append(line)
if name is not None: records[name]=''.join(seq)
lengths=set(map(len, records.values())) if records else set()
if not records or len(lengths)!=1:
    open(summary,'w').write('No usable equal-length FASTA alignment was available for IQ-TREE.\n')
    open(outfile,'w').write('')
    sys.exit(12)
L=next(iter(lengths))
inc=OrderedDict(); exc=[]
for n,s in records.items():
    su=s.upper(); acgt=sum(1 for b in su if b in 'ACGT'); miss=L-acgt; frac=miss/L if L else 1.0
    if L==0 or acgt==0 or frac>=max_missing:
        exc.append((n,L,acgt,miss,frac,'missing_or_no_acgt'))
    else:
        inc[n]=s
with open(excluded,'w') as out:
    out.write('sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\tthreshold\treason\texclusion_note\n')
    for n,L,acgt,miss,frac,reason in exc:
        out.write(f'{n}\t{L}\t{acgt}\t{miss}\t{frac:.6f}\t{max_missing:.6f}\t{reason}\tExcluded only from IQ-TREE because of high missing/ambiguous/gap content or no usable ACGT bases.\n')
with open(included_tsv,'w') as out:
    out.write('sample\talignment_length\tacgt_count\tmissing_count\tmissing_fraction\n')
    for n,s in inc.items():
        su=s.upper(); acgt=sum(1 for b in su if b in 'ACGT'); miss=L-acgt; frac=miss/L if L else 1.0
        out.write(f'{n}\t{L}\t{acgt}\t{miss}\t{frac:.6f}\n')
with open(outfile,'w') as out:
    for n,s in inc.items():
        out.write(f'>{n}\n')
        for i in range(0,len(s),80): out.write(s[i:i+80]+'\n')
with open(summary,'w') as out:
    out.write(f'Original sequences: {len(records)}\nIncluded sequences: {len(inc)}\nExcluded sequences: {len(exc)}\nAlignment length: {L}\n')
if len(inc) < min_samples:
    sys.exit(12)
PY

    FILTER_RC=$?
    echo "Alignment filtering exit code: ${FILTER_RC}" > logs/iqtree.command.log
    cat iqtree/alignment_filtering_summary.txt >> logs/iqtree.command.log || true

    if command -v iqtree2 >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree2)"
      IQTREE_MODE="host"
    elif command -v iqtree >/dev/null 2>&1; then
      IQTREE_BIN="$(command -v iqtree)"
      IQTREE_MODE="host"
    elif command -v docker >/dev/null 2>&1; then
      IQTREE_BIN="iqtree2"
      IQTREE_MODE="docker"
    else
      IQTREE_BIN=""
      IQTREE_MODE=""
    fi

    if [ "${FILTER_RC}" -eq 0 ] && [ -n "${IQTREE_BIN}" ]; then
      if [ "${IQTREE_MODE}" = "docker" ]; then
        # Run IQ-TREE inside Docker without relying on Cromwell-managed docker runtime.
        # The container entrypoint is overridden because some local backends/images
        # otherwise try to execute the wrong binary.
        IQTREE_CMD='set -euo pipefail
          if command -v iqtree2 >/dev/null 2>&1; then BIN="$(command -v iqtree2)";
          elif command -v iqtree >/dev/null 2>&1; then BIN="$(command -v iqtree)";
          else echo "ERROR: iqtree2/iqtree not found inside container" >&2; exit 127; fi
          if [ "~{bootstrap_replicates}" -gt 0 ]; then
            "$BIN" -s "iqtree/'"${SLUG}"'.core.filtered.aln" -m "~{model}" -B ~{bootstrap_replicates} -T ~{cpu} --prefix "iqtree/'"${SLUG}"'" -redo;
          else
            "$BIN" -s "iqtree/'"${SLUG}"'.core.filtered.aln" -m "~{model}" -T ~{cpu} --prefix "iqtree/'"${SLUG}"'" -redo;
          fi'
        docker run --rm --entrypoint /bin/bash \
          -v "$PWD:$PWD" \
          -w "$PWD" \
          "~{docker_image}" \
          -lc "${IQTREE_CMD}" > logs/iqtree.run.log 2>&1
      else
        if [ "~{bootstrap_replicates}" -gt 0 ]; then
          "${IQTREE_BIN}" -s "iqtree/${SLUG}.core.filtered.aln" -m "~{model}" -B ~{bootstrap_replicates} -T ~{cpu} --prefix "iqtree/${SLUG}" -redo > logs/iqtree.run.log 2>&1
        else
          "${IQTREE_BIN}" -s "iqtree/${SLUG}.core.filtered.aln" -m "~{model}" -T ~{cpu} --prefix "iqtree/${SLUG}" -redo > logs/iqtree.run.log 2>&1
        fi
      fi
      IQ_RC=$?
      if [ "${IQ_RC}" -eq 0 ] && [ -s "iqtree/${SLUG}.treefile" ]; then
        cp "iqtree/${SLUG}.treefile" iqtree/final.treefile
        [ -s "iqtree/${SLUG}.iqtree" ] && cp "iqtree/${SLUG}.iqtree" iqtree/iqtree.report
        [ -s "iqtree/${SLUG}.log" ] && cp "iqtree/${SLUG}.log" iqtree/iqtree.log
        echo "success; iqtree_bootstrap_replicates=~{bootstrap_replicates}" > iqtree/iqtree_status.txt
      else
        echo "iqtree_failed_after_filtering" > iqtree/iqtree_status.txt
      fi
    else
      echo "too_few_or_invalid_samples_after_filtering" > iqtree/iqtree_status.txt
      echo "IQ-TREE was skipped because too few valid samples remained after filtering, or IQ-TREE executable was unavailable." > logs/iqtree.run.log
    fi

    # Stable fallback/validation:
    # - Do not pass placeholder Newick such as (IQTREE_not_started:0.0) to the renderer.
    # - If IQ-TREE failed/skipped but a filtered alignment exists, build a small deterministic
    #   distance-based fallback Newick from the filtered alignment so the report always has a
    #   real sample tree instead of a broken image.
    python3 <<'PY'
from pathlib import Path
from collections import OrderedDict
import re, math

slug = Path('iqtree/group_slug.txt').read_text().strip()
final = Path('iqtree/final.treefile')
aln = Path(f'iqtree/{slug}.core.filtered.aln')
status = Path('iqtree/iqtree_status.txt')
runlog = Path('logs/iqtree.run.log')

def sanitize_name(name):
    name = str(name).strip().split()[0]
    name = re.sub(r'[^A-Za-z0-9_.-]+', '_', name)
    return name or "sample"

def read_fasta(path):
    records = OrderedDict()
    name = None
    seq = []
    if not path.exists():
        return records
    for line in path.read_text(errors='replace').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                records[sanitize_name(name)] = ''.join(seq).upper()
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records[sanitize_name(name)] = ''.join(seq).upper()
    return records

def looks_like_real_tree(txt):
    t = (txt or '').strip()
    if not t or not t.endswith(';'):
        return False
    bad_tokens = [
        'IQTREE_not_started',
        'IQTREE_no_tree_generated',
        'IQTREE_failed',
        'not_started',
        'no_tree_generated'
    ]
    if any(x in t for x in bad_tokens):
        return False
    # Require at least two named tips. This deliberately avoids full Newick parsing.
    labels = re.findall(r'([A-Za-z0-9_.-]+)\s*:', t)
    labels = [x for x in labels if not x.replace('.', '', 1).isdigit()]
    return len(set(labels)) >= 2

def hamming_distance(a, b):
    n = min(len(a), len(b))
    comparable = 0
    diff = 0
    for x, y in zip(a[:n], b[:n]):
        if x in 'ACGT' and y in 'ACGT':
            comparable += 1
            if x != y:
                diff += 1
    if comparable == 0:
        return 0.0
    return diff / comparable

def quote_newick_label(label):
    # Labels are sanitized, so quoting is not required.
    return sanitize_name(label)

def build_upgma(records):
    names = list(records.keys())
    if len(names) < 2:
        return None

    clusters = {}
    for n in names:
        clusters[n] = {
            'members': [n],
            'newick': quote_newick_label(n),
            'height': 0.0
        }

    # Precompute pairwise distances.
    pairdist = {}
    for i, a in enumerate(names):
        for b in names[i+1:]:
            pairdist[frozenset([a, b])] = hamming_distance(records[a], records[b])

    def avg_dist(c1, c2):
        vals = []
        for a in clusters[c1]['members']:
            for b in clusters[c2]['members']:
                if a == b:
                    continue
                vals.append(pairdist.get(frozenset([a, b]), 0.0))
        return sum(vals) / len(vals) if vals else 0.0

    counter = 0
    while len(clusters) > 1:
        keys = sorted(clusters.keys())
        best = None
        for i, a in enumerate(keys):
            for b in keys[i+1:]:
                d = avg_dist(a, b)
                candidate = (d, a, b)
                if best is None or candidate < best:
                    best = candidate
        d, a, b = best
        counter += 1
        parent_height = max(d / 2.0, clusters[a]['height'], clusters[b]['height'])
        la = max(parent_height - clusters[a]['height'], 0.0)
        lb = max(parent_height - clusters[b]['height'], 0.0)
        newick = f"({clusters[a]['newick']}:{la:.8f},{clusters[b]['newick']}:{lb:.8f})"
        members = clusters[a]['members'] + clusters[b]['members']
        del clusters[a]
        del clusters[b]
        clusters[f'cluster_{counter}'] = {
            'members': members,
            'newick': newick,
            'height': parent_height
        }

    root = next(iter(clusters.values()))
    return root['newick'] + ';\n'

current = final.read_text(errors='replace') if final.exists() else ''
if looks_like_real_tree(current):
    raise SystemExit(0)

records = read_fasta(aln)
fallback = build_upgma(records)
if fallback:
    final.write_text(fallback)
    previous_status = status.read_text(errors='replace').strip() if status.exists() else 'unknown'
    status.write_text(previous_status + '; fallback_distance_tree_generated\n')
    with runlog.open('a') as log:
        log.write('\nIQ-TREE final.treefile was missing, invalid, or placeholder.\n')
        log.write('Generated deterministic UPGMA-like fallback Newick from filtered core alignment.\n')
        log.write('Bootstrap labels are not available for fallback trees; rerun IQ-TREE successfully to obtain bootstrap values.\n')
        log.write(f'Fallback samples: {len(records)}\n')
else:
    # Last resort: keep a syntactically valid two-tip diagnostic tree so renderer does not break.
    final.write_text('(No_valid_tree_A:0.0,No_valid_tree_B:0.0);\n')
    previous_status = status.read_text(errors='replace').strip() if status.exists() else 'unknown'
    status.write_text(previous_status + '; diagnostic_two_tip_tree_generated\n')
    with runlog.open('a') as log:
        log.write('\nCould not generate fallback tree from alignment; wrote diagnostic two-tip tree.\n')
PY

    [ -s iqtree/iqtree.report ] || echo "IQ-TREE report not available; see iqtree_status.txt and logs/iqtree.run.log." > iqtree/iqtree.report
    [ -s iqtree/iqtree.log ] || echo "IQ-TREE log not available; see iqtree_status.txt and logs/iqtree.run.log." > iqtree/iqtree.log
  >>>

  output {
    File final_tree = "iqtree/final.treefile"
    File iqtree_report = "iqtree/iqtree.report"
    File iqtree_log = "iqtree/iqtree.log"
    File iqtree_status = "iqtree/iqtree_status.txt"
    File excluded_from_iqtree = "iqtree/excluded_from_iqtree.tsv"
    File included_in_iqtree = "iqtree/included_in_iqtree.tsv"
    File alignment_filtering_summary = "iqtree/alignment_filtering_summary.txt"
  }

  runtime {
    docker: "~{docker_image}"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


task MERGE_MYC_REPORTS {
  input {
    Array[String] sample_names
    File? multiqc_report
    Array[File] trimming_htmls
    Array[File] species_top_tsvs
    Array[File] kraken_reports
    Array[File] bracken_reports
    Array[File] assembly_summaries
    Array[File] quast_reports
    Array[File]? busco_summaries
    Array[File]? busco_summary_tsvs
    Array[File]? compleasm_summaries
    Array[File]? compleasm_summary_tsvs
    Array[File] amr_summaries
    Array[File] amr_htmls
    File? phylogeny_group_summary
    Array[File] phylogeny_core_alignments
    Array[File] phylogeny_newick_trees
    Array[File] phylogeny_iqtree_reports
    Array[File] phylogeny_tree_images
    File? surveillance_metadata_tsv
    Int cpu = 1
    Int memory_gb = 32
    Int disk_gb = 1000
  }

  command <<<
    set -euo pipefail

    mkdir -p report_inputs report_links
    cp ~{sep=' ' species_top_tsvs} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' kraken_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' bracken_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' assembly_summaries} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' quast_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' select_first([busco_summaries, []])} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' select_first([busco_summary_tsvs, []])} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' select_first([compleasm_summaries, []])} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' select_first([compleasm_summary_tsvs, []])} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' amr_summaries} report_inputs/ 2>/dev/null || true
    PHYLO_SUMMARY="~{default='' phylogeny_group_summary}"
    if [ -n "$PHYLO_SUMMARY" ] && [ -f "$PHYLO_SUMMARY" ]; then cp "$PHYLO_SUMMARY" report_inputs/ 2>/dev/null || true; fi
    cp ~{sep=' ' phylogeny_core_alignments} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' phylogeny_newick_trees} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' phylogeny_iqtree_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' phylogeny_tree_images} report_inputs/ 2>/dev/null || true
    METADATA_TSV="~{default='' surveillance_metadata_tsv}"
    if [ -n "$METADATA_TSV" ] && [ -f "$METADATA_TSV" ]; then cp "$METADATA_TSV" report_inputs/surveillance_metadata.tsv 2>/dev/null || true; fi
    cp ~{sep=' ' trimming_htmls} report_links/ 2>/dev/null || true
    cp ~{sep=' ' amr_htmls} report_links/ 2>/dev/null || true

    rm -f sample_names.txt
    for s in ~{sep=' ' sample_names}; do
      echo "$s" >> sample_names.txt
    done

    python3 <<'PY'
import csv, html, pathlib, re, datetime, statistics, base64
from collections import Counter

root = pathlib.Path("report_inputs")
samples = [x.strip() for x in pathlib.Path("sample_names.txt").read_text().splitlines() if x.strip()]

def esc(x):
    return html.escape("" if x is None else str(x))

def safe_float(x):
    try:
        if x is None or str(x).strip().upper() in {"", "NA", "NAN"}:
            return None
        return float(str(x).replace("%","").replace(",",""))
    except Exception:
        return None

def fmt_num(x):
    try:
        if x is None or str(x).strip().upper() in {"", "NA"}:
            return "NA"
        f = float(str(x).replace(",",""))
        return f"{int(f):,}" if f.is_integer() else f"{f:,.2f}"
    except Exception:
        return esc(x)

def pct_from_fraction(x):
    f = safe_float(x)
    if f is None:
        return "NA"
    if f <= 1:
        f *= 100
    return f"{f:.2f}"

def infer_sample(filename, suffixes):
    name = pathlib.Path(filename).name
    for s in suffixes:
        if name.endswith(s):
            return name[:-len(s)]
    for sample in samples:
        if sample in name:
            return sample
    return name.split(".")[0]

def read_dict_tsvs(pattern):
    out=[]
    for p in sorted(root.glob(pattern)):
        try:
            with p.open(errors="ignore", newline="") as f:
                r = csv.DictReader(f, delimiter="\t")
                for row in r:
                    row["_source_file"] = p.name
                    out.append(row)
        except Exception:
            pass
    return out

def first(row, keys, default="NA"):
    norm = {str(k).strip().lower().replace(" ","_").replace("-","_"): v for k,v in row.items()}
    for k in keys:
        kk = k.strip().lower().replace(" ","_").replace("-","_")
        v = norm.get(kk)
        if v is not None and str(v).strip() != "":
            return str(v).strip()
    return default

def parse_species():
    rows=[]
    for r in read_dict_tsvs("*.top_species.tsv"):
        sample = first(r, ["sample_id","sample"], infer_sample(r.get("_source_file",""), [".top_species.tsv"]))
        species = first(r, ["top_species","species","organism"], "Not determined")
        est = first(r, ["estimated_reads","new_est_reads","clade_reads","taxon_reads"], "NA")
        frac = first(r, ["fraction_total_reads","percent_reads","reads_percent","abundance"], "NA")
        rows.append({
            "sample_id": sample,
            "species": species,
            "percent_reads": pct_from_fraction(frac),
            "clade_reads": fmt_num(est),
            "taxon_reads": fmt_num(est),
            "taxid": first(r, ["taxid","taxonomy_id"], "NA"),
            "evidence": first(r, ["evidence"], "Bracken species-level abundance")
        })
    by_sample={r["sample_id"]:r for r in rows}
    for s in samples:
        by_sample.setdefault(s, {"sample_id":s,"species":"Not determined","percent_reads":"NA","clade_reads":"NA","taxon_reads":"NA","taxid":"NA","evidence":"No species TSV found"})
    return [by_sample[s] for s in samples if s in by_sample] + [r for r in rows if r["sample_id"] not in samples]

def parse_assembly():
    rows=[]
    for r in read_dict_tsvs("*.assembly_summary.tsv"):
        sample = first(r, ["sample_id","sample"], infer_sample(r.get("_source_file",""), [".assembly_summary.tsv"]))
        rows.append({
            "sample_id": sample,
            "contigs": fmt_num(first(r, ["contigs","num_contigs","#_contigs"], "NA")),
            "total_bp": fmt_num(first(r, ["total_bp","total_length","total_length_bp"], "NA")),
            "n50": fmt_num(first(r, ["n50","N50"], "NA")),
            "largest_contig": fmt_num(first(r, ["largest_contig","largest"], "NA"))
        })
    return rows

def parse_quast():
    rows=[]
    for p in sorted(root.glob("*.quast.report.tsv")):
        sample = infer_sample(p.name, [".quast.report.tsv"])
        metrics={}
        try:
            with p.open(errors="ignore", newline="") as f:
                table = [row for row in csv.reader(f, delimiter="\t") if row]
            if not table:
                continue
            header = [h.strip() for h in table[0]]
            header_norm = [h.lower().replace(" ","_") for h in header]
            if "sample_id" in header_norm:
                with p.open(errors="ignore", newline="") as f:
                    for r in csv.DictReader(f, delimiter="\t"):
                        sample = first(r, ["sample_id","sample"], sample)
                        metrics.update(r)
                        break
            elif header and header[0].strip().lower() == "assembly":
                value_col = 1 if len(header) > 1 else None
                for row in table[1:]:
                    if value_col is not None and len(row) > value_col:
                        metrics[row[0].strip()] = row[value_col].strip()
            else:
                # Last-resort metric/value parser
                for row in table:
                    if len(row) > 1:
                        metrics[row[0].strip()] = row[1].strip()
        except Exception:
            pass
        rows.append({
            "sample_id": sample,
            "# contigs": fmt_num(metrics.get("# contigs") or metrics.get("# contigs (>= 0 bp)") or metrics.get("# contigs (>= 500 bp)") or "NA"),
            "Largest contig": fmt_num(metrics.get("Largest contig") or "NA"),
            "Total length": fmt_num(metrics.get("Total length") or metrics.get("Total length (>= 0 bp)") or metrics.get("Total length (>= 500 bp)") or "NA"),
            "GC (%)": metrics.get("GC (%)", "NA"),
            "N50": fmt_num(metrics.get("N50") or "NA")
        })
    by={r["sample_id"]:r for r in rows}
    return [by.get(s,{"sample_id":s,"# contigs":"NA","Largest contig":"NA","Total length":"NA","GC (%)":"NA","N50":"NA"}) for s in samples]

def parse_busco():
    """
    Parse normalized Compleasm TSVs emitted with BUSCO-compatible column names.

    Fix applied:
    - accepts both <sample>.busco_summary.tsv and <sample>.busco.summary.tsv
    - accepts lowercase BUSCO task fields such as complete_busco_pct
    - falls back to parsing native BUSCO short_summary text
    """
    by={}
    patterns = ["*.busco_summary.tsv", "*.busco.summary.tsv", "*busco*summary*.tsv"]
    seen=set()
    rows=[]
    for pat in patterns:
        for r in read_dict_tsvs(pat):
            src = r.get("_source_file", "")
            key = (src, r.get("sample_id",""), r.get("sample",""))
            if key in seen:
                continue
            seen.add(key)
            rows.append(r)

    for r in rows:
        sample = first(r, ["sample_id", "sample"], infer_sample(r.get("_source_file", ""), [".busco_summary.tsv",".busco.summary.tsv"]))
        status = first(r, ["busco_status", "status"], "parsed")
        note = first(r, ["busco_note", "busco_details", "details"], "")
        by[sample]={
            "sample_id": sample,
            "Complete_BUSCO_%": first(r, ["complete_busco_pct","Complete_BUSCO_%","complete_busco","C"], "NA"),
            "Single_copy_BUSCO_%": first(r, ["single_copy_busco_pct","Single_copy_BUSCO_%","single_copy_busco","S"], "NA"),
            "Duplicated_BUSCO_%": first(r, ["duplicated_busco_pct","Duplicated_BUSCO_%","duplicated_busco","D"], "NA"),
            "Fragmented_BUSCO_%": first(r, ["fragmented_busco_pct","Fragmented_BUSCO_%","fragmented_busco","F"], "NA"),
            "Missing_BUSCO_%": first(r, ["missing_busco_pct","Missing_BUSCO_%","missing_busco","M"], "NA"),
            "BUSCO_n": first(r, ["busco_markers","BUSCO_n","n"], "NA"),
            "busco_status": status,
            "busco_details": note
        }

    if not by:
        for p in sorted(root.glob("*.busco.short_summary.txt")):
            sample = infer_sample(p.name, [".busco.short_summary.txt"])
            txt = p.read_text(errors="ignore")
            row = {
                "sample_id": sample, "Complete_BUSCO_%":"NA",
                "Single_copy_BUSCO_%":"NA", "Duplicated_BUSCO_%":"NA",
                "Fragmented_BUSCO_%":"NA", "Missing_BUSCO_%":"NA",
                "BUSCO_n":"NA", "busco_status":"not_parsed_or_failed",
                "busco_details":"No normalized Compleasm TSV was found."
            }
            m = re.search(r"C:([0-9.]+)%\s*\[\s*S:([0-9.]+)%\s*,\s*D:([0-9.]+)%\s*\]\s*,\s*F:([0-9.]+)%\s*,\s*M:([0-9.]+)%\s*,\s*n:(\d+)", txt)
            if m:
                row.update({
                    "Complete_BUSCO_%": m.group(1),
                    "Single_copy_BUSCO_%": m.group(2),
                    "Duplicated_BUSCO_%": m.group(3),
                    "Fragmented_BUSCO_%": m.group(4),
                    "Missing_BUSCO_%": m.group(5),
                    "BUSCO_n": m.group(6),
                    "busco_status":"parsed",
                    "busco_details":"Parsed from native BUSCO/Compleasm summary text."
                })
            by[sample]=row

    return [by.get(s,{"sample_id":s,"Complete_BUSCO_%":"NA","Single_copy_BUSCO_%":"NA","Duplicated_BUSCO_%":"NA","Fragmented_BUSCO_%":"NA","Missing_BUSCO_%":"NA","BUSCO_n":"NA","busco_status":"missing","busco_details":"No Compleasm output was found for this sample."}) for s in samples]

def parse_amr():
    rows=[]
    generic_nohit_values={"NA","N/A","","NO HIT","NO MARKER","NO MARKER DETECTED","NO CURATED MARKER DETECTED","CANDIDATE AMR MARKER","SCANNER FAILED"}

    def rescue_marker_from_interpretation(gene, mutation, effect, interp):
        # rc152: MERGE must not throw away a true ChroQueTas marker merely
        # because the per-sample summary lacks an explicit hit_status column.
        # Some earlier summaries put the exact marker only in interpretation:
        # "Curated ChroQueTas/FungAMR marker detected: Cyp51 K143R (FungAMR MUTATION) ..."
        g=(gene or "").strip()
        m=(mutation or "").strip()
        e=(effect or "").strip()
        text=interp or ""
        mm=re.search(r'detected:\s*([A-Za-z0-9_./-]+)\s+([A-Z*][0-9]+[A-Z*])\s*\(([^)]*)\)', text, flags=re.I)
        if mm:
            if g.upper() in generic_nohit_values or "NO MARKER" in g.upper() or "CANDIDATE" in g.upper():
                g=mm.group(1)
            if m.upper() in {"NA","N/A","","SEE RAW TSV"}:
                m=mm.group(2)
            if e.upper() in {"NA","N/A",""}:
                e=mm.group(3)
        return g, m, e

    for r in read_dict_tsvs("*.fungal_amr_summary.tsv") + read_dict_tsvs("*.fungal_amr.summary.tsv"):
        sample = first(r, ["sample_id","sample"], infer_sample(r.get("_source_file",""), [".fungal_amr_summary.tsv",".fungal_amr.summary.tsv"]))
        raw_hit_status = first(r, ["hit_status","status"], "")
        hit_status = raw_hit_status.lower().strip()
        gene = first(r, ["gene","amr_gene","resistance_gene","gene_or_status"], "NA")
        mutation = first(r, ["mutation","variant"], "NA")
        interp = first(r, ["interpretation","prediction","message"], "NA")
        evidence = first(r, ["evidence_level","evidence"], "NA")
        effect = first(r, ["effect"], "NA")

        # Rescue exact markers from interpretation when earlier AMR summaries
        # carried the marker there but left gene_or_status/mutation generic.
        gene, mutation, effect = rescue_marker_from_interpretation(gene, mutation, effect, interp)

        blob = " ".join([hit_status, gene, mutation, effect, evidence, interp]).lower()
        failed_terms = ["failed", "scanner_failed", "exit_1", "exit_127", "unknown option", "not available", "summary generation failed", "positive_control_failed", "dependency_missing"]
        nohit_terms = ["no marker", "no hit", "not detected", "none", "no curated marker", "no amr mutations found"]

        is_positive_control_failure = ("positive_control_failed" in blob or "positive control failed" in blob)
        is_failure = (is_positive_control_failure or hit_status in {"scanner_failed","failed","error"} or any(t in blob for t in failed_terms))

        curated_evidence = (
            "fungamr curated marker detected" in evidence.lower() or
            "curated chroquetas/fungamr marker detected" in interp.lower() or
            "fungamr mutation" in effect.lower()
        )
        exact_gene_mutation = (
            gene.upper() not in generic_nohit_values and
            "NO MARKER" not in gene.upper() and
            mutation.upper() not in {"NA","N/A","","SEE RAW TSV"}
        )

        is_hit = (not is_failure and (
            curated_evidence or
            hit_status in {"hit","detected","positive","reported","marker_reported"} or
            exact_gene_mutation
        ))

        # Do not let the historical default no_hit/status-less summaries override
        # a curated marker. Only call no-hit after excluding curated/exact hits.
        is_nohit = (not is_hit and not is_failure and (
            hit_status in {"no_hit","negative"} or any(t in blob for t in nohit_terms)
        ))

        display_status = "hit" if is_hit else ("positive_control_failed" if is_positive_control_failure else ("scanner_failed" if is_failure else "no_hit"))
        rows.append({
            "sample_id": sample,
            "species": first(r, ["species"], "NA"),
            "drug_class": first(r, ["drug_class","class"], "NA"),
            "drug": first(r, ["drug","antifungal"], "NA"),
            "gene": gene if (is_hit or is_failure) else "No curated marker detected",
            "mutation": mutation if is_hit else "NA",
            "effect": effect if is_hit else "NA",
            "evidence_level": evidence,
            "interpretation": interp,
            "hit_status": display_status
        })
    by_sample={}
    for r in rows:
        by_sample.setdefault(r["sample_id"], []).append(r)
    for s in samples:
        if s not in by_sample:
            rows.append({"sample_id":s,"species":"NA","drug_class":"NA","drug":"NA","gene":"No hit","mutation":"NA","effect":"NA","evidence_level":"NA","interpretation":"No AMR summary file found.","hit_status":"no_hit"})
    return rows

species_rows=parse_species()
assembly_rows=parse_assembly()
quast_rows=parse_quast()
busco_rows=parse_busco()
amr_rows=parse_amr()

species_by={r["sample_id"]:r for r in species_rows}
assembly_by={r["sample_id"]:r for r in assembly_rows}
quast_by={r.get("sample_id",""):r for r in quast_rows}
busco_by={r["sample_id"]:r for r in busco_rows}

def has_amr_hit(r):
    return r.get("hit_status") == "hit"

amr_by_count=Counter()
amr_fail_by_count=Counter()
for r in amr_rows:
    if has_amr_hit(r):
        amr_by_count[r["sample_id"]] += 1
    if r.get("hit_status") in {"scanner_failed", "positive_control_failed"}:
        amr_fail_by_count[r["sample_id"]] += 1

def badge(label, kind="info"):
    return f'<span class="badge {kind}">{label}</span>'

def species_badge(sp):
    if not sp or sp == "Not determined":
        return badge("Not determined", "muted")
    return badge(f"<i>{esc(sp)}</i>", "species")

def nohit_badge():
    return badge("No curated genomic AMR marker detected — not susceptible", "muted")

def amr_card_badge(s):
    if amr_by_count[s] > 0:
        return badge(str(amr_by_count[s])+" AMR hit(s)", "warn")
    if amr_fail_by_count[s] > 0:
        return badge("AMR validation/scanner warning", "warn")
    return nohit_badge()

def bar(value):
    f=safe_float(value)
    if f is None:
        return "NA"
    width=max(0,min(100,f))
    return f'<div class="bar"><span style="width:{width:.1f}%"></span></div><small>{f:.2f}%</small>'

def table(rows, cols, headers=None, fmt=None):
    headers=headers or {}
    fmt=fmt or {}
    out='<div class="table-wrap"><table><thead><tr>'
    out+=''.join(f'<th>{esc(headers.get(c,c))}</th>' for c in cols)
    out+='</tr></thead><tbody>'
    for r in rows:
        out+='<tr>'
        for c in cols:
            v=r.get(c,"NA")
            out += f'<td>{fmt[c](v,r) if c in fmt else esc(v)}</td>'
        out+='</tr>'
    out+='</tbody></table></div>'
    return out

def sample_cards():
    cards=[]
    for s in samples:
        sp=species_by.get(s,{})
        asm=assembly_by.get(s,{})
        qst=quast_by.get(s,{})
        conf_label, conf_kind, conf_note = species_confidence(s)
        overall, overall_kind, overall_note = surveillance_status(s)
        phy_label, phy_kind, phy_note = phylogeny_status(s)
        cards.append(f"""
        <details class="sample-card" open>
          <summary class="sample-head"><h3>{esc(s)}</h3><span>{species_badge(sp.get("species","Not determined"))} {badge(overall, overall_kind)}</span></summary>
          <div class="mini-grid">
            <div><small>Species reads</small><strong>{esc(sp.get("percent_reads","NA"))}%</strong></div>
            <div><small>Contigs</small><strong>{esc(asm.get("contigs", qst.get("# contigs","NA")))}</strong></div>
            <div><small>N50</small><strong>{esc(asm.get("n50", qst.get("N50","NA")))}</strong></div>
            <div><small>GC (%)</small><strong>{esc(qst.get("GC (%)","NA"))}</strong></div>
          </div>
          <p>{amr_card_badge(s)} {badge('Species confidence: ' + conf_label, conf_kind)} {badge('Phylogeny: ' + phy_label, phy_kind)}</p>
          <div class="details-note"><strong>Interpretation:</strong> {esc(overall_note)}. <strong>Species note:</strong> {esc(conf_note)}. <strong>Phylogeny note:</strong> {esc(phy_note)}</div>
        </details>""")
    return ''.join(cards)

top_species = "Not determined"
species_counts=Counter(r["species"] for r in species_rows)
if species_counts:
    top_species=species_counts.most_common(1)[0][0]

species_dist=''.join(
    f'<div class="dist-row">{species_badge(sp)}<div class="bar"><span style="width:{(count/max(species_counts.values()))*100 if species_counts else 0:.1f}%"></span></div><strong>{count}</strong></div>'
    for sp,count in species_counts.most_common()
)

n50_vals=[safe_float(r.get("n50")) for r in assembly_rows if safe_float(r.get("n50")) is not None]
median_n50=f"{statistics.median(n50_vals):,.0f}" if n50_vals else "NA"
total_hits=sum(1 for r in amr_rows if has_amr_hit(r))

summary_cols=["sample_id","species","percent_reads","contigs","n50","gc_percent","amr_hits"]
summary_rows=[]
for s in samples:
    summary_rows.append({
        "sample_id":s,
        "species":species_by.get(s,{}).get("species","Not determined"),
        "percent_reads":species_by.get(s,{}).get("percent_reads","NA"),
        "contigs":assembly_by.get(s,{}).get("contigs","NA"),
        "n50":assembly_by.get(s,{}).get("n50","NA"),
        "gc_percent":quast_by.get(s,{}).get("GC (%)","NA"),
        "amr_hits":str(amr_by_count[s])
    })

with open("rMAP-Myc-Candida-Candida_summary.tsv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=summary_cols, delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

species_table=table(species_rows, ["sample_id","species","percent_reads","clade_reads","taxon_reads","taxid","evidence"], {"sample_id":"Sample","species":"Top species","percent_reads":"Reads (%)"}, {"species":lambda v,r: species_badge(v), "percent_reads":lambda v,r: bar(v)})
assembly_table=table(assembly_rows, ["sample_id","contigs","total_bp","n50","largest_contig"], {"sample_id":"Sample","total_bp":"Total bp","n50":"N50","largest_contig":"Largest contig"})
quast_table=table(quast_rows, ["sample_id","# contigs","Largest contig","Total length","GC (%)","N50"], {"sample_id":"Sample"})
busco_table=table(busco_rows, ["sample_id","Complete_BUSCO_%","Single_copy_BUSCO_%","Duplicated_BUSCO_%","Fragmented_BUSCO_%","Missing_BUSCO_%","BUSCO_n","busco_status","busco_details"], {"sample_id":"Sample","BUSCO_n":"Markers","busco_status":"Status","busco_details":"Compleasm note"}, {"Complete_BUSCO_%":lambda v,r: bar(v), "Missing_BUSCO_%":lambda v,r: bar(v), "busco_status":lambda v,r: badge(v, "success" if str(v).upper() in ["PASS","PARSED"] else "warn")})
amr_table=table(amr_rows, ["sample_id","species","drug_class","drug","gene","mutation","effect","evidence_level","interpretation"], {"sample_id":"Sample","gene":"Gene / status"}, {"gene":lambda v,r: badge("No marker — not susceptible","muted") if r.get("hit_status")=="no_hit" else (badge("Positive control failed","warn") if r.get("hit_status")=="positive_control_failed" else (badge("Scanner failed","muted") if r.get("hit_status")=="scanner_failed" else badge(esc(v),"amrhit")))})

# -------------------------------------------------------------------------
# rc173 surveillance-report additions:
# 1) optional surveillance metadata
# 2) integrated readiness dashboard
# 3) species-confidence/mixed-sample warnings
# 4) phylogeny eligibility summary
# 5) pairwise SNP distance table from core-SNP alignments
# -------------------------------------------------------------------------

def parse_metadata():
    p = root / "surveillance_metadata.tsv"
    if not p.exists() or p.stat().st_size == 0:
        return [], {}
    rows=[]
    try:
        with p.open(errors="ignore", newline="") as f:
            reader=csv.DictReader(f, delimiter="\t")
            for row in reader:
                sid = first(row, ["sample_id","sample","sample_name","isolate_id","isolate"], "")
                if sid:
                    row["sample_id"] = sid
                rows.append(row)
    except Exception:
        rows=[]
    return rows, {r.get("sample_id",""):r for r in rows if r.get("sample_id")}

metadata_rows, metadata_by = parse_metadata()

def species_confidence(sample):
    pct = safe_float(species_by.get(sample,{}).get("percent_reads"))
    if pct is None:
        return "Not determined", "muted", "No species abundance value found"
    if pct >= 95:
        return "High", "success", "Top species abundance ≥95%"
    if pct >= 80:
        return "Moderate", "warn", "Top species abundance 80–95%; review for possible mixed signal"
    return "Low / review", "warn", "Top species abundance <80%; possible mixed sample, contamination, or weak species assignment"

def assembly_status(sample):
    asm = assembly_by.get(sample,{})
    n50 = safe_float(asm.get("n50"))
    contigs = safe_float(str(asm.get("contigs","NA")).replace(",",""))
    if n50 is None and contigs is None:
        return "Not run", "muted", "No assembly summary found"
    warnings=[]
    if n50 is not None and n50 < 20000:
        warnings.append("low N50")
    if contigs is not None and contigs > 2500:
        warnings.append("high contig count")
    if warnings:
        return "Review", "warn", "; ".join(warnings)
    return "Pass", "success", "Assembly continuity acceptable for surveillance screening"

def quast_contiguity_status(sample):
    label, kind, note = assembly_status(sample)
    if label == "Pass":
        return "Pass", "success", "QUAST contiguity metrics are acceptable for surveillance screening"
    if label == "Review":
        return "Review", "warn", "QUAST contiguity metrics suggest the assembly should be reviewed"
    return "Not available", "muted", "QUAST/assembly metrics were not available"

def completeness_status(sample):
    # rc175 QUAST-only mode: no gene-completeness tool is run.
    return "Not assessed", "muted", "Gene completeness was not assessed in QUAST-only mode"

def parse_phylogeny_summary():
    rows = []
    for p in sorted(root.glob("*phylogeny_group_summary.tsv")):
        try:
            with p.open(errors="ignore", newline="") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    rows.append(row)
        except Exception:
            pass
    return rows

def phylogeny_by_species():
    out={}
    for row in parse_phylogeny_summary():
        sp = str(row.get("species","")).strip()
        if sp:
            out[sp] = row
    return out

phylo_by_sp = phylogeny_by_species()

def phylogeny_status(sample):
    sp = species_by.get(sample,{}).get("species","Not determined")
    row = phylo_by_sp.get(sp)
    if not row:
        return "Not assessed", "muted", "No phylogeny group summary for this species"
    status = str(row.get("status","NA"))
    if status.upper() == "PASS":
        return "Included / eligible", "success", f"Species group passed phylogeny stage ({row.get('sample_count','NA')} samples)"
    if status.startswith("SKIPPED"):
        return "Excluded", "warn", row.get("notes", status)
    return "Review", "warn", row.get("notes", status)

def surveillance_status(sample):
    sp_label, sp_kind, sp_note = species_confidence(sample)
    asm_label, asm_kind, asm_note = assembly_status(sample)
    phy_label, phy_kind, phy_note = phylogeny_status(sample)
    hits = amr_by_count[sample]
    failures = amr_fail_by_count[sample]
    reasons=[]
    if sp_kind == "warn" or sp_label == "Not determined": reasons.append("species confidence")
    if asm_kind == "warn" or asm_label == "Not run": reasons.append("QUAST assembly QC")
    if failures > 0: reasons.append("AMR scanner")
    if hits > 0: reasons.append("AMR marker detected")
    if not reasons:
        return "Ready for surveillance interpretation", "success", "No major automated warning flags using species typing, QUAST assembly QC, AMR screening, and phylogeny status"
    return "Review", "warn", "; ".join(reasons)

def metadata_value(sample, keys):
    row = metadata_by.get(sample,{})
    return first(row, keys, "NA") if row else "NA"

def build_surveillance_rows():
    rows=[]
    for s in samples:
        sp_label, sp_kind, sp_note = species_confidence(s)
        asm_label, asm_kind, asm_note = assembly_status(s)
        quast_label, quast_kind, quast_note = quast_contiguity_status(s)
        phy_label, phy_kind, phy_note = phylogeny_status(s)
        overall, overall_kind, overall_note = surveillance_status(s)
        rows.append({
            "sample_id": s,
            "collection_date": metadata_value(s, ["collection_date","date","sampling_date"]),
            "site": metadata_value(s, ["site","facility","ward_or_facility","location"]),
            "species": species_by.get(s,{}).get("species","Not determined"),
            "species_confidence": sp_label,
            "assembly_qc": asm_label,
            "quast_contiguity_qc": quast_label,
            "amr_screen": f"{amr_by_count[s]} marker(s)" if amr_by_count[s] else ("AMR scanner warning" if amr_fail_by_count[s] else "No curated marker detected"),
            "phylogeny_status": phy_label,
            "surveillance_status": overall,
            "interpretation_note": overall_note
        })
    return rows

surveillance_rows = build_surveillance_rows()
with open("rMAP_Candida_surveillance_summary.tsv","w",newline="") as f:
    fields=["sample_id","collection_date","site","species","species_confidence","assembly_qc","quast_contiguity_qc","amr_screen","phylogeny_status","surveillance_status","interpretation_note"]
    w=csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(surveillance_rows)

surveillance_table = table(surveillance_rows,
    ["sample_id","collection_date","site","species","species_confidence","assembly_qc","quast_contiguity_qc","amr_screen","phylogeny_status","surveillance_status","interpretation_note"],
    {"sample_id":"Sample","collection_date":"Collection date","species_confidence":"Species confidence","assembly_qc":"Assembly QC","quast_contiguity_qc":"QUAST contiguity QC","amr_screen":"AMR screen","phylogeny_status":"Phylogeny","surveillance_status":"Surveillance status","interpretation_note":"Reason / note"},
    {
        "species": lambda v,r: species_badge(v),
        "species_confidence": lambda v,r: badge(v, "success" if v=="High" else ("warn" if "Low" in v or v=="Moderate" else "muted")),
        "assembly_qc": lambda v,r: badge(v, "success" if v=="Pass" else ("warn" if "Review" in v or "Fail" in v else "muted")),
        "quast_contiguity_qc": lambda v,r: badge(v, "success" if v=="Pass" else ("warn" if "Review" in v or "Fail" in v else "muted")),
        "amr_screen": lambda v,r: badge(v, "amrhit" if "marker" in str(v) and not str(v).startswith("No") else ("warn" if "warning" in str(v).lower() else "muted")),
        "phylogeny_status": lambda v,r: badge(v, "success" if str(v).startswith("Included") else ("warn" if v in {"Excluded","Review"} else "muted")),
        "surveillance_status": lambda v,r: badge(v, "success" if str(v).startswith("Ready") else "warn")
    })

if metadata_rows:
    # Show a compact metadata table with the most surveillance-relevant columns available.
    preferred=["sample_id","country","site","collection_date","specimen_type","patient_group","ward_or_facility","sequencing_platform"]
    available=[]
    for c in preferred:
        if any(c in row for row in metadata_rows):
            available.append(c)
    if "sample_id" not in available:
        available.insert(0,"sample_id")
    metadata_table = table(metadata_rows, available, {"sample_id":"Sample"})
else:
    metadata_table = '<div class="note"><strong>No surveillance metadata TSV was provided.</strong> Add <code>rMAP_Candida.surveillance_metadata_tsv</code> to your input JSON to display collection date, country/site, specimen type, patient group, facility/ward, and sequencing platform.</div>'

def read_alignment(path):
    seqs={}; name=None; chunks=[]
    try:
        with open(path, errors="ignore") as fh:
            for line in fh:
                line=line.strip()
                if not line: continue
                if line.startswith(">"):
                    if name is not None: seqs[name]="".join(chunks).upper()
                    name=line[1:].split()[0]; chunks=[]
                else:
                    chunks.append(line)
            if name is not None: seqs[name]="".join(chunks).upper()
    except Exception:
        return {}
    return seqs

def pairwise_distances_from_alignment(path):
    seqs=read_alignment(path)
    names=list(seqs)
    rows=[]
    valid=set("ACGT")
    label=path.name.replace(".core_snps.aln","").replace("_"," ")
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            a,b=names[i],names[j]
            sa,sb=seqs[a],seqs[b]
            L=min(len(sa),len(sb))
            compared=0; dist=0
            for k in range(L):
                ca,cb=sa[k],sb[k]
                if ca in valid and cb in valid:
                    compared += 1
                    if ca != cb: dist += 1
            rows.append({"species":label,"sample_a":a,"sample_b":b,"snp_distance":str(dist),"compared_sites":str(compared)})
    return rows

def build_snp_distance_section():
    distance_rows=[]
    for p in sorted(root.glob("*.core_snps.aln")):
        distance_rows.extend(pairwise_distances_from_alignment(p))
    if not distance_rows:
        pathlib.Path("rMAP_Candida_pairwise_snp_distances.tsv").write_text("species\tsample_a\tsample_b\tsnp_distance\tcompared_sites\n")
        return '<p>No core-SNP alignment was available for pairwise SNP-distance calculation.</p>'
    with open("rMAP_Candida_pairwise_snp_distances.tsv","w",newline="") as f:
        fields=["species","sample_a","sample_b","snp_distance","compared_sites"]
        w=csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(distance_rows)
    nearest={}
    for r in distance_rows:
        try: d=int(r["snp_distance"])
        except Exception: continue
        for a,b in [(r["sample_a"],r["sample_b"]),(r["sample_b"],r["sample_a"] )]:
            if a not in nearest or d < nearest[a]["snp_distance_num"]:
                nearest[a] = {"sample_id":a,"closest_sample":b,"snp_distance":str(d),"snp_distance_num":d,"species":r["species"],"cluster_flag":"Possible close genetic relationship" if d <= 25 else "Distinct / review with metadata"}
    nearest_rows=[{k:v for k,v in row.items() if k!="snp_distance_num"} for row in nearest.values()]
    nearest_rows=sorted(nearest_rows, key=lambda x:(x.get("species",""), x.get("sample_id","")))
    nearest_table=table(nearest_rows, ["species","sample_id","closest_sample","snp_distance","cluster_flag"], {"sample_id":"Sample","closest_sample":"Closest sample","snp_distance":"SNP distance","cluster_flag":"Conservative cluster flag"}, {"cluster_flag":lambda v,r: badge(v,"warn" if str(v).startswith("Possible") else "muted")})
    pair_table=table(distance_rows, ["species","sample_a","sample_b","snp_distance","compared_sites"], {"sample_a":"Sample A","sample_b":"Sample B","snp_distance":"SNP distance","compared_sites":"Compared core-SNP sites"})
    return '<h3>Closest-neighbor summary</h3>' + nearest_table + '<h3>Pairwise SNP-distance matrix/table</h3>' + pair_table + '<div class="note"><strong>SNP-distance interpretation:</strong> low SNP distances suggest close genetic relatedness but should not be interpreted as transmission without epidemiological metadata, recombination-aware analysis, and species-specific validation.</div>'

snp_distance_section = build_snp_distance_section()

def parse_phylogeny_summary():
    rows = []
    for p in sorted(root.glob("*phylogeny_group_summary.tsv")):
        try:
            with p.open(errors="ignore", newline="") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    rows.append(row)
        except Exception:
            pass
    return rows

def build_phylogeny_section():
    phylo_rows = parse_phylogeny_summary()
    image_paths = sorted(list(root.glob("*.core_snp_tree.png")) + list(root.glob("*.core_snp_tree.svg")))

    parts = []
    if phylo_rows:
        cols = ["species","status","sample_count","branch","ploidy","variable_sites","notes"]
        parts.append(table(phylo_rows, cols, {
            "species":"Species",
            "status":"Status",
            "sample_count":"Samples",
            "branch":"Variant-calling branch",
            "ploidy":"Ploidy model",
            "variable_sites":"Core variable sites",
            "notes":"Notes"
        }, {
            "species":lambda v,r: species_badge(v),
            "status":lambda v,r: badge(v, "success" if str(v).upper() in {"BUILT","TREE_READY","PASS"} else "warn")
        }))

    if image_paths:
        for img in image_paths:
            label = img.name.replace(".core_snp_tree.png","").replace(".core_snp_tree.svg","").replace("_"," ")
            try:
                if img.suffix.lower() == ".svg":
                    data = img.read_text(errors="replace")
                    parts.append(f'<h3>{esc(label)}</h3><div class="tree-panel">{data}</div>')
                else:
                    encoded = base64.b64encode(img.read_bytes()).decode("ascii")
                    parts.append(f'<h3>{esc(label)}</h3><div class="tree-panel"><img class="tree-img" src="data:image/png;base64,{encoded}" alt="{esc(label)} core-SNP phylogenetic tree"></div>')
            except Exception as exc:
                parts.append(f'<p>Could not embed tree image {esc(img.name)}: {esc(exc)}</p>')

    skipped_rows = [r for r in phylo_rows if str(r.get("status","")).upper() not in {"BUILT","TREE_READY","PASS"}]
    if skipped_rows:
        skipped_bits = []
        for r in skipped_rows:
            sp = r.get("species", "Unknown species")
            st = r.get("status", "Not built")
            sc = r.get("sample_count", "NA")
            br = r.get("branch", "NA")
            nt = r.get("notes", "")
            skipped_bits.append(
                f'<li><strong>{esc(sp)}</strong>: {esc(st)}; samples after species-specific consensus/eligibility = {esc(sc)}; branch = {esc(br)}; note = {esc(nt)}</li>'
            )
        parts.append('<div class="note"><strong>Species without a displayed tree:</strong><ul>' + ''.join(skipped_bits) + '</ul></div>')

    if not parts:
        return '<p>Phylogeny was not enabled or no eligible species group met the minimum sample/reference requirements.</p>'

    parts.append('<div class="note"><strong>Phylogeny interpretation note:</strong> rMAP-Candida builds one tree per species. Mixed-species Candida phylogenies are intentionally avoided. Treat these trees as species-level genomic relatedness visualizations unless recombination filtering and epidemiologic metadata support transmission interpretation.</div>')
    return ''.join(parts)

phylogeny_section = build_phylogeny_section()

run_dt = datetime.datetime.utcnow()
run_generated_utc = run_dt.strftime("%Y-%m-%d %H:%M:%S UTC")
run_stamp = run_dt.strftime("%Y%m%d_%H%M%S_UTC")


css = """
body{margin:0;background:#f3f6fb;color:#13242d;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Arial,sans-serif}
.hero{background:linear-gradient(125deg,#00524d,#146c8f,#2563eb);color:white;padding:44px 0 72px}
.wrap{max-width:1180px;margin:0 auto;padding:0 24px}
.kicker{font-size:12px;letter-spacing:.12em;text-transform:uppercase;font-weight:800}
h1{font-size:38px;margin:10px 0 14px}
.hero p{max-width:760px;line-height:1.55}
.pill{display:inline-block;padding:8px 12px;border-radius:999px;background:#dbeafe;color:#1d4ed8;font-weight:800;font-size:12px;margin-right:8px}
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:18px;margin-top:24px;margin-bottom:24px}
.metric{background:white;border-radius:18px;padding:22px;box-shadow:0 12px 30px #0f172a18}
.metric small{display:block;text-transform:uppercase;color:#64748b;letter-spacing:.1em;font-weight:800}
.metric strong{display:block;font-size:30px;margin-top:8px}
.card{background:white;border-radius:18px;padding:24px;margin:22px 0;box-shadow:0 8px 25px #0f172a12;border:1px solid #e5e7eb}
.two{display:grid;grid-template-columns:1.25fr .85fr;gap:24px}
.sample-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:18px}
.sample-card{border:1px solid #e5e7eb;border-radius:18px;padding:20px;background:#fff}.sample-card summary{cursor:pointer;list-style:none}.sample-card summary::-webkit-details-marker{display:none}.details-note{background:#f8fafc;border:1px solid #e5e7eb;border-radius:12px;padding:12px;margin-top:12px;color:#334155}.toc{display:flex;flex-wrap:wrap;gap:12px;align-items:center;background:white;border:1px solid #e5e7eb;border-radius:16px;padding:14px 16px;margin:24px 0 0;box-shadow:0 8px 25px #0f172a12}.toc a{display:inline-block;margin:6px 10px 6px 0;font-weight:800;color:#1d4ed8;text-decoration:none}.downloads{display:grid;grid-template-columns:repeat(3,1fr);gap:10px}.download-pill{background:#f8fafc;border:1px solid #e5e7eb;border-radius:12px;padding:12px;font-weight:800;color:#334155}
.sample-head{display:flex;justify-content:space-between;gap:12px;align-items:center}
.sample-head h3{margin:0}
.mini-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:16px 0}
.mini-grid div{background:#f8fafc;border:1px solid #e5e7eb;border-radius:12px;padding:12px}
.mini-grid small{display:block;text-transform:uppercase;color:#64748b;font-size:11px;font-weight:800}
.mini-grid strong{display:block;margin-top:4px}
.badge{display:inline-block;border-radius:999px;padding:6px 10px;font-weight:800;font-size:12px}
.badge.species{background:#ccfbf1;color:#0f766e}.badge.success{background:#dcfce7;color:#166534}.badge.warn{background:#ffedd5;color:#9a3412}.badge.amrhit{background:#8B0000;color:#ffffff}.badge.muted{background:#e5e7eb;color:#475569}.badge.info{background:#dbeafe;color:#1d4ed8}
.table-wrap{overflow:auto;border:1px solid #e5e7eb;border-radius:14px}
table{width:100%;border-collapse:collapse;font-size:14px}th,td{padding:12px;border-bottom:1px solid #e5e7eb;text-align:left;vertical-align:top}th{background:#f8fafc;text-transform:uppercase;font-size:12px;letter-spacing:.06em}
.bar{height:10px;background:#e5e7eb;border-radius:999px;overflow:hidden;min-width:130px}.bar span{display:block;height:100%;background:linear-gradient(90deg,#0f766e,#2563eb)}
.dist-row{display:grid;grid-template-columns:180px 1fr 30px;gap:12px;align-items:center;margin:12px 0}
.note{border-left:5px solid #0f766e;background:#ecfdf5;border-radius:12px;padding:14px}
.tree-panel{overflow:auto;border:1px solid #e5e7eb;border-radius:16px;background:#fff;padding:18px 22px;margin:14px 0;min-height:180px}
.tree-panel h3{margin-top:0}.tree-img{max-width:100%;height:auto;display:block;margin:0 auto;image-rendering:auto}
.footer{text-align:center;color:#64748b;padding:40px 0}
@media(max-width:900px){.metrics,.sample-grid,.two{grid-template-columns:1fr}.mini-grid{grid-template-columns:repeat(2,1fr)}}
"""

html_doc=f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>rMAP-Myc-Candida report</title><style>{css}</style></head>
<body>
<section class="hero"><div class="wrap">
<div class="kicker">rMAP-Myc-Candida surveillance report</div>
<h1>Integrated Candida fungal genomics report</h1>
<p>This report summarizes paired-end fungal genome analysis for Candida-focused surveillance, combining read QC, Kraken2/Bracken species typing, MEGAHIT assembly, QUAST-only assembly contiguity assessment, and genomic antifungal-resistance screening.</p>
<span class="pill">Run generated: {esc(run_generated_utc)}</span><span class="pill">Run stamp: {esc(run_stamp)}</span><span class="pill">Custom Candida Kraken2/Bracken DB</span><span class="pill">MEGAHIT assembly</span>
</div></section>
<div class="wrap">
<nav class="toc"><a href="#executive">Executive summary</a><a href="#surveillance">Surveillance dashboard</a><a href="#metadata">Metadata</a><a href="#samples">Samples</a><a href="#species">Species</a><a href="#assembly">Assembly QC</a><a href="#amr">AMR</a><a href="#phylogeny">Phylogeny</a><a href="#snpdist">SNP distances</a><a href="#provenance">Outputs</a></nav>
<section class="metrics">
<div class="metric"><small>Samples analyzed</small><strong>{len(samples)}</strong></div>
<div class="metric"><small>Top species groups</small><strong>{len(species_counts)}</strong></div>
<div class="metric"><small>Total AMR hits</small><strong>{total_hits}</strong></div>
<div class="metric"><small>Median N50</small><strong>{esc(median_n50)}</strong></div>
</section>
<section class="card two" id="executive"><div><h2>1. Executive summary</h2>
<p>The primary detected fungal species group was {species_badge(top_species)}. The sample-level cards below provide a compact interpretation of species assignment, QUAST assembly-contiguity metrics, phylogeny status, and genomic AMR screening status.</p>
<div class="note"><strong>Interpretation note:</strong> genomic antifungal-resistance findings should be treated as screening evidence. Clinically important results should be interpreted with isolate metadata, species identity, validated mutation catalogues, and phenotypic antifungal susceptibility testing where required.</div></div>
<div><h3>Species distribution</h3>{species_dist or '<p>No species calls available.</p>'}</div></section>
<section class="card" id="surveillance"><h2>2. Surveillance readiness and interpretation dashboard</h2><p>This integrated table combines species confidence, QUAST assembly-contiguity metrics, AMR marker status, phylogeny eligibility, and optional metadata into a practical surveillance-readiness view.</p>{surveillance_table}</section>
<section class="card" id="metadata"><h2>3. Surveillance metadata</h2>{metadata_table}</section>
<section class="card" id="samples"><h2>4. Sample-level surveillance summary</h2><div class="sample-grid">{sample_cards()}</div></section>
<section class="card" id="species"><h2>5. Candida species typing using Kraken2/Bracken</h2><p>Top species calls are derived from the custom Candida-focused Kraken2/Bracken database bundled in the species-typing Docker image.</p>{species_table}</section>
<section class="card" id="assembly"><h2>6. MEGAHIT assembly summary</h2>{assembly_table}</section>
<section class="card"><h2>7. Assembly quality assessment with QUAST</h2><p>QUAST values are parsed from the native per-sample <code>report.tsv</code> format and transposed into one row per sample. This default surveillance mode reports QUAST assembly contiguity metrics such as contig count, total length, N50, largest contig, and GC percentage. It is very fast. BUSCO and Compleasm tasks are available as optional modules, but they are disabled by default in the recommended JSON and are not required for this QUAST-based assembly assessment.</p>{quast_table}</section>
<section class="card" id="amr"><h2>8. Fungal antifungal-resistance characterization</h2><p>This section reports mutation/gene-level evidence emitted by the configured fungal AMR container. A “No marker detected” result is not a susceptible call. Fluconazole resistance can be caused by ERG11 alterations, TAC1/UPC2/MRR1/PDR1-mediated efflux, aneuploidy/LOH, copy-number changes, species-specific mechanisms, or markers absent from the current AMR database.</p>{amr_table}</section>
<section class="card" id="phylogeny"><h2>9. Species-aware core-SNP phylogeny</h2><p>When enabled, rMAP-Candida builds phylogenies separately for each species with sufficient samples and a matching reference. Mixed-species phylogenies are intentionally avoided. Outputs include species-group summaries, core-SNP alignments, and IQ-TREE Newick trees.</p>{phylogeny_section if 'phylogeny_section' in globals() else '<p>Phylogeny was not enabled or no eligible species group met the minimum sample/reference requirements.</p>'}</section>
<section class="card" id="snpdist"><h2>10. Species-aware pairwise SNP distances and closest-neighbor summary</h2>{snp_distance_section}</section>
<section class="card" id="provenance"><h2>11. Output navigation and provenance</h2><p>The workflow emits downloadable tabular outputs in addition to this integrated HTML report.</p><div class="downloads"><div class="download-pill">rMAP-Myc-Candida-Candida_summary.tsv</div><div class="download-pill">rMAP_Candida_surveillance_summary.tsv</div><div class="download-pill">rMAP_Candida_pairwise_snp_distances.tsv</div></div><ul><li>Per-sample Kraken2, Bracken, FASTQ QC, assembly, QUAST, AMR, and optional phylogeny files are available in the Cromwell execution outputs.</li><li>Dockerized execution supports reproducibility across local Cromwell and cloud environments.</li></ul></section>
<div class="footer">rMAP-Myc-Candida | Rapid Mycological Analysis Pipeline for Candida</div>
</div></body></html>"""

pathlib.Path("rMAP_Candida_report.html").write_text(html_doc)
PY
  >>>

  output {
    File html_report = "rMAP_Candida_report.html"
    File summary_tsv = "rMAP-Myc-Candida-Candida_summary.tsv"
    File surveillance_summary_tsv = "rMAP_Candida_surveillance_summary.tsv"
    File pairwise_snp_distances_tsv = "rMAP_Candida_pairwise_snp_distances.tsv"
  }

  runtime {
    docker: "python:3.11-slim"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

