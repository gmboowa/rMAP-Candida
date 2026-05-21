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
    Boolean do_busco = true
    Boolean do_fungal_amr = true

    String fungal_kraken2_bracken_docker = "gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db"
    String fungamr_docker = "gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed"
    String megahit_docker = "quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6"
    String quast_docker = "staphb/quast:5.2.0"
    String busco_docker = "ezlabgva/busco:v5.7.1_cv1"
    String fastp_docker = "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    String fastqc_docker = "staphb/fastqc:0.12.1"
    String multiqc_docker = "multiqc/multiqc:v1.24"

    Int max_cpus = 8
    Int max_memory_gb = 32
    Int min_read_length = 50
    Int bracken_read_length = 150
    String busco_lineage = "saccharomycetes_odb10"
    String kraken_db_path = "/opt/kraken2_db/candida"
    String bracken_level = "S"
  }

  Int n = length(sample_names)
  Int cpu_4 = if max_cpus < 4 then max_cpus else 4
  Int cpu_8 = if max_cpus < 8 then max_cpus else 8

  #
  # Stage 1: per-sample read-level work and assembly.
  #
  # The first scatter deliberately contains only tasks that should run directly from
  # the paired FASTQ inputs.  Assembly-dependent tasks are launched in a second,
  # explicit scatter below over the complete Array[File] of assembly outputs.  This
  # prevents QUAST, BUSCO, and AMR from accidentally receiving only shard-0 or a
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
          cpu = cpu_4,
          memory_gb = max_memory_gb
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
          cpu = cpu_4,
          memory_gb = max_memory_gb
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
          cpu = cpu_8,
          memory_gb = max_memory_gb
      }
    }

    if (do_assembly) {
      call FUNGAL_ASSEMBLY as ASSEMBLY {
        input:
          sample_name = sample_names[i],
          read1 = analysis_read1,
          read2 = analysis_read2,
          docker_image = megahit_docker,
          cpu = cpu_8,
          memory_gb = max_memory_gb
      }
    }
  }

  #
  # Stage 2: per-assembly downstream work.
  #
  # This second scatter is intentionally outside the assembly scatter.  Cromwell
  # must first collect all ASSEMBLY.contigs_fasta outputs, then QUAST, BUSCO, and
  # AMR are scattered over the full array.  With two samples, this creates:
  #   call-QUAST/shard-0 and call-QUAST/shard-1
  #   call-BUSCO/shard-0 and call-BUSCO/shard-1
  #   call-AMR/shard-0   and call-AMR/shard-1
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
            cpu = cpu_4,
            memory_gb = max_memory_gb
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
            cpu = cpu_8,
            memory_gb = max_memory_gb
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
            threads = cpu_4
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
        docker_image = multiqc_docker
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
      busco_summaries = select_first([BUSCO.busco_short_summary_txt, []]),
      busco_summary_tsvs = select_first([BUSCO.busco_summary_tsv, []]),
      amr_summaries = select_first([AMR.amr_summary_tsv, []]),
      amr_htmls = select_first([AMR.amr_report_html, []])
  }

  output {
    File rmap_myc_html_report = MERGE_MYC_REPORTS.html_report
    File rmap_myc_summary_tsv = MERGE_MYC_REPORTS.summary_tsv
    Array[File] trimmed_reads_1 = select_all(TRIM.trimmed_read1)
    Array[File] trimmed_reads_2 = select_all(TRIM.trimmed_read2)
    File? multiqc_report = MULTIQC_REPORT.multiqc_report
    Array[File] fungal_species_summaries = select_all(SPECIES.top_species_tsv)
    Array[File] fungal_assemblies = select_all(ASSEMBLY.contigs_fasta)
    Array[File] fungal_amr_summaries = select_first([AMR.amr_summary_tsv, []])
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
  }
}

task MULTIQC_REPORT {
  input {
    Array[File] fastqc_reports
    Array[File] trimming_json
    Array[File] trimming_html
    String docker_image
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
    cpu: 2
    memory: "8 GB"
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
    docker: busco_docker
    cpu: cpu
    memory: memory_gb + " GB"
  }
}



task FUNGAL_AMR_CHARACTERIZATION {
  input {
    String sample_id
    File assembly_fasta
    File species_summary_tsv
    String fungal_amr_docker_image
    Int threads = 4
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
    # ChroQueTas requires miniprot. Older AMR images can contain ChroQueTas.sh but
    # lack miniprot, causing:
    #   ERROR: miniprot is required and not installed
    # Previous rc146 logic could then collapse this technical failure into
    # "No curated marker detected". This block tries a best-effort runtime install
    # and, if still missing, reports SCANNER_DEPENDENCY_MISSING instead of a no-hit.
    MINIPROT_MISSING=0
    if ! command -v miniprot >/dev/null 2>&1; then
      echo "WARNING: miniprot is missing before scanner execution; attempting runtime bootstrap." >> "${LOG_OUT}"
      if command -v micromamba >/dev/null 2>&1; then
        micromamba install -y -c bioconda -c conda-forge miniprot >> "${LOG_OUT}" 2>&1 || true
      elif command -v mamba >/dev/null 2>&1; then
        mamba install -y -c bioconda -c conda-forge miniprot >> "${LOG_OUT}" 2>&1 || true
      elif command -v conda >/dev/null 2>&1; then
        conda install -y -c bioconda -c conda-forge miniprot >> "${LOG_OUT}" 2>&1 || true
      else
        echo "No micromamba/mamba/conda executable found for runtime miniprot bootstrap." >> "${LOG_OUT}"
      fi
      hash -r 2>/dev/null || true
    fi
    if command -v miniprot >/dev/null 2>&1; then
      echo "Confirmed dependency: miniprot at $(command -v miniprot)" >> "${LOG_OUT}"
      miniprot --version >> "${LOG_OUT}" 2>&1 || true
    else
      MINIPROT_MISSING=1
      echo "DEPENDENCY_MISSING: miniprot is not available after bootstrap attempt." >> "${LOG_OUT}"
    fi

    # -------------------------------------------------------------------------
    # v5 FIX: use a ChroQueTas-enabled AMR Docker image and keep shell-only parsing.
    # The previous WDL patch protected Cromwell outputs but the old Docker image
    # did not contain ChroQueTas in PATH. This WDL defaults to the rebuilt image
    # gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v6. Parsing remains POSIX shell/awk.
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
      echo "DEPENDENCY_MISSING: miniprot is required by ChroQueTas but is not installed in the AMR Docker image and runtime bootstrap failed. Build/pull an AMR image that includes miniprot, then rerun." > amr_out/scanner.stderr
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
    docker: fungal_amr_docker_image
    cpu: threads
    memory: "8 GB"
    disks: "local-disk 40 HDD"
    continueOnReturnCode: 0
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
    Array[File] busco_summaries
    Array[File] busco_summary_tsvs
    Array[File] amr_summaries
    Array[File] amr_htmls
  }

  command <<<
    set -euo pipefail

    mkdir -p report_inputs report_links
    cp ~{sep=' ' species_top_tsvs} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' kraken_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' bracken_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' assembly_summaries} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' quast_reports} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' busco_summaries} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' busco_summary_tsvs} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' amr_summaries} report_inputs/ 2>/dev/null || true
    cp ~{sep=' ' trimming_htmls} report_links/ 2>/dev/null || true
    cp ~{sep=' ' amr_htmls} report_links/ 2>/dev/null || true

    rm -f sample_names.txt
    for s in ~{sep=' ' sample_names}; do
      echo "$s" >> sample_names.txt
    done

    python3 <<'PY'
import csv, html, pathlib, re, datetime, statistics
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
    Parse normalized BUSCO TSVs emitted by the BUSCO task.

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
                "busco_details":"No normalized BUSCO TSV was found."
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
                    "busco_details":"Parsed from native BUSCO short-summary text."
                })
            by[sample]=row

    return [by.get(s,{"sample_id":s,"Complete_BUSCO_%":"NA","Single_copy_BUSCO_%":"NA","Duplicated_BUSCO_%":"NA","Fragmented_BUSCO_%":"NA","Missing_BUSCO_%":"NA","BUSCO_n":"NA","busco_status":"missing","busco_details":"No BUSCO output was found for this sample."}) for s in samples]

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
        bus=busco_by.get(s,{})
        cards.append(f"""
        <div class="sample-card">
          <div class="sample-head"><h3>{esc(s)}</h3>{species_badge(sp.get("species","Not determined"))}</div>
          <div class="mini-grid">
            <div><small>Species reads</small><strong>{esc(sp.get("percent_reads","NA"))}%</strong></div>
            <div><small>Contigs</small><strong>{esc(asm.get("contigs","NA"))}</strong></div>
            <div><small>N50</small><strong>{esc(asm.get("n50","NA"))}</strong></div>
            <div><small>BUSCO complete</small><strong>{esc(bus.get("Complete_BUSCO_%","NA"))}%</strong></div>
          </div>
          <p>{amr_card_badge(s)}</p>
        </div>""")
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

summary_cols=["sample_id","species","percent_reads","contigs","n50","Complete_BUSCO_%","amr_hits"]
summary_rows=[]
for s in samples:
    summary_rows.append({
        "sample_id":s,
        "species":species_by.get(s,{}).get("species","Not determined"),
        "percent_reads":species_by.get(s,{}).get("percent_reads","NA"),
        "contigs":assembly_by.get(s,{}).get("contigs","NA"),
        "n50":assembly_by.get(s,{}).get("n50","NA"),
        "Complete_BUSCO_%":busco_by.get(s,{}).get("Complete_BUSCO_%","NA"),
        "amr_hits":str(amr_by_count[s])
    })

with open("rMAP-Myc-Candida-Candida_summary.tsv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=summary_cols, delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

species_table=table(species_rows, ["sample_id","species","percent_reads","clade_reads","taxon_reads","taxid","evidence"], {"sample_id":"Sample","species":"Top species","percent_reads":"Reads (%)"}, {"species":lambda v,r: species_badge(v), "percent_reads":lambda v,r: bar(v)})
assembly_table=table(assembly_rows, ["sample_id","contigs","total_bp","n50","largest_contig"], {"sample_id":"Sample","total_bp":"Total bp","n50":"N50","largest_contig":"Largest contig"})
quast_table=table(quast_rows, ["sample_id","# contigs","Largest contig","Total length","GC (%)","N50"], {"sample_id":"Sample"})
busco_table=table(busco_rows, ["sample_id","Complete_BUSCO_%","Single_copy_BUSCO_%","Duplicated_BUSCO_%","Fragmented_BUSCO_%","Missing_BUSCO_%","BUSCO_n","busco_status","busco_details"], {"sample_id":"Sample","BUSCO_n":"BUSCO markers","busco_status":"Status","busco_details":"BUSCO note"}, {"Complete_BUSCO_%":lambda v,r: bar(v), "Missing_BUSCO_%":lambda v,r: bar(v), "busco_status":lambda v,r: badge(v, "success" if str(v).upper() in ["PASS","PARSED"] else "warn")})
amr_table=table(amr_rows, ["sample_id","species","drug_class","drug","gene","mutation","effect","evidence_level","interpretation"], {"sample_id":"Sample","gene":"Gene / status"}, {"gene":lambda v,r: badge("No marker — not susceptible","muted") if r.get("hit_status")=="no_hit" else (badge("Positive control failed","warn") if r.get("hit_status")=="positive_control_failed" else (badge("Scanner failed","muted") if r.get("hit_status")=="scanner_failed" else badge(esc(v),"amrhit")))})

now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

css = """
body{margin:0;background:#f3f6fb;color:#13242d;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Arial,sans-serif}
.hero{background:linear-gradient(125deg,#00524d,#146c8f,#2563eb);color:white;padding:44px 0 108px}
.wrap{max-width:1180px;margin:0 auto;padding:0 24px}
.kicker{font-size:12px;letter-spacing:.12em;text-transform:uppercase;font-weight:800}
h1{font-size:38px;margin:10px 0 14px}
.hero p{max-width:760px;line-height:1.55}
.pill{display:inline-block;padding:8px 12px;border-radius:999px;background:#dbeafe;color:#1d4ed8;font-weight:800;font-size:12px;margin-right:8px}
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:18px;margin-top:-52px}
.metric{background:white;border-radius:18px;padding:22px;box-shadow:0 12px 30px #0f172a18}
.metric small{display:block;text-transform:uppercase;color:#64748b;letter-spacing:.1em;font-weight:800}
.metric strong{display:block;font-size:30px;margin-top:8px}
.card{background:white;border-radius:18px;padding:24px;margin:22px 0;box-shadow:0 8px 25px #0f172a12;border:1px solid #e5e7eb}
.two{display:grid;grid-template-columns:1.25fr .85fr;gap:24px}
.sample-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:18px}
.sample-card{border:1px solid #e5e7eb;border-radius:18px;padding:20px;background:#fff}
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
.footer{text-align:center;color:#64748b;padding:40px 0}
@media(max-width:900px){.metrics,.sample-grid,.two{grid-template-columns:1fr}.mini-grid{grid-template-columns:repeat(2,1fr)}}
"""

html_doc=f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>rMAP-Myc-Candida report</title><style>{css}</style></head>
<body>
<section class="hero"><div class="wrap">
<div class="kicker">rMAP-Myc-Candida surveillance report</div>
<h1>Integrated Candida fungal genomics report</h1>
<p>This report summarizes paired-end fungal genome analysis for Candida-focused surveillance, combining read QC, Kraken2/Bracken species typing, MEGAHIT assembly, assembly quality assessment, BUSCO completeness, and genomic antifungal-resistance screening.</p>
<span class="pill">Generated {esc(now)}</span><span class="pill">Custom Candida Kraken2/Bracken DB</span><span class="pill">MEGAHIT assembly</span>
</div></section>
<div class="wrap">
<section class="metrics">
<div class="metric"><small>Samples analyzed</small><strong>{len(samples)}</strong></div>
<div class="metric"><small>Top species groups</small><strong>{len(species_counts)}</strong></div>
<div class="metric"><small>Total AMR hits</small><strong>{total_hits}</strong></div>
<div class="metric"><small>Median N50</small><strong>{esc(median_n50)}</strong></div>
</section>
<section class="card two"><div><h2>1. Executive summary</h2>
<p>The primary detected fungal species group was {species_badge(top_species)}. The sample-level cards below provide a compact interpretation of species assignment, assembly continuity, BUSCO completeness, and genomic AMR screening status.</p>
<div class="note"><strong>Interpretation note:</strong> genomic antifungal-resistance findings should be treated as screening evidence. Clinically important results should be interpreted with isolate metadata, species identity, validated mutation catalogues, and phenotypic antifungal susceptibility testing where required.</div></div>
<div><h3>Species distribution</h3>{species_dist or '<p>No species calls available.</p>'}</div></section>
<section class="card"><h2>2. Sample-level surveillance summary</h2><div class="sample-grid">{sample_cards()}</div></section>
<section class="card"><h2>3. Candida species typing using Kraken2/Bracken</h2><p>Top species calls are derived from the custom Candida-focused Kraken2/Bracken database bundled in the species-typing Docker image.</p>{species_table}</section>
<section class="card"><h2>4. MEGAHIT assembly summary</h2>{assembly_table}</section>
<section class="card"><h2>5. Assembly quality assessment with QUAST</h2><p>QUAST values are parsed from the native per-sample <code>report.tsv</code> format and transposed into one row per sample.</p>{quast_table}</section>
<section class="card"><h2>6. BUSCO fungal completeness</h2><p>BUSCO is run per assembly. For Candida, this patched workflow defaults to saccharomycetes_odb10 and retries ascomycota_odb10/fungi_odb10 if needed. If values remain NA, inspect the BUSCO_STATUS and BUSCO_NOTE columns because the task now reports the real reason instead of silently masking failures.</p>{busco_table}</section>
<section class="card"><h2>7. Fungal antifungal-resistance characterization</h2><p>This section reports mutation/gene-level evidence emitted by the configured fungal AMR container. A “No marker detected” result is not a susceptible call. Fluconazole resistance can be caused by ERG11 alterations, TAC1/UPC2/MRR1/PDR1-mediated efflux, aneuploidy/LOH, copy-number changes, species-specific mechanisms, or markers absent from the current AMR database.</p>{amr_table}</section>
<section class="card"><h2>8. Output navigation and provenance</h2><ul><li>Per-sample Kraken2, Bracken, FASTQ QC, assembly, QUAST, BUSCO, and AMR files are available in the Cromwell execution outputs.</li><li>The workflow also emits <code>rMAP-Myc-Candida-Candida_summary.tsv</code> for downstream tabular review.</li><li>Dockerized execution supports reproducibility across local Cromwell and cloud environments.</li></ul></section>
<div class="footer">rMAP-Myc-Candida | Rapid Mycological Analysis Pipeline for Candida</div>
</div></body></html>"""

pathlib.Path("rMAP_Candida_report.html").write_text(html_doc)
PY
  >>>

  output {
    File html_report = "rMAP_Candida_report.html"
    File summary_tsv = "rMAP-Myc-Candida-Candida_summary.tsv"
  }

  runtime {
    docker: "python:3.11-slim"
    cpu: 1
    memory: "4 GB"
  }
}

