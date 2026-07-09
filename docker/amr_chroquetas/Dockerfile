FROM gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed

LABEL maintainer="Gerald Mboowa"
LABEL description="Production AMR container for rMAP-Candida ChroQueTas/FungAMR workflow"
LABEL version="2026.07-chroquetas-v9"
LABEL base_image="gmboowa/rmap-myc-candida-amr:2026.05-chroquetas-v7-fixed"

SHELL ["/bin/bash", "-lc"]

USER root

ENV MINIPROT_ENV=/opt/miniprot_env

RUN set -euo pipefail; CONDA_SOLVER=classic conda create -y -p "${MINIPROT_ENV}" -c conda-forge -c bioconda miniprot; CONDA_SOLVER=classic conda clean -a -y; ln -sf "${MINIPROT_ENV}/bin/miniprot" /usr/local/bin/miniprot

RUN set -euo pipefail; which run_fungamr_scan; which ChroQueTas.sh; which miniprot; miniprot --version

WORKDIR /data

ENTRYPOINT []

CMD ["bash"]
