# Candida reference-genome container

Container tag used by the workflow:

```text
gmboowa/rmap-candida-refs:2026.05
```

Reference manifest path:

```text
/opt/rmap_candida_refs/references.tsv
```

Reviewer access check:

```bash
docker pull gmboowa/rmap-candida-refs:2026.05
docker run --rm gmboowa/rmap-candida-refs:2026.05 \
  bash -lc 'cat /opt/rmap_candida_refs/references.tsv'
```

The workflow uses the manifest to select species-matched references for species-aware phylogenomics.
