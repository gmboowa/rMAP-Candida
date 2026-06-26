# Candida-focused Kraken2/Bracken database container

Container tag used by the workflow:

```text
gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
```

Database path inside the container:

```text
/opt/kraken2_db/candida
```

Reviewer access check:

```bash
docker pull gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db
docker run --rm gmboowa/rmap-myc-candida-kraken2-bracken:2026.05-db \
  bash -lc 'ls -lh /opt/kraken2_db/candida && ls /opt/kraken2_db/candida'
```

Expected files include `hash.k2d`, `opts.k2d`, and `taxo.k2d`.

Maintain the database build recipe, source genome/taxonomy list, and Bracken build read length here so that reviewers can rebuild or audit the database if needed.
