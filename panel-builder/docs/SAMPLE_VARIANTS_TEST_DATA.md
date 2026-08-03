# Sample Variant Test Data

Generator for **synthetic** PanelBuilder "sample variant" input (Purple + Linx output). Used to exercise the sample variant
code paths end to end, producing a few probes per probe type. All data is invented; no patient data is involved.

## What it produces

For a sample id (default `FAKE01T`), under an output directory:

```
purple/       FAKE01T.purple.somatic.vcf.gz     somatic SNV/INDEL
              FAKE01T.purple.germline.vcf.gz    germline SNV/INDEL
              FAKE01T.purple.sv.vcf.gz          somatic SVs (paired breakends)
              FAKE01T.purple.cnv.gene.tsv       gene copy number (for the DEL driver)
linx/         FAKE01T.linx.svs.tsv, .breakend.tsv, .fusion.tsv, .drivers.tsv, .clusters.tsv
linx_germline/ FAKE01T.linx.germline.disruption.tsv, .linx.germline.breakend.tsv
```

One variant per probe type: somatic SNV driver, somatic SNV/INDEL non-driver, germline SNV driver, and somatic SV
fusion / amplification / deletion / disruption drivers, plus a germline SV driver.

## Reference genome coupling

Variants sit at real **GRCh38** coordinates (chr-prefixed contigs). Reference bases are read from the genome you pass in, so the
data is always valid against that genome. Run PanelBuilder with the **same** reference genome.

## Generate against a real genome

`GenerateSampleVariantsTestData` is a command line tool (`ConfigBuilder` arguments) in the test sources. Run it on the test
classpath via the exec plugin:

```sh
mvn -pl panel-builder org.codehaus.mojo:exec-maven-plugin:3.1.0:java \
    -Dexec.classpathScope=test \
    -Dexec.mainClass=com.hartwig.hmftools.panelbuilder.samplevariants.testdata.GenerateSampleVariantsTestData \
    -Dexec.args="-ref_genome /path/to/GRCh38.fa -output_dir /path/to/out -sample FAKE01T"
```

Arguments (`-ref_genome`, `-output_dir`, `-sample`) are parsed by the tool via `ConfigBuilder`; the `-Dexec.*` values are just how
Maven passes the command line. The tool only writes the input files.

The FASTA needs its `.fai` index (the same one PanelBuilder uses). Coordinates live in
`GenerateSampleVariantsTestData.curatedGrch38Specs()` — adjust them to hit genes in your panel.

## Feed to PanelBuilder

```sh
java -jar panel-builder.jar \
  -ref_genome GRCh38.fa -bwa_index_image GRCh38.fa.img -probe_quality_profile probe_quality_profile.38.tsv.gz \
  -sample FAKE01T -purple_dir out/purple -linx_dir out/linx -linx_germline_dir out/linx_germline \
  -genes genes.tsv -output_dir panel_out
```

Pair with any other real inputs (e.g. `-genes`). Sample variant probes appear in `probes.tsv` with target types
`SAMPLE_*`; `sample_variant_info.tsv` lists per-variant filtering.

## Verifying

Format only (no genome) — `SampleVariantsTestDataTest` writes files against a mock genome and reloads them through the loaders,
asserting each probe type is produced:

```sh
mvn -pl panel-builder -Dtest=SampleVariantsTestDataTest test
```

Against a real genome — run the real PanelBuilder jar with the generated directories (see "Feed to PanelBuilder"); its probes and
`sample_variant_info.tsv` are the real-genome check.
