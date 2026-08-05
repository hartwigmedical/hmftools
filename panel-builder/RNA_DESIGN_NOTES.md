# PanelBuilder RNA Support — Design Notes

Working design notes for adding RNA panel support to PanelBuilder. Captures the plan,
decisions made, and open issues to revisit. Not user-facing documentation.

## Outstanding follow-ups

Single source of truth for what remains. Each is tagged in code where practical (search the tag). Detail for several
is in the sections further down. Everything in the batch/merge records below (A1–B6, M1–M5) is DONE.

| # | Item | Code tag / location | Notes |
|---|------|---------------------|-------|
| 1 | RNA generation performance (end-to-end run spends minutes on RNA) | *(untagged)* | Root-caused: BWA-MEM alignment of junction-crossing candidates. Option A done (1.45x); A+ in progress. See "Performance — RNA generation". |
| 2 | Quality score of non-contiguous (spliced / SV) probes distorted | `FIXME` `ProbeQualityModel.computeFromAlignments` | `targetScore` assumes a full-length on-target match. Pre-existing (affects SV too). See "Open issue — quality score…". |
| 3 | RefSeq/NM transcript resolution disabled | `TODO` `GenesRna.resolveTranscript` | Validate the non-1:1 Ensembl↔RefSeq mapping, then re-enable. |
| 4 | Part-coding exon classified as fully coding | `TODO` `GenesRna.createTargets` | A long exon with few coding bases is tiled entirely as coding. Reconsider. |
| 5 | `PanelData` getters return live internal lists | `TODO` `PanelData` | No live aliasing bug found; defensively copy (also `ProbeGenerationResult` ctor). Own commit. |
| 6 | Rejected spliced (multi-region) probe reporting | `TODO(RNA)` `ProbeGenerationResult` | Decide how a rejected non-variant multi-region probe is reported. |
| 7 | `SingleProbe` requires exact `PROBE_LENGTH` | *(untagged)* | Define the can't-fill fallback for a short-exon padded probe. |
| 8 | RNA strandedness | *(design decision, deferred)* | RNA is single-stranded; currently emit genome-forward only. Decide one strand vs both. |
| 9 | README "exon aware tiling algorithm" section | *(doc)* | Complete the user-facing description. |

## Implementation status

- **A1 [DONE]** — `SequenceDefinition` segment-list model + segment types + `BasicProbeLayout` + generalized `probeTargetedRegions` +
  region-based ordering.
- **A2 [DONE, folded into A1]** — multi-region probes treated as complex by `sequenceIndelSize` (→ alignment model, skip coverage check).
- **B1 [DONE]** — `RegionMapping` (probe-space ↔ genome conversions) + `SequenceDefinition.spliced` converter.
- **B2 [DONE]** — RNA probe-space tiler. **B2a** = tiling math extracted to shared `ProbeTiling.calculateOptimalProbeTiling` (pure move, DNA
  behaviour preserved); **B2a.2** = added `ProbeTiling.calculateContainedTiling(region, pinStart, pinEnd)` — contained tiling with
  edge-pinning flush to boundaries (the RNA/final-design math; DNA path untouched); **B2b** = `RnaProbeGenerator.coverExonRange` —
  probe-space candidate generation + acceptable-subregion split + shifting + size-rule tiling (short→spliced padding across the junction,
  mid→single centred, long→edge-pinned tiling), plus `RegionMapping.isRegionBoundary` for pin selection. All unit-tested (
  `RnaProbeGeneratorTest`, `RegionMappingTest`). Not yet wired (B5). Candidate target regions are added by the caller (mirrors
  `coverUncoveredRegion`).
- **B3 [DONE]** — `GenesRna` (mirrors `Genes`): loads the RNA genes file (enum `Column` with `DelimFileReader`), resolves transcripts
  (default all merged; else Ensembl `TransName` — **RefSeq/NM resolution disabled for now**, needs validation of the non-1:1 mapping),
  merges
  exons via the new shared `GeneUtils.mergeExons`, builds the per-gene exon `RegionMapping`, and computes strand-aware target ranges
  (`createTargets`). **Whole-exon classification:** an exon with any coding base is one coding target (whole exon); a fully noncoding exon
  is a
  single 5'/3' UTR target (by position vs coding start, strand-aware). Coding always included; UTR optional per flag. TODO left to
  reconsider
  the part-coding-as-fully-coding rule (long exon, few coding bases). Also contains the generation glue driving
  `RnaProbeGenerator.coverExonRange` per exon target range into a `PanelData` + gene stats. Unit tested (`GenesRnaTest`: forward/reverse
  target
  ranges over coding/noncoding/part-coding exons, UTR gating, transcript resolution). Input file columns: `GeneName`, `Include5UTR`,
  `Include3UTR`, `TransNames`. RNA QS/GC criteria + `GENE_RNA` target type added. Not yet wired to config/app/output.
- **B4 [DONE]** — RNA output writers, separate from the DNA output. New files, all `rna_` prefixed: `rna_probes.tsv`/`.fasta`/`.bed`,
  `rna_targets.bed`, `rna_panel.bed`, `rna_rejections.tsv`/`.bed`, `rna_candidate_targets.bed.gz`, `rna_gene_stats.tsv`. The RNA probes TSV
  lists all segments generically in a `Segments` column (`;`-joined; ref segments as `region:orientation`, insert segments as their bases)
  since a spliced probe can exceed two regions; probes sorted by full `SequenceDefinition`.
  **Factored out** the per-panel file set into `ProbeOutputWriter` (probes TSV/FASTA/BED, covered target regions, covered regions, candidate
  target regions, rejected features, gene stats). `OutputWriter` now composes one DNA instance + one optional RNA instance (`rnaOutput`
  flag), plus the DNA-only verbose `candidate_probes` and `sample_variant_info`. The DNA/RNA difference is a single boolean: the probes TSV
  layout (`BasicProbeLayout` start/end, ≤2 regions vs generic region list) and the FASTA id prefix. `GeneStats` was unified into one shared
  record (was duplicated in `Genes`/`GenesRna`). Probe FASTA ids are now per-panel with separate counters — `dna_{i}` / `rna_{j}` — so ids are
  unambiguous across panels (this changes the DNA `probes.fasta` label from `probe{i}` to `dna_{i}`; DNA file names are otherwise unchanged).
  Unit-tested (`OutputWriterRnaTest`).
- **B5 [DONE]** — wired end-to-end. `-rna_genes` config (`PanelBuilderConfig`); separate RNA `PanelData`; `RnaProbeGenerator.construct`
  (own `ProbeEvaluator`); `PanelBuilderApplication.generateRnaGeneProbes()` guarded by config, driving `GenesRna.generateProbes` into the RNA
  panel, then the B4 writers. RNA is a fully separate panel — no overlap interaction with DNA. Verbose RNA candidate-probe output
  (`rna_candidate_probes.tsv.gz`) uses the same segment-based format as `rna_probes.tsv` (the DNA candidate format's `BasicProbeLayout` can't
  represent spliced probes). The Ensembl cache is loaded once and shared between DNA and RNA generation.
- **B6 [DONE]** — README RNA section + these notes updated. End-to-end validated externally as byte-for-byte identical DNA output (the
  intentional `probes.fasta` id-label change, `probe{i}` -> `dna_{i}`, is the only DNA difference). Remaining follow-ups are consolidated in
  the "Outstanding follow-ups" table at the top (tagged in code where practical). Ensembl double-load, RNA verbose candidate-probe output, the
  `dna_`/`rna_` file rename, `Segments` orientation, and gene-name-only target extra are all now done.
- **M1–M5 [DONE]** — unified DNA + RNA probe generation into one `ProbeGenerator` path over a `RegionMapping` (identity for DNA, exons for
  RNA); `RnaProbeGenerator` deleted. See the "Merge plan" section. DNA output validated byte-identical; RNA output changed only at unpinned
  rejection-split sub-ranges (uncovered-budget fix + DNA-consistent extension → more coverage).

Everything above is additive/behaviour-preserving for DNA; all DNA and RNA unit tests green.

## Decisions (B4/B5)

- **Forward strand only.** RNA probes are emitted genome-forward (matching the DNA FASTA convention); single-strand output of a chosen strand
  is deferred (see the strandedness open issue).
- **Per-panel file prefixes.** DNA and RNA produce the same set of files, prefixed `dna_` / `rna_` (e.g. `dna_probes.tsv` /
  `rna_probes.tsv`). Originally the DNA files were left unprefixed to keep DNA output byte-identical during validation; once validated, the
  `dna_` prefix was applied. `sample_variant_info.tsv` is not a per-panel file and stays unprefixed.

## Goal

Add a new "RNA transcript" region type. User supplies a list of genes; per gene chooses
coding / 5'UTR / 3'UTR and optionally a subset of transcripts (NM or Ensembl). Exon-aware
probe design. Separate RNA config file and separate RNA output (fasta, bed).

## Exon-aware design rules (from spec)

- Exon shorter than probe length: centre one probe in the exon, pad out to probe length
  with sequence from the adjacent exon (spliced probe).
- Exon between probe length and probe length + X: single centred probe.
- Exon longer than probe length + X: more than one probe, tiled evenly.
- Rejected region inside an exon: treat each side as a separate exon and repeat the above.
- Probes must only cover exons (never intronic bases).

## Architectural approach

**One unified probe model, not a DNA/RNA fork.** The pipeline currently funnels through
`SequenceDefinition` = (startRegion, startOrient, insert, endRegion, endOrient) — at most two
genome regions plus one insert. RNA probes contained within a single exon are still a single
contiguous genome region and reuse the existing machinery directly. Only two things are
genuinely new: short-exon spliced probes (2–3 disjoint regions) and exon-bounded tiling.

### Implementation batches

Small, independently reviewable batches. Phase A = refactoring of shared code, each
behaviour-preserving: **validate by regenerating an existing DNA panel and diffing to byte-identical
output**, plus green DNA tests. Phase B = RNA, added as new code that does not touch DNA paths;
validate by unit tests, then a synthetic end-to-end panel. Nothing DNA-facing changes in phase B.

**Phase A — refactoring up front (DNA output must stay identical)**

- **A1. Generalize `SequenceDefinition` to a segment list. [DONE]** As built:
    - `record SequenceDefinition(List<SequenceSegment> segments)`. Segment types are their own files:
      `SequenceSegment` (sealed interface, `baseLength()`, `Comparable`), `RefSegment(region,
    orientation)`, `InsertSeqSegment(sequence)`. Each segment validates itself in its own
      constructor.
    - No public positional constructor. Scenario factories only: `singleRegion`, `forwardSgl`,
      `reverseSgl`, `variant(start, startOrient, insert, end, endOrient)`. RNA builds the canonical
      `List` constructor directly.
    - Canonical constructor invariants: ≥1 region; no empty inserts; no two consecutive inserts (the
      caller must collapse them); no consecutive regions that are genome-adjacent with the same
      orientation (should be one region).
    - `SequenceDefinition implements Comparable`: compares the segment lists as tuples (element-wise,
      shorter-if-prefix). `SequenceSegment` is `Comparable` — `RefSegment` by region then orientation,
      `InsertSeqSegment` by sequence, cross-type by an arbitrary stable type rank. No toString compare.
    - The start/insert/end view is extracted to a separate `BasicProbeLayout.from(def)` (explicit shape
      matching + validation; throws for >2 regions / multiple inserts). Used by `OutputWriter` and
      `sequenceIndelSize`. Not a method on `SequenceDefinition`.
    - `SequenceUtils.buildSequence` maps each segment to its `SequenceData` then concatenates.
    - `ProbeUtils.probeTargetedRegions` generalized to walk all segments (no ≤2 assumption).
    - `OutputWriter` sorts probes by `SequenceDefinition` (region-based) and rejects >2-region probes
      from the DNA panel output. Rename: `isSpliced()` → `isMultiRegion()`.
    - *Validation:* DNA output content unchanged; only multi-region probe **row ordering** in
      `probes.tsv`/`.fasta` changes (now region-sorted, was grouped-last) — a deliberate improvement.
      Single-region-only panels stay byte-identical. All existing DNA tests green.
- **A2. QS routing for non-contiguous sequences. [mostly folded into A1]** `sequenceIndelSize`
  returns empty for anything not expressible as a `BasicProbeLayout` (>2 regions / complex), so
  `canUseProfile` and `PanelCoverage.needsCoverageCheck` already treat multi-region probes as complex
  (→ alignment model, skip coverage check). Remaining explicit-`isMultiRegion()` short-circuiting, if
  any, rides along with the RNA work when spliced probes actually exist.

**Phase B — RNA, additive (no DNA path touched)**

- **B1. `RegionMapping` + converters. [DONE]** Concrete `RegionMapping` (position/region conversions
  both directions; genome→probe partial→`OptionalInt`, probe→genome total→throws only on bad input) +
  `SequenceDefinition.spliced` converter. Unit-tested, not wired.
- **B2. RNA probe-space tiler. [IN PROGRESS]** Separate from the DNA tiler; drives off `RegionMapping`
  in probe-space; candidate generation, rejected-region split/recurse, forced junction-crossing.
    - **B2a [DONE]** — extracted the tiling math to shared `ProbeTiling.calculateOptimalProbeTiling`
      (pure move; DNA tests gate; unit-tested directly). Both DNA `ProbeGenerator` and the future RNA
      tiler use it.
    - **B2b [DONE]** — `RnaProbeGenerator.coverExonRange(mapping, rangeStart, rangeEnd, ...)`: probe-space
      candidate generation, acceptable-subregion split, and probe shifting (the probe-space analogues of
      `ProbeGenerator.coverUncoveredRegion` / `coverAcceptableSubregion`), building candidate
      `SequenceDefinition`s via `RegionMapping.toGenomeRegions` and evaluating them. Applies the exon size
      rules per acceptable sub-range: short (< PROBE_LENGTH) → single spliced probe padded across the
      junction; mid (≤ PROBE_LENGTH + `RNA_EXON_SINGLE_PROBE_SLACK`) → single centred probe; long →
      **edge-pinned tiling** via `calculateContainedTiling`, pinning each edge iff it coincides with a
      `RegionMapping` boundary (`isRegionBoundary`). Input range must lie within a single exon (asserted);
      junction crossing happens only in the padding branch. Kept as a separate class from `ProbeGenerator`
      for now (composition over the shared `ProbeEvaluator`); DNA paths untouched.
- **B3. Transcript resolution + target ranges. [DONE]** `GenesRna`: loads RNA genes file; resolves
  transcripts (all-transcripts default, else Ensembl `TransName` subset — RefSeq/NM disabled for now);
  merges exons (shared `GeneUtils.mergeExons`); builds the exon `RegionMapping`; computes strand-aware
  target ranges via `createTargets` with **whole-exon classification** (any coding base → whole exon is
  a coding target; fully noncoding exon → single 5'/3' UTR target by position + strand). Coding always,
  UTR optional. Includes the generation glue over `RnaProbeGenerator.coverExonRange`. Unit-tested
  (`GenesRnaTest`); not wired to config/app/output.
- **B4. RNA output writers.** Separate RNA fasta / bed / TSV (new `RnaOutputWriter` or additions to
  `OutputWriter`). Takes `List<Probe>`; independently testable.
- **B5. Wire end-to-end.** `GenesRna.generateProbes` (already built, B3) ties B2+B3 into a separate
  `PanelData`; construct an `RnaProbeGenerator` from the ref genome / quality model; add `-rna_genes`
  config; call from `PanelBuilderApplication.run()` guarded by config; invoke B4. First point RNA output
  exists. *Validate: synthetic end-to-end panel; IGV check of the RNA bed.*
- **B6. Docs + finalise.** README RNA section; resolve this notes file; log the QS open issue as a
  tracked task.

**Later (not in initial delivery):** merge the RNA tiler into the DNA tiler by having DNA pass a
whole-genome identity mapping, with existing DNA tests as the gate.

Unification is cleaner than a mode fork. The DNA/RNA difference collapses to **which `RegionMapping`
is passed**, given one key rule: **pin a tiling edge iff it coincides with a `RegionMapping`
boundary.**

- DNA identity mapping (one whole-chromosome region): target edges are mid-chromosome, never at a
  mapping boundary → never pinned → centred with `REGION_UNCOVERED_MAX` edge gap (today's behaviour).
- RNA exon mapping: a target edge at an exon boundary is pinned; a rejection-split interior edge is
  not a mapping boundary → not pinned. This also gives the "pin only true exon boundaries" nuance for
  free.
- Junction crossing is not a separate concern: edge-pinning ends the outermost probe exactly at the
  boundary, so within-exon tiling never spills across a junction; crossing only occurs in the
  short-region padding branch.
  So `ProbeTiling` would gain `pinStart`/`pinEnd` booleans (derived from boundary coincidence; when set,
  `probeCount = ceil(length / PROBE_LENGTH)` and the pinned end is flush), and candidate/orchestration
  logic becomes mapping-driven and shared. Deferred only because it makes the well-tested DNA tiler
  mapping-aware up front (regression risk); staged so RNA is proven first.

### Model (as built in A1)

```java
record SequenceDefinition(List<SequenceSegment> segments)   // implements Comparable
sealed

interface SequenceSegment permits RefSegment, InsertSeqSegment   // Comparable, baseLength()
record RefSegment(ChrBaseRegion region,Orientation orientation)

record InsertSeqSegment(String sequence)
```

See the A1 batch entry above for the full details (factories, invariants, ordering, `BasicProbeLayout`,
generalized `probeTargetedRegions`, `isMultiRegion()`). RNA spliced probes are just multi-`RefSegment`
definitions; single-exon RNA probes remain single-region.

### Generation

- New `RnaTranscripts` generator mirroring `Genes.java`. Per gene → per transcript → per exon,
  applying the size rules. Writes into a **separate `PanelData` instance** so RNA and DNA
  coverage/overlap never interact.
- Tiling is driven by the region mapping (see below). Introns are absent from the RNA mapping, so
  "no intron spill" is automatic — the tiler can only extend into mapped (exon) sequence. Reuse
  `coverUncoveredRegion`'s acceptable-subregion split for the "rejected region → split and repeat"
  rule.
- Spliced probes come from **short regions**, decided by length — not only from whole short exons.
  A long exon (≥ PROBE_LENGTH) split at an internal rejected region can leave a short acceptable
  sub-range, and padding it to PROBE_LENGTH must cross the exon boundary (the rejection blocks the
  inward side). So the tiling recursion also triggers splicing.

- **Edge-pinned tiling (RNA requirement).** When placing multiple probes within one exon, the
  outermost probes must line up exactly with the exon boundaries (first probe start = exon start,
  last probe end = exon end), rather than leaving a gap between the outermost probe and the boundary.
  This gives good coverage of the sequence right at the splice junction. The existing DNA tiling
  (`ProbeTiling.calculateOptimalProbeTiling`) does NOT do this — it centres the tiling on the region
  and deliberately permits up to `REGION_UNCOVERED_MAX` uncovered bases at each edge (fine for DNA,
  where sequencing overhang captures the edges). So RNA needs a distinct edge-pinned tiling:
    - `probeCount = ceil(length / PROBE_LENGTH)`; first probe pinned to the start boundary, last to the
      end boundary, the rest evenly spaced (spacing `(length − PROBE_LENGTH) / (probeCount − 1)`), so
      overlap is spread evenly and there is no edge gap or extension beyond the boundary.
    - Pin an edge only when it is a true exon boundary. After a rejected-region split, a sub-range edge
      that abuts the rejected region is NOT an exon boundary — don't pin there (that edge behaves like a
      normal interior tiling edge; only the exon-adjacent outer edges are pinned).
    - Single-probe case (exon between PROBE_LENGTH and PROBE_LENGTH + X) can't pin both ends — it stays
      a single centred probe per the size rules.
    - Implemented as `ProbeTiling.calculateContainedTiling(region, pinStart, pinEnd)` [DONE] —
      contained within the region, flush at pinned edges, small gap allowed at unpinned edges. The RNA
      tiler selects `pinStart`/`pinEnd` per whether each edge coincides with an exon boundary. DNA keeps
      the centred `calculateOptimalProbeTiling`; the two converge at the later merge.

### Region mapping (the unifying abstraction)

A **region mapping** is an ordered, non-overlapping list of ref-genome regions treated as
contiguous in probe-sequence space. Everywhere probe-generation math used to reach outside the
input region via raw chromosome coordinates, it instead consults the mapping.

- **DNA**: mapping = one region per chromosome `[1, chromosomeLength]`. Identity — extension just
  continues along the chromosome, i.e. current behaviour.
- **RNA**: mapping = the transcript's ordered exons. Extension beyond the current exon lands in the
  adjacent exon; a probe window crossing a boundary maps to a spliced (multi-region)
  `SequenceDefinition`.

`RegionMapping` is probe-agnostic — it deals in genome regions/positions, never `SequenceDefinition`
(which is a probe concept). A separate converter at the probe layer,
`SequenceDefinition.spliced(List<ChrBaseRegion>)`, turns regions into a genome-forward
`SequenceDefinition` (strand handling, if any, is the caller's concern given the genome-forward
output decision).

As built (`RegionMapping.java`) — symmetric position/region conversions, no baked-in "single region"
assumption; the caller checks any condition it needs:

```java
class RegionMapping
{                                            // built from ordered, non-overlapping, ascending regions
    int length();
    OptionalInt toProbeSpacePosition(String chromosome, int position);  // genome pos -> probe-space pos; empty if unmapped
    BasePosition toGenomePosition(int spacePosition);                   // probe-space pos -> genome pos
    List<ChrBaseRegion> toGenomeRegions(int spaceStart, int spaceEnd);  // probe-space range -> genome regions, split at boundaries
}
```

**Tiler input is a probe-space range + the mapping, never a genome `ChrBaseRegion`.** Probe-space is
the spliced-contiguous coordinate system, so a probe-space range cannot express an inconsistent
junction crossing — mapping back to genome (`toGenomeRegions`, which splits at boundaries) is always
well-defined. To enter probe-space, the caller converts a single-exon target region via
`toProbeSpacePosition` on its start and end and checks the two map to a contiguous span
(`end − start + 1 == baseLength`) — that check is how the caller enforces "maps to a single region",
rather than the mapping throwing. So the "input range already crosses a junction" bug is caught at the
boundary, and the tiler itself only ever sees probe-space integers.

**Padding is subsumed into the tiler — no separate padding constructor.** The tiler and candidate
generation run in probe-space integers, then map each probe-space window to genome regions via the
mapping. The probe-space analogues of `allOverlappingProbes`, `calculateOptimalProbeTiling`, and the
shift/centre logic replace chromosome-clamped math. Junction-crossing happens only when the tiling
is *forced* to spread out — by rejected regions or tiling constraints — not for normal edges. The
tiling math already prefers minimal probes and centre-on-tie, so a long region won't gratuitously
cross junctions; a short region forces the crossing. This satisfies "centre in exon" while also
selecting the best-QS/GC candidate. The earlier `SplicedTranscript` / `buildPaddedProbe` idea is
subsumed and dropped.

The rejected-region recursion is unchanged: split into acceptable sub-ranges, tile each; short ones
get probes that extend across the mapping if a valid + acceptable candidate exists.

### Edge cases (probe needs adjacent-exon sequence but none available)

All reduce to "no valid + acceptable candidate → no probe":

1. Inward side blocked by a rejected region, outward side has no adjacent exon (transcript end) →
   every candidate is rejected or runs off the mapping → no probe.
2. Ran out of mapped sequence on one side, other side is open exon → the tiler's clamp-and-shift
   pushes the deficit to the side with available sequence → asymmetric pad. (Verified:
   `calculateOptimalProbeTiling` clamps `tilingStart` to `tilingBounds` and shifts on far-end
   overflow.)
3. Out of mapped region and out of acceptable region → no probe.

### Trade-off / risk

Making the **core tiler mapping-aware** is a more invasive refactor of well-tested DNA code (not
"tiler untouched"). Hard requirement: the whole-genome mapping must be an exact identity so all
existing DNA tests stay green with zero behaviour change — this is the main validation gate for the
refactor. Spliced candidates always route to the BWA alignment model (window profile invalid), so
short-exon RNA generation is slower per candidate; bounded because short regions have few candidate
positions.

- Exon classification is **whole-exon, not a coding/UTR split within an exon** (revised from the earlier
  per-exon split idea): an exon with any coding base is one coding target covering the whole exon; a fully
  noncoding exon is one 5'/3' UTR target. 5' vs 3' is by exon position relative to the gene coding start,
  **strand-dependent** (below coding start → 5' on forward, 3' on reverse). This keeps both target edges on
  true splice junctions (pinned flush) and lets 5'/3' UTR be toggled independently, at the cost that a
  part-coding exon's UTR bases are covered as coding and not attributed to a UTR feature. TODO in
  `GenesRna.createTargets` to reconsider (long exon, few coding bases). Contrast: DNA `Genes` clamps coding
  to `[CodingStart, CodingEnd]` (±`GENE_CODING_REGION_EXPAND`) and gives a fully noncoding exon just one
  centred UTR probe, with no 5'/3' distinction.
- RNA constants in `PanelBuilderConstants` (`RNA_EXON_SINGLE_PROBE_SLACK` = the size-rule X;
  `GENE_RNA_QUALITY_MIN` / `GENE_RNA_GC_TARGET` / `GENE_RNA_GC_TOLERANCE` criteria).

### Transcript resolution (Ensembl; NM deferred)

**Default = all transcripts for the gene, merged.** Unlike DNA `Genes` (which defaults to the
canonical transcript plus optional extras), RNA defaults to every transcript of the gene with exons
merged via `GeneUtils.mergeExons`. An optional input field (`TransNames`) specifies an exact subset of
transcripts instead.

- No subset given → `EnsemblDataCache.getTranscripts(geneId)` → merge all.
- Subset given → resolve each listed ID against `TranscriptData.TransName` (Ensembl).
- **NM/RefSeq resolution is disabled for now** (`GenesRna.resolveTranscript` matches `TransName` only). It
  needs validation first: the Ensembl↔RefSeq mapping is not 1:1 and `RefSeqId` may be null or multi-valued.
  Re-enable via a `RefSeqId` fallback with clear not-found / ambiguous errors once validated.

## Decisions

- **QS for spliced probes**: treat as a pre-existing issue (see below), document, revisit later.
  Do not block RNA on it.
- **FASTA orientation**: genome forward (matches current DNA convention). Do NOT reverse-complement
  reverse-strand transcripts.
- **Output separation**: separate DNA and RNA output files where appropriate. Keep the DNA probe
  output format the same — it will no longer match the generalized `SequenceDefinition` model, so
  assert in code at output time that DNA probes have ≤2 region segments.
- **Phasing**: model refactor is a standalone first step (pure refactor, DNA tests green).

## Open issue — quality score of non-contiguous constructed sequences

**This is a pre-existing issue, not RNA-specific.** `ProbeQualityModel.computeFromAlignments`
assumes `alignments.get(0)` is the full-length on-target exact self-match and uses its score as the
`targetScore` denominator that normalises the quality score.

For a sequence that is **not contiguous in the reference genome**, BWA's best hit is only a partial
(single-fragment) alignment with a lower score, so `targetScore` is understated and the resulting
quality score is distorted (off-target risk overstated).

- **Already affects** SV probes that join distant or cross-chromosome loci (the constructed
  sequence's best alignment is ~half-length).
- **Small INDELs escape** it because the sequence stays near-contiguous, so BWA still finds a
  near-full-length alignment.
- **RNA spliced probes** are the same class as SV probes — the issue just triggers more often.

Routing to the profile vs the alignment model is decided in `ProbeQualityScorer`:

- `canUseProfile` (`isSequenceSimilarToRef`) sends anything not near-contiguous with the ref to the model:
  >2 regions, cross-chromosome, or a same-chromosome gap/insert larger than `PROBE_QUALITY_PROFILE_MAX_REF_DIFF`.
  So a spliced probe across a normal (large) intron already goes to the model.
- **Short-region guard (generic, not RNA-specific):** the profile scores each constituent region against
  fixed-length windows and cannot score a region shorter than one window (`ProbeQualityProfile.baseWindowLength()`,
  40b) — that similarity measure is length-dependent and unvalidated below the window length. So `ProbeQualityScorer`
  routes any probe with a region shorter than the window to the model (which scores the full-length probe sequence).
  This is what makes RNA junction-crossing candidates (a few bases in one exon) and short-exon spliced probes safe,
  and also covers any DNA variant probe with a very short region. Previously this threw
  `probe length must be >= 40` from `ProbeQualityProfile.computeQualityScore`.

These model-routed spliced probes are exactly the class hit by the `targetScore` distortion above.

Candidate fixes to revisit:

1. Use a theoretical `targetScore = probeLength * matchScore` for non-contiguous sequences instead
   of reading it from alignment 0. Small, local.
2. Align against an exon-junction-aware transcriptome index so the on-target match is contiguous.
   More accurate, heavier (new resource + index).

Action: prototype (1), compare against hand-checked probes, decide. Track as its own task.

## Performance — RNA generation (investigation + fixes)

Investigated on branch `AUS-422-panel-builder-rna-v1` (2026-08-04), validated on the public `panelbuilder/validate` set (489 RNA genes,
GRCh38). Timings are the isolated RNA-generation phase (27s fixed load overhead subtracted), 13 threads.

- **Root cause.** RNA gene generation runs single-threaded on `main`; the only parallelism is inside BWA-MEM's own thread pool. Junction-
  crossing (spliced) and short-region candidate probes route to the BWA **alignment model** (`ProbeQualityModel`) rather than the profile,
  because the profile cannot score a non-contiguous sequence or a region shorter than one window (40b). `ProbeGenerator` generated a
  candidate at every probe-space position across each of ~10,713 merged-exon targets and evaluated all of them — ~2.4M of ~5.7M candidates
  routed to BWA, to keep only ~3,516 spliced probes (>99.8% discarded).
- **Profiling (JFR, full RNA run).** `BwaMemIndex.createAlignments` = ~74% of main-thread runnable samples; the rest is one-time profile-file
  loading (~11%, float-parse + gz-inflate). Candidate generation and profile scoring are negligible. So **BWA query volume is the only
  meaningful lever** — batching across exons (limited, threads already ~68% saturated) and profile-path micro-optimisation are not worth it.
- **Option A [DONE].** Two-pass candidate generation in `ProbeGenerator.coverMappedRange`: pass 1 = within-exon (single-region, profile-scored)
  candidates; pass 2 = spliced candidates only for target positions pass 1 leaves uncovered (short exons; short sub-ranges abutting a junction
  after a rejection split), bounded to the junction neighbourhood. DNA identity mapping never splices, so pass 2 is inert (DNA unaffected).
  Result: RNA phase 326s -> 225s (**1.45x**). Output byte-identical (`rna_probes`/`panel`/`rejections`/`gene_stats`/`candidate_targets`;
  `rna_targets.bed` differs in tie-order only). All unit tests green.
- **Option A+ [DONE].** A short exon (< `PROBE_LENGTH`) still BWA-evaluated every junction-crossing position to keep one centred padded probe.
  The single-probe placement only ever selects within `PROBE_SHIFT_MAX` of the centred tiling position, so pass 2 evaluates just that window
  for short-exon targets (no single-region candidate), falling back to the full sweep only when the centred candidate is unacceptable (the
  placement may then shift and depend on the full acceptance pattern). Result: RNA phase 225s -> **189s** (cumulative **1.73x** vs baseline
  326s; total run 353s -> 216s). CPU 2221s -> 1823s (fewer BWA queries). Output byte-identical; all unit tests green.
- **Cross-exon alignment batching [tried, reverted - no gain].** After A+, wall is BWA-bound. Hypothesis: pooling each exon's junction-crossing
  alignments into one big batch across all exons would raise thread saturation. Implemented and measured against A+ (averaged, alternating, to
  remove thermal bias):
    - No batching at all (align one query at a time): **~1680s, fully serial** (user/real ~1.0) - i.e. **batching is worth ~7x**, because BWA
      only parallelises across the queries within a single align call.
    - A+ per-exon batching: **~226s** (10 threads).
    - Cross-exon batching: **~227s** (10 threads) - statistically identical to A+.
  Conclusion: A+'s per-exon batches (~11-200 queries) already exceed the core count, so they already saturate the alignment threads; pooling
  further adds nothing at ~10 cores (and was marginally *worse* at an oversubscribed 13 threads). Cross-exon batching would only help when cores
  far exceed the per-exon candidate count (many-core servers); revisit there. **Kept:** RNA generation now goes through the `CoverExonRange`
  spec / `ProbeGenerator.generateBatch` interface for consistency with the other generators (per-exon alignment batching retained). **Reverted:**
  the cross-exon pooling machinery. `ProbeGenerator.Batch` documents which spec types are pooled (alignment-model: single-probe, exon) and which
  are not (profile-scored generic specs).
- **DNA `coverRegion` batching [investigated, not worth it].** JFR of a DNA-only gene run shows **zero** alignment-model calls - DNA gene tiling
  is 100% profile-scored (contiguous single-region), dominated by the one-time profile-file load, with negligible candidate scoring. Batching
  only accelerates the alignment model (BWA cross-query parallelism), so pooling `coverRegion` candidates gives no speedup at any thread count.
  The only DNA path that reaches the model is SV/custom-variant single probes, already pooled by `singleProbeBatch`.

### Rejected — profile aggregation as a spliced-probe quality score

Side analysis over 84,648 real spliced candidates: a profile estimate = min over the constituent region fragments' profile scores vs the
alignment-model score. **Not usable.** Pearson r = 0.75 but systematically optimistic (mean 0.78 vs model 0.59); at the accept threshold it
**falsely accepts 9.5% of model-rejected spliced probes** — i.e. it overestimates quality, which is disallowed. The error concentrates where
a fragment is substantial but shorter than the 40b profile window (false-accept rate peaks at 30–39b = 3.3%; <20b or ≥40b are safe ~0.3%):
the profile is blind to that sub-window fragment, and to the artificial junction's off-targets entirely — both fundamental, not tunable. So
the alignment model stays for spliced probes; the win is fewer BWA queries (Option A/A+), not a cheaper score. (The model QS for spliced is
itself distorted low by the `targetScore` issue, follow-up #2 — comparing against it is imperfect, but does not change the conclusion.)

Do **not** change the profile window aggregation (`ProbeQualityProfile.aggregateQualityScore`, a length-weighted soft-min) — it is
empirically validated against whole-probe alignment and real probe performance, and must not be made to overestimate. A safe cheap
accept/reject for the *single-region* (profile) path is available via bounds `LB = min raw window over all overlapping windows <= aggregate <=
UB = min over fully-contained windows` (auto-accept when `LB >= threshold`, auto-reject when `UB < threshold`, exact soft-min only in the
narrow ambiguous band) — never overestimates — but the profile path is not the bottleneck, so it is not pursued.

## Other issues to keep in mind

Background/context (the actionable items are in the "Outstanding follow-ups" table at the top; detail here).

- Determinism: recent work (AUS461, AUS434) fought nondeterminism. RNA iteration/collection must use
  sorted, fixed order (transcript order, exon order).
- Testing without sample data: build fake `TranscriptData` / `ExonData` fixtures + a small synthetic
  ref genome. Validate BED in IGV. No real sample IDs anywhere.
- **Follow-up #6 detail** — `ProbeGenerationResult.rejectProbe` (`TODO(RNA)`): a rejected multi-region non-variant
  (spliced) probe currently splits single-region into rejected regions, else reports the whole probe.
- **Follow-up #5 detail** — defensive-copying: an audit found no live aliasing bug, but several classes store/return
  collections by reference while peers (`SequenceDefinition`, `RegionMapping`, `ProbeGenerationResult` factories)
  defensively copy. Priority: `PanelData` getters and the `ProbeGenerationResult` canonical constructor. Not worth it
  for read-only sites like `ProbeOutputWriter`. Own commit.

## Done in follow-up work (post initial delivery)

- **Output file naming:** DNA and RNA share one set of base names, prefixed per panel (`dna_` / `rna_`). One set of
  file-name constants in `PanelBuilderConstants` plus `DNA_OUTPUT_PREFIX` / `RNA_OUTPUT_PREFIX`; the duplicated `RNA_*`
  constants removed. `sample_variant_info.tsv` (not a per-panel file) stays unprefixed.
- `Segments` orientation emitted as `1`/`-1`; target extra info is the gene name unless a transcript subset was
  specified; Ensembl cache loaded once and shared; RNA verbose candidate-probe output added.
- **Merge probe generation:** DNA and RNA unified into one `ProbeGenerator` path (M1–M5); `RnaProbeGenerator` deleted.

Still open: see the "Outstanding follow-ups" table at the top (performance, README tiling section, plus the
code-tagged items).

## Merge plan: unify DNA + RNA probe generation

Goal: one probe-generation path parameterised by a `RegionMapping`, deleting the DNA/RNA fork
(`RnaProbeGenerator`). RNA is DNA generalised: DNA uses an identity mapping (one whole-chromosome
region); RNA uses the transcript's ordered exons. No mode boolean — the differences are parameters.

**Unifying rules**

- Everything runs in probe-space over a `RegionMapping`. Candidate windows map back to genome via
  `toGenomeRegions` (a single region for DNA, possibly spliced for RNA).
- Pinning is derived, not passed: pin a tiling edge iff it coincides with a mapping boundary
  (`isRegionBoundary`). Identity mapping has no interior boundaries → never pins → today's centred DNA
  tiling. Exon boundaries pin → flush (good splice-junction coverage).
- Rejected regions carry no special behaviour: the acceptable-subrange split treats each side as an
  ordinary region; the only constraint is that tiling cannot overlap the rejection (guaranteed because
  `probeBounds` = the acceptable subrange). `REGION_UNCOVERED_MAX` is the total uncovered budget per
  acceptable region, independent of adjacency (the old per-edge doubling in `calculateContainedTiling`
  was a bug, now fixed).
- Short-range padding is emergent, not special-cased: an acceptable subrange spanning a junction already
  contains spliced candidate windows, so the general tiler fits a probe by extending into adjacent mapped
  space (chromosome for DNA, next exon = spliced for RNA).

**Parameters replacing the fork**

- the `RegionMapping` (identity vs exon) — the core parameter; drives extension.
- the pin-boundary predicate (`IntPredicate`) — DNA `pos -> false`, RNA `mapping::isRegionBoundary`.
- eval criteria, select strategy (already parameters).

**Interface (as built)**

Core: `ProbeGenerator.coverMappedRange(RegionMapping mapping, int rangeStart, int rangeEnd, TargetMetadata,
ProbeEvaluator.Criteria, ProbeSelector.Strategy, IntPredicate pinBoundary)`. Two convenience wrappers capture
the pin policy so callers never pass the predicate: `coverRegion` (DNA — builds an identity mapping, pins
nothing; also does the existing genome-space coverage subtraction, then calls the core per uncovered
subregion) and `coverExonRange` (RNA — pins at exon boundaries). `GenesRna` calls `coverExonRange` per exon.
The other spec types (`CoverOneSubregion`, `CoverOnePosition`, `SingleProbe`) stay single-region — RNA never
needs them. Coverage-avoidance stayed in `coverRegion` (genome-space), not a core parameter, since RNA is a
separate panel with no coverage to subtract.

**Correctness gate**

DNA tiling behaviour was not covered by unit tests originally (`ProbeGeneratorTest` only tests candidate
generation), so `ProbeGeneratorTilingTest` (8 DNA tiling characterisation cases) was added first, then the
DNA-space rewrite validated by regenerating a DNA panel and diffing to byte-identical output (confirmed:
all DNA probe/rejection/panel/fasta output byte-identical; `targets.bed` content-identical, tie-order only).
Existing DNA and RNA unit tests stayed green throughout.

**Pinning as a parameter, not mapping data.** `RegionMapping` stays pure geometry (`isRegionBoundary`).
The pin policy is an `IntPredicate pinBoundary` passed to `coverMappedRange`: DNA passes `pos -> false`
(pin nothing — chromosome ends and target edges are arbitrary), RNA passes `mapping::isRegionBoundary`
(pin every exon boundary). Pinning is also suppressed for a single-probe (short/mid) range, which cannot
be flush to both boundaries and is centred (short ranges pad across the junction) — mirroring the old RNA
short/mid branches, which are otherwise emergent from the shared path.

**Phasing (each independently gated) — ALL DONE.**

- **M1 [DONE].** Unified tiler entry `ProbeTiling.calculateProbeTiling(region, probeBounds, pinStart,
  pinEnd)` — unpinned delegates to `calculateOptimalProbeTiling`, pinned to the contained-flush placement.
  Also fixed the per-edge `REGION_UNCOVERED_MAX` bug in `calculateContainedTiling`.
- **M2 [DONE].** `RegionMapping.wholeChromosome` identity factory.
- **M3 [DONE].** `ProbeGenerator.coverMappedRange` / `coverAcceptableMappedSubrange` (probe-space port of
  the DNA `coverUncoveredRegion` / `coverAcceptableSubregion`, pin-aware via the `pinBoundary` predicate).
  DNA `coverRegion` routes through it via an identity mapping (no pinning). Gated by the new
  `ProbeGeneratorTilingTest` (8 DNA tiling characterisation cases) + existing DNA tests.
- **M4 [DONE].** `ProbeGenerator.coverExonRange` convenience wrapper (pins at `isRegionBoundary`). `GenesRna`
  and the app use `ProbeGenerator` (a second instance with the RNA candidate callback). All 15
  `RnaProbeGeneratorTest` cases pass through the unified path (one synthetic interior-target characterisation
  re-baselined — it now extends into surrounding exon like DNA; not reachable via `GenesRna`, which targets
  whole exons).
- **M5 [DONE].** `RnaProbeGenerator` deleted; dead DNA `coverUncoveredRegion` / `coverAcceptableSubregion`
  removed. (`calculateContainedTiling` kept — it is the pinned tiling implementation and is directly tested.)

**Still requires your byte-identical DNA panel run** to confirm the DNA-space rewrite (M3): the 8
characterisation tests are a partial gate, not exhaustive. Convenience wrappers: `coverRegion` (DNA,
existing) and `coverExonRange` (RNA, new) capture the pin policy so callers never pass the predicate.
