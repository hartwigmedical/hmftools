# TARS liftback policy

How liftback decides between a bwa reference (ref) placement and a transcript-contig (tx) placement, and
how it then modifies a read on the way to the genomic output BAM. Values below are the current in-code
defaults; the parenthesised name is the constant they come from.

## 0. CONSTANTS IN PLAY

### Junctions / anchors
```text
MIN_JUNCTION_ANCHOR = 8
    : WHERE  -> SpliceCommon
    : WHY    -> an N is trusted as a real splice only when BOTH flanking M runs are >= 8 bp.
    : GOTCHA -> this is the discriminator/collapse TRUST floor; it is NOT the floor for keeping an annotated
                tx junction (that is ANNOTATED_JUNCTION_MIN_ANCHOR_BP = 1). Easy to conflate the two.

ANNOTATED_JUNCTION_MIN_ANCHOR_BP = 1
    : WHERE  -> LiftBackResolver
    : WHY    -> every tx-contig junction is annotated by construction, so even a 1 bp anchor is a
                real exon base.
    : GOTCHA -> applies only to a BARE anchor; an anchor next to a softclip uses the 3 bp floor below instead.

ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP = 3
    : WHERE  -> LiftBackResolver
    : WHY    -> a softclip next to the anchor means bwa likely over-ran the boundary, so demand more evidence.

SPLICE_FLANKING_DELETION_MAX_BP = 5
    : WHERE  -> ContigTranslator
    : WHY    -> a small D straddling an exon boundary is a tx-FASTA off-by-N artefact; fold it into
                the intron N.
    : GOTCHA -> bound is inclusive (<= 5 absorbed, 6+ preserved).
```

### Scoring (bwa-mem defaults, used to re-score lifted boundaries)
```text
MATCH / MISMATCH = +1 / -4
    : WHERE  -> BwaMemScore (shared by LiftBackResolver + BoundaryReclaim)
    : WHY    -> the one bwa-mem affine model used by both the base walk and the CIGAR reconstruction.
    : GOTCHA -> XA alts carry no AS tag, so their score is RECONSTRUCTED from the CIGAR, not read
                off the record.

GAP_OPEN / GAP_EXTEND = -6 / -1
    : WHERE  -> BwaMemScore (used by LiftBackResolver's CIGAR reconstruction)
    : WHY    -> standard affine gap cost for the reconstructed score.
    : GOTCHA -> soft-clips are NOT penalised (clipped bases score 0), so a clipped alt can still win on score.
```

### Rescue (-rescue_via_supp)
```text
DEFAULT_MIN_INTRON_LENGTH = 21
    : WHERE  -> RescueConfig
    : WHY    -> shortest plausible spliceosomal intron; below it an N gap is more likely a deletion.

DEFAULT_MAX_INTRON_LENGTH = 1,000,000
    : WHERE  -> RescueConfig
    : WHY    -> upper bound on the span a primary+supp merge may bridge.

DEFAULT_MIN_ANCHOR_OVERHANG = 3
    : WHERE  -> RescueConfig
    : WHY    -> each side of the merged junction needs >= 3 matched bases to be credible.

DEFAULT_MAX_CHAIN_DEPTH = 4
    : WHERE  -> RescueConfig
    : WHY    -> cap on supplementary pieces merged into one primary (a multi-exon split read).

DEFAULT_SOFTCLIP_TOLERANCE = 5
    : WHERE  -> RescueConfig
    : WHY    -> primary/supp boundaries often disagree by 1-5 bp; snap to an annotated junction
                within this window.

DEFAULT_MAX_BOUNDARY_SHIFT = 8
    : WHERE  -> RescueConfig
    : WHY    -> bwa over-extends past the true exon boundary by up to ~8 bp; trim that far when
                probing for a donor/acceptor.

DEFAULT_MIN_PARTIAL_MATCH_RUN = 11
    : WHERE  -> RescueConfig
    : WHY    -> exon-proximal matched run required for a partial ref-verify rescue (clip with an
                unmatched outer tail).
    : GOTCHA -> longer than a full-match clip needs ON PURPOSE: the unexplained residual raises
                false-positive risk.

ref-verify match run = bwa-mem score walk (match +1, mismatch -4)
    : WHERE  -> BoundaryReclaim, via JunctionRescueResolver
    : WHY    -> the clip-vs-candidate-exon match length is the score-maximising proximal prefix, the same
                model collapse + tail-extend use; a mismatch is reclaimed only if later matches
                recover the score.
    : GOTCHA -> outer residual past the score-max prefix stays soft-clipped (e.g. 20S -> 18M 2S, not 20M);
                a longer floor (DEFAULT_MIN_PARTIAL_MATCH_RUN = 11) still gates partial-match rescues.
```

### Tail extension (-extend_softclip_tails)
```text
DEFAULT_MIN_SOFTCLIP_LENGTH = 3
    : WHERE  -> TailExtensionConfig
    : WHY    -> shorter terminal clips are not worth re-walking into the genome.

DEFAULT_MIN_EXTENSION = 3
    : WHERE  -> TailExtensionConfig
    : WHY    -> reclaim only when at least 3 bases match contiguously.

DEFAULT_MAX_EXTENSION = 30
    : WHERE  -> TailExtensionConfig
    : WHY    -> cap on bases reclaimed into the genome.
    : GOTCHA -> also bounds the annotated-junction guard probe window; caps even when more bases would match.
```

### Terminal micro-junction collapse
```text
anchor threshold = MIN_JUNCTION_ANCHOR (8)
    : WHERE  -> TerminalMicroJunctionCollapser (via SpliceCommon)
    : WHY    -> a sub-threshold terminal anchor across an intron is a candidate fabricated junction.
    : GOTCHA -> shares the SpliceCommon constant with the discriminator; changing it moves both behaviours.

reclaim scoring = bwa-mem score (match +1, mismatch -4)
    : WHERE  -> BoundaryReclaim, via TerminalMicroJunctionCollapser
    : WHY    -> the reclaimed prefix is the one maximising cumulative bwa-mem score walking from the near-exon
                boundary outward; bwa's split is kept only when the anchor scores STRICTLY higher on the far
                exon than the whole terminal window does extended contiguously. No fixed
                mismatch-count threshold.
    : GOTCHA -> NOT a "max(1, length/10)" tolerance; an internal mismatch is reclaimed when later matches
                lift the cumulative score back to its running max.
```

### Canonicalize
```text
DEFAULT_MAX_SHIFT = 5
    : WHERE  -> JunctionCanonicalizer
    : WHY    -> an intron may slide at most 5 bp onto a higher splice motif.
    : GOTCHA -> CIGAR-only: the intron position moves but the alignment START never does.

motif tiers: NONE 0  <  SEMI_CANONICAL 1  <  CANONICAL 2  <  ANNOTATED 3
    : WHERE  -> SpliceMotif
    : WHY    -> rank used to pick the best donor/acceptor; a slide must STRICTLY improve the tier.
    : GOTCHA -> classify() is strand-agnostic: it matches GT-AG and its reverse-complement CT-AC alike.
```

### MAPQ / NH
```text
RESCUE_MAPQ = 60
    : WHERE  -> LiftBackResolver
    : WHY    -> bwa's unique-placement MAPQ, and the value liftback assigns when it resolves a placement
                (discriminator swap, or rescued MAPQ-0 unique). One constant: the recognise and
                assign values coincide.

SUPP_RESCUE_MAPQ_CAP = 55
    : WHERE  -> LiftBackGroupProcessor
    : WHY    -> ceiling on a primary+supp merge, deliberately < 60 to flag a constructed alignment.
    : GOTCHA -> capping only ever lowers MAPQ; it never raises it.

NH = max(numLoci, 1)
    : WHERE  -> LiftBackGroupProcessor / LiftBackRecordOps
    : WHY    -> number of distinct genomic loci the read lifts back to.
    : GOTCHA -> counts LOCI, not emitted records; one junction across many tx contigs would
                otherwise inflate NH.

over-cap unmap = (inputMapq == 0 AND numXaAlts == 0 AND comp == REF_ONLY)
    : WHERE  -> LiftBackRecordOps.exceedsMappingCap
    : WHY    -> bwa was run with -h 75 (XA cap); a read mapping to more loci than the cap is
                emitted MAPQ 0 with NO XA (the alt list is suppressed, not truncated). Too many
                genomic places to trust, so the primary is unmapped.
    : GOTCHA -> the REF_ONLY gate matters: a tx-contig primary hits 75+ transcript contigs of one
                gene (a shared exon) that all lift to ONE genomic locus, so its suppressed XA must
                NOT be read as too many genomic places. Only a genomic (REF_ONLY) primary is
                unmapped; TX_ONLY / REF_AND_TX lift normally.
    : GOTCHA -> keyed on the INPUT state, not the post-lift MAPQ: with no XA the resolver sees a
                single locus and decidePrimaryMapq would otherwise rescue the MAPQ to 60. A unique
                read (MAPQ 60, no XA) and a few-way multimapper (MAPQ 0 WITH XA) are both excluded.
```

### Excluded regions (-rna_unmap_regions)
```text
excluded region overlap  (checked POST-lift on genomic coords)
    : WHERE  -> ExcludedRegions.excludes, via LiftBackGroupProcessor
    : WHY    -> curated rRNA / 7SL / acrocentric / multi-map contamination zones. A read lifting into one is
                contamination: the primary is UNMAPPED (REDUX-style, kept not dropped), a
                supplementary is DROPPED.
    : GOTCHA -> POST-lift, not pre-lift: a tx-contig read's input coords are chrN_tx and can't be tested against
                the genomic region list, and some excluded zones (acrocentric p-arms) ARE in the
                transcriptome, so
                contamination reaches them via tx contigs and is only visible once lifted to genomic coords.
                Primary unmap is done by flipping its LiftBackResult to UNMAPPED, so the mate is
                coordinated via the
                cache (willBeUnmapped) and the dropped supps' SA entries are removed from the primary's SA tag.
```

### Primary AS floor
```text
PRIMARY_AS_UNMAP_THRESHOLD = 30
    : WHERE  -> LiftBackGroupProcessor
    : WHY    -> bwa runs at -T 19 to surface short-anchor supps for rescue; a primary still below
                the default -T 30 AS floor after liftback (not rescued/extended/collapsed) is
                residual short-anchor noise, so it is unmapped.
    : GOTCHA -> gated on rescue running; AS is never recomputed, so post-processed primaries keep a stale-low AS
                and are exempt. Unmapped (not dropped) to keep the pair + SA references intact.

SUPP_AS_DROP_THRESHOLD = 30
    : WHERE  -> LiftBackGroupProcessor
    : WHY    -> a supplementary in the [19, 30) AS band that survived rescue + tail-extend is
                residual noise; drop it.
```

## 1. REF vs TX DECISION TREE  (LiftBackDiscriminator.categorize, then .apply)
```text
Inputs per lifted alignment set:
  hasRef / hasTx ....... at least one ref / tx alignment present
  loci ................. count of distinct genomic chrom:pos loci
  RefFullMatch ......... a ref alt matches cleanly with no softclip
  RefSoftClipped ....... a ref alt is softclipped
  TxHasNCigar .......... a tx alt carries a real N junction (both flanks >= MIN_JUNCTION_ANCHOR = 8)
  TxSoftClipAtBoundary . a tx alt softclips exactly at an exon boundary (no N)

  loci >= 2 ?
  |
  +-- YES (multi-locus)
  |     hasRef AND hasTx ?
  |     |  YES:  TxHasNCigar AND NOT RefHasNCigar -> BOTH_MULTI_TX_JUNCTION   [tx faithful: ref alts are
  |     |        |                                                             intronless paralogs/pseudogenes]
  |     |        else                             -> BOTH_MULTI               [genuine multimapper, no swap]
  |     |  tx only                                -> TX_MULTI
  |     |  ref only                               -> REF_MULTI
  |     (REF_MULTI implies zero tx alts; TX_MULTI implies zero ref alts.)
  |
  +-- NO (single locus)
        only ref                                  -> REF_SINGLE
        only tx                                   -> TX_SINGLE
        both ref and tx:
           same CIGAR and no N                    -> BOTH_AGREE                       [no contest]
           TxHasNCigar AND RefSoftClipped         -> BOTH_TX_JUNCTION_REF_SOFTCLIP    *** TX WINS ***
           TxHasNCigar AND RefFullMatch           -> BOTH_TX_JUNCTION_REF_MATCH       *** REF WINS ***
             (no softclip; clean ref through the "intron" = unspliced / retained intron / DNA contamination)
           TxSoftClipAtBoundary AND RefFullMatch  -> BOTH_TX_SOFTCLIP_REF_MATCH       *** REF WINS ***
             (tx only clipped, never spanned a junction; ref explains the read = intron retention)
           otherwise                              -> BOTH_AMBIGUOUS                   [no rule fires, no swap]

  apply() then executes the category:
    TX WINS   (BOTH_TX_JUNCTION_REF_SOFTCLIP, BOTH_MULTI_TX_JUNCTION):
                if the current primary is the ref placement, swap primary to the best tx alt (best =
                N-junction alt with fewest mismatches, else any tx alt); the old ref self stays as a
                non-primary alt, any other ref alt is dropped. note = "swapped_ref_to_tx". If self was
                already tx, drop the ref alts, note = "".
    REF WINS  (BOTH_TX_JUNCTION_REF_MATCH, BOTH_TX_SOFTCLIP_REF_MATCH):
                if the current primary is tx, swap primary to the first ref alt; note = "swapped_tx_to_ref".
                If self was already ref, drop the tx alts, note = "". (Both categories carry a ref full-match
                by construction, so there is always a ref to promote to.)
    no-swap   : REF_SINGLE, TX_SINGLE, BOTH_AGREE, BOTH_AMBIGUOUS, REF_MULTI, TX_MULTI, BOTH_MULTI leave the
                primary as-is and drop nothing.
```

## 2. HOW LIFTBACK MODIFIES A READ
```text
Passes run per name-group in this order; each feeds the next. Ref-dependent passes (rescue ref-verify,
collapse, tail-extend, canonicalize) only fire when a reference genome is loaded.

  STEP 0  TRANSLATE tx-contig -> genome           (ContigTranslator, always)
          - walk the tx CIGAR through the transcript's exons; insert an N at each exon boundary crossed.
          - a D <= SPLICE_FLANKING_DELETION_MAX_BP (5) straddling a boundary is absorbed into the N;
            a larger D is preserved after the N.
          - terminal micro-anchors below the anchor floors are folded into a softclip
            (bare anchor floor 1 bp; softclip-adjacent floor 3 bp).
          - read overhang past the first/last exon span is clamped to a leading/trailing softclip.
          - position outside any transcript span (inter-transcript spacer / past contig end) -> LIFT_FAILED.

  STEP 1  REF vs TX DISCRIMINATOR                 (see Part 1) - swaps/drops the primary at a shared locus.

  STEP 2  RESCUE via supplementary                (JunctionRescueResolver, -rescue_via_supp)
          - merge a primary's terminal softclip with a short-anchor supplementary across an annotated (or
            motif-supported) junction into one "M N M" primary; drop the merged supplementary.
          - gated on softclip complementarity, intron length [21 .. 1,000,000], anchor overhang >= 3,
            coverage overlap, and chain depth <= 4; snap point chosen by motif tier.
          - merged primary MAPQ capped at SUPP_RESCUE_MAPQ_CAP (55).

  STEP 3  TERMINAL MICRO-JUNCTION COLLAPSE        (TerminalMicroJunctionCollapser)
          - a sub-threshold terminal anchor (< MIN_JUNCTION_ANCHOR = 8) across an intron is dropped when
            the contiguous genome explains the tail at least as well as bwa's split (head-to-head bwa-mem
            score: anchor-on-far-exon vs whole-window-extended-contiguously). The score-maximising prefix is
            reclaimed into the near exon as M, the rest soft-clipped. A real short junction (anchor scores
            strictly higher on the far exon than the contiguous walk) is kept.

  STEP 4  TAIL EXTENSION                          (SoftclipTailExtender, -extend_softclip_tails)
          - walk a terminal softclip into contiguous genome (intron retention, no N); softclip >= 3,
            extension in [3 .. 30]. Skipped when an annotated intron starts inside the window and the walk
            stalls short (that range belongs to rescue).

  STEP 5  CANONICALIZE                            (JunctionCanonicalizer)
          - slide an intron up to DEFAULT_MAX_SHIFT (5) bp onto a higher splice-motif tier
            (NONE < SEMI < CANONICAL < ANNOTATED), smallest shift that strictly improves the tier with every
            moved base still matching. CIGAR only; alignment start never moves.

  RECORD MUTATIONS                                (LiftBackRecordOps / MateFieldPatcher)
          On a successful lift:
          - reference name, alignment start, CIGAR, strand and MAPQ set to the resolved genomic placement.
          - XA rebuilt to genomic alt coordinates; SA rewritten to genomic; MC set to the mate's lifted CIGAR.
          - XS:A:+/- written only when the lifted CIGAR has an N AND the transcript strand is known
            (ref-only N-cigars from rescue/tail-extend ship without XS).
          - MD always dropped (never rebuilt); NM recomputed against the genomic reference.
          - NH = max(numLoci, 1).
          - mate fields patched: mapped mate -> mate coords + TLEN (signed 5'-to-5' distance); unmapped mate
            -> mate parked on this read's locus; an unmapped read with a mapped mate adopts the mate's coords.
          On LIFT_FAILED:
          - read marked unmapped, ref/start/CIGAR/MAPQ cleared, XA and SA stripped (a supplementary whose
            own lift failed is instead mirrored onto its primary's lifted coords, keeping the 0x800 flag).

  MAPQ POLICY                                     (LiftBackResolver.decidePrimaryMapq)
          - no tx alignment in the set            -> keep bwa's input MAPQ (MAPQ-0 here is not a tx artefact).
          - discriminator swapped the primary     -> RESCUE_MAPQ (60).
          - single locus, input MAPQ 0, no unresolved hidden tie -> RESCUE_MAPQ (60).
            (hidden tie = XS == AS on a ref-only primary outside any annotated exon; it blocks the rescue.)
          - input MAPQ == 60                       -> 60 (reaffirmed).
          - otherwise                              -> keep input MAPQ.
          MAPQ is never raised by collapse / tail-extend / canonicalize; a primary+supp merge caps it at 55.

  POST-FILTERS (unmap, never drop - the pair + SA references are preserved)
          - over the XA cap (genomic primary, input MAPQ 0, no XA) -> primary unmapped: maps to too
            many genomic loci.
          - lifts into an excluded region          -> primary UNMAPPED (REDUX-style: kept, flagged unmapped, no
            cigar, placed on the mate's coords if the mate stays mapped, else *:0); supplementary
            DROPPED and its entry removed from the primary's SA. Checked post-lift on genomic
            coords (see Excluded regions above).
          - residual short-anchor primary          -> unmapped when rescue ran and post-lift AS is still
            below PRIMARY_AS_UNMAP_THRESHOLD (30) and the read was not rescued/extended/collapsed.
          - residual short-anchor supplementary     -> dropped when its AS is in the [19, 30) band after rescue.
          - any dropped supplementary (excluded / orphan / low-AS) has its SA entry stripped from the primary's
            SA tag (SaTagRewriter excludeKeys), so the primary never references a supp that is not emitted.
```

## 3. RULE PRECEDENCE  (what overrides what, highest priority first)
```text
When several rules could touch one record, they resolve in this order:

  1. UNMAPPED / LIFT_FAILED
       Input unmapped, or the position lifts outside any transcript span. Always wins: the record is unmapped
       and its tags stripped; no later rule runs on it. (A supplementary whose own lift failed is mirrored onto
       its primary's coords instead, to keep the 0x800 flag valid.)

  2. OVER-CAP UNMAP  (genomic primary, input MAPQ 0 + no XA)
       Overrides the MAPQ rescue below: even though the missing XA makes the read look single-locus (which would
       rescue MAPQ to 60), the primary is unmapped. Beats the discriminator result too - it does
       not matter which placement won if the read maps to 75+ genomic loci. Only REF_ONLY primaries
       qualify; a tx-contig primary's 75+ hits collapse to one genomic locus, so it lifts normally.

  3. REF vs TX DISCRIMINATOR  (Part 1)
       Decides which placement is primary BEFORE the refinement passes run. A swap fixes the MAPQ at 60.
       Within the discriminator, a clean ref full-match lets REF swap out a tx primary in both REF-WINS cases:
       tx spans a real N "intron" the ref reads straight through (BOTH_TX_JUNCTION_REF_MATCH = unspliced /
       retained intron / DNA), or tx only softclips at the boundary (BOTH_TX_SOFTCLIP_REF_MATCH =
       intron retention).

  4. REFINEMENT PASSES, in order: rescue -> collapse -> tail-extend -> canonicalize  (Part 2, steps 2-5)
       Each consumes the previous pass's CIGAR. Rescue caps MAPQ at 55 (overriding a 60). Collapse, tail-extend
       and canonicalize NEVER raise MAPQ and never move past what the prior pass produced.

  5. AS FLOOR UNMAP  (primary) / AS DROP (supplementary)
       Last: applied after all refinement, only when rescue ran, to records still below the -T 30 AS floor.

The full end-to-end behaviour is pinned by LiftBackEndToEndTest (full per-group engine). Its scenarios:

  exonSpanningReadLiftsToJunctionCigar          - tx read spanning two exons -> M N M, mate fields patched.
  unliftableReadIsMarkedUnmapped                - position past the contig end -> unmapped, mate flagged.
  supplementarySaTagRewrittenToGenomicCoords    - split read; supp + SA lifted to genomic coords.
  terminalSoftclipExtendedIntoContiguousGenome  - intron-retention tail walked into the genome (tail-extend).
  splitReadRescuedAcrossAnnotatedJunction       - primary softclip + short-anchor supp merged into one M N M.
  offMotifJunctionSlidToCanonical               - junction slid 2 bp onto a GT-AG motif (canonicalize).
  nativeGenomicReadPassesThroughUnchanged       - a read already on the genome is left untouched.

Component layers below the engine:

  LiftBackResolverTest                          - categorize + the resolved LiftBackResult.
  SpliceLiftBackApplyTest                       - resolve + apply record mutation: CIGAR / XA / XS / NH / unmap.
  LiftBackDiscriminatorTest                     - the full category matrix + swap/drop.
  rescue / tailextend suites                    - the individual refinement passes.
```
