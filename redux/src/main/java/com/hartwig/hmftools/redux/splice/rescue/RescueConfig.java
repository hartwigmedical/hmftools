package com.hartwig.hmftools.redux.splice.rescue;

// Config for JunctionRescueResolver. Defaults mirror STAR's annotated-junction policy:
// alignSJoverhangMin=3, alignIntronMin=21, alignIntronMax=1_000_000.
public class RescueConfig
{
    public static final int DEFAULT_MIN_INTRON_LENGTH = 21;
    public static final int DEFAULT_MAX_INTRON_LENGTH = 1_000_000;
    public static final int DEFAULT_MIN_ANCHOR_OVERHANG = 3;
    public static final int DEFAULT_MAX_CHAIN_DEPTH = 4;
    // Tolerance (in read bases) for overlap between primary's matched region and supp's matched
    // region. BWA's primary and supplementary often disagree by 1-5 bases on where the junction
    // sits because they use different seed/extend choices. Allowing a small overlap and snapping
    // to an annotated junction within the disputed window catches the real splice without merging
    // garbage. STAR-style annotated-junction snap is preferred; ties resolve via "trust primary".
    public static final int DEFAULT_SOFTCLIP_TOLERANCE = 5;
    // For the no-supp ref-verify path: how many over-extended matched bases to trim back into the
    // softclip when probing for an annotated donor/acceptor. BWA frequently extends past the true exon
    // boundary into the intron (a coincidental ref match), so an exact boundary probe misses the
    // junction. Measured on exp8: the over-extension caps at ~8 bases (a wider snap recovers no more),
    // so 8 covers the real range while staying bounded.
    public static final int DEFAULT_MAX_BOUNDARY_SHIFT = 8;
    // Minimum matched exon-proximal run for a PARTIAL ref-verify rescue (where bwa's softclip carries
    // an outer adapter/low-quality tail that is left clipped). A full-match clip keeps the lenient
    // MinAnchorOverhang floor; a partial match leaves an unexplained residual, so it needs a longer
    // anchor before we trust the spliced alignment over bwa's contiguous/soft-clipped one. ~11 bases
    // makes a coincidental match implausible; sweeps to 15 if the measured false-positive cost warrants.
    public static final int DEFAULT_MIN_PARTIAL_MATCH_RUN = 11;
    // Supp-merge gate: don't build a spliced read from a primary whose own MAPQ is below this. bwa
    // MAPQ 0 means >=2 equally-good placements (no trustworthy anchor); MAPQ >=1 has a clear best, so
    // the floor is 1 — reject only the fully-ambiguous primaries.
    public static final int DEFAULT_MIN_PRIMARY_MAPQ = 1;

    public final boolean Enabled;
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MinAnchorOverhang;
    public final int MaxChainDepth;
    public final boolean AnnotatedOnly;
    public final int SoftclipTolerance;
    public final int MaxBoundaryShift;
    public final int MinPrimaryMapq;
    public final int MinPartialMatchRun;

    public RescueConfig(
            final boolean enabled, final int minIntronLength, final int maxIntronLength,
            final int minAnchorOverhang, final int maxChainDepth, final boolean annotatedOnly,
            final int softclipTolerance, final int maxBoundaryShift, final int minPrimaryMapq,
            final int minPartialMatchRun)
    {
        Enabled = enabled;
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MinAnchorOverhang = minAnchorOverhang;
        MaxChainDepth = maxChainDepth;
        AnnotatedOnly = annotatedOnly;
        SoftclipTolerance = softclipTolerance;
        MaxBoundaryShift = maxBoundaryShift;
        MinPrimaryMapq = minPrimaryMapq;
        MinPartialMatchRun = minPartialMatchRun;
    }

    public static RescueConfig defaults()
    {
        // tolerance=0 preserves strict-complementary semantics for tests that pre-date the snap
        // logic; production paths use enabledDefaults().
        return new RescueConfig(
                false,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                true, 0, 0, DEFAULT_MIN_PRIMARY_MAPQ, DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }

    public static RescueConfig enabledDefaults()
    {
        // AnnotatedOnly=false in production: the merge constructs a better alignment from existing
        // BWA records (primary + supp), so annotation absence shouldn't block it. When annotation
        // does match a candidate snap point we still prefer that one (lets us recover STAR's exact
        // junction when annotated); otherwise trust-primary as fallback.
        return new RescueConfig(
                true,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                false, DEFAULT_SOFTCLIP_TOLERANCE, DEFAULT_MAX_BOUNDARY_SHIFT, DEFAULT_MIN_PRIMARY_MAPQ,
                DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }
}
