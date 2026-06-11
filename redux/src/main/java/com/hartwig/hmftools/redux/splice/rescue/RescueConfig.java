package com.hartwig.hmftools.redux.splice.rescue;

// Config for JunctionRescueResolver. Defaults mirror STAR's annotated-junction policy:
// alignSJoverhangMin=3, alignIntronMin=21, alignIntronMax=1_000_000.
public class RescueConfig
{
    public static final int DEFAULT_MIN_INTRON_LENGTH = 21;
    public static final int DEFAULT_MAX_INTRON_LENGTH = 1_000_000;
    public static final int DEFAULT_MIN_ANCHOR_OVERHANG = 3;
    public static final int DEFAULT_MAX_CHAIN_DEPTH = 4;
    // Overlap tolerance between primary and supp matched regions. BWA primary/supp often disagree
    // by 1-5 bases on junction placement; snapping to an annotated junction within this window recovers
    // the real splice. STAR-style snap preferred; ties resolve via trust-primary.
    public static final int DEFAULT_SOFTCLIP_TOLERANCE = 5;
    // How many over-extended matched bases to trim when probing for an annotated donor/acceptor
    // (no-supp ref-verify path). BWA frequently extends past the true exon boundary into the intron;
    // measured over-extension caps at ~8 bases.
    public static final int DEFAULT_MAX_BOUNDARY_SHIFT = 8;
    // Minimum exon-proximal matched run for partial ref-verify rescue (softclip with an outer
    // unmatched tail). Partial matches need a longer anchor than full-match clips because the
    // unexplained residual raises false-positive risk; 11 makes coincidental matches implausible.
    public static final int DEFAULT_MIN_PARTIAL_MATCH_RUN = 11;
    // Floor is 0 because the merge validates its own anchors; MAPQ-0 primaries are exactly the
    // tx-contig duplicate artifact the lift is designed to rescue. Knob retained for callers.
    public static final int DEFAULT_MIN_PRIMARY_MAPQ = 0;

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
        // tolerance=0 preserves strict-complementary semantics for older tests; production uses enabledDefaults().
        return new RescueConfig(
                false,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                true, 0, 0, DEFAULT_MIN_PRIMARY_MAPQ, DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }

    public static RescueConfig enabledDefaults()
    {
        // AnnotatedOnly=false: merge builds a better alignment from existing BWA records, so missing
        // annotation shouldn't block it. Annotated snap preferred when available; otherwise trust-primary.
        return new RescueConfig(
                true,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                false, DEFAULT_SOFTCLIP_TOLERANCE, DEFAULT_MAX_BOUNDARY_SHIFT, DEFAULT_MIN_PRIMARY_MAPQ,
                DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }
}
