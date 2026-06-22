package com.hartwig.hmftools.tars.liftback.rescue;

import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_BOUNDARY_SHIFT;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_CHAIN_DEPTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MIN_ANCHOR_OVERHANG;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MIN_PARTIAL_MATCH_RUN;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_SOFTCLIP_TOLERANCE;

// Config for JunctionRescueResolver. Defaults for the annotated-junction policy:
// min anchor overhang 3, min intron 21, max intron 1_000_000.
public class RescueConfig
{
    public final boolean Enabled;
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MinAnchorOverhang;
    public final int MaxChainDepth;
    public final boolean AnnotatedOnly;
    public final int SoftclipTolerance;
    public final int MaxBoundaryShift;
    public final int MinPartialMatchRun;

    public RescueConfig(
            final boolean enabled, final int minIntronLength, final int maxIntronLength,
            final int minAnchorOverhang, final int maxChainDepth, final boolean annotatedOnly,
            final int softclipTolerance, final int maxBoundaryShift, final int minPartialMatchRun)
    {
        Enabled = enabled;
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MinAnchorOverhang = minAnchorOverhang;
        MaxChainDepth = maxChainDepth;
        AnnotatedOnly = annotatedOnly;
        SoftclipTolerance = softclipTolerance;
        MaxBoundaryShift = maxBoundaryShift;
        MinPartialMatchRun = minPartialMatchRun;
    }

    public static RescueConfig defaults()
    {
        // tolerance=0 preserves strict-complementary semantics for older tests; production uses enabledDefaults().
        return new RescueConfig(
                false,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                true, 0, 0, DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }

    public static RescueConfig enabledDefaults()
    {
        // AnnotatedOnly=false: merge builds a better alignment from existing BWA records, so missing
        // annotation shouldn't block it. Annotated snap preferred when available; otherwise trust-primary.
        return new RescueConfig(
                true,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                false, DEFAULT_SOFTCLIP_TOLERANCE, DEFAULT_MAX_BOUNDARY_SHIFT, DEFAULT_MIN_PARTIAL_MATCH_RUN);
    }
}
