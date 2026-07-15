package com.hartwig.hmftools.tars.liftback.supplementary;

import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_BOUNDARY_SHIFT;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_SUPP_MERGES;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_SOFTCLIP_TOLERANCE;

// Config for SupplementaryResolver. Defaults for the annotated-junction policy:
// min intron 21, max intron 1_000_000.
public class SupplementaryConfig
{
    public final boolean Enabled;
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MaxSuppMerges;
    public final boolean AnnotatedOnly;
    public final int SoftclipTolerance;
    public final int MaxBoundaryShift;

    public SupplementaryConfig(
            final boolean enabled, final int minIntronLength, final int maxIntronLength,
            final int maxSuppMerges, final boolean annotatedOnly,
            final int softclipTolerance, final int maxBoundaryShift)
    {
        Enabled = enabled;
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MaxSuppMerges = maxSuppMerges;
        AnnotatedOnly = annotatedOnly;
        SoftclipTolerance = softclipTolerance;
        MaxBoundaryShift = maxBoundaryShift;
    }

    public static SupplementaryConfig defaults()
    {
        // tolerance=0 preserves strict-complementary semantics for older tests; production uses enabledDefaults().
        return new SupplementaryConfig(
                false,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MAX_SUPP_MERGES,
                true, 0, 0);
    }

    public static SupplementaryConfig enabledDefaults()
    {
        // AnnotatedOnly=false: merge builds a better alignment from existing BWA records, so missing
        // annotation shouldn't block it. Annotated snap preferred when available; otherwise trust-primary.
        return new SupplementaryConfig(
                true,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MAX_SUPP_MERGES,
                false, DEFAULT_SOFTCLIP_TOLERANCE, DEFAULT_MAX_BOUNDARY_SHIFT);
    }
}
