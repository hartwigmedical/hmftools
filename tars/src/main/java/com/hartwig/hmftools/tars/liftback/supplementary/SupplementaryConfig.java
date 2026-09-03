package com.hartwig.hmftools.tars.liftback.supplementary;

import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_REF_VERIFY_BOUNDARY_SHIFT;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_MERGES;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_READ_OVERLAP;

// Config for SupplementaryResolver; defaults implement the annotated-junction policy.
public class SupplementaryConfig
{
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MaxSuppMerges;
    public final boolean AnnotatedOnly;
    public final int MaxSuppReadOverlap;
    public final int MaxRefVerifyBoundaryShift;

    public SupplementaryConfig(
            final int minIntronLength, final int maxIntronLength, final int maxSuppMerges,
            final boolean annotatedOnly, final int maxSuppReadOverlap, final int maxRefVerifyBoundaryShift)
    {
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MaxSuppMerges = maxSuppMerges;
        AnnotatedOnly = annotatedOnly;
        MaxSuppReadOverlap = maxSuppReadOverlap;
        MaxRefVerifyBoundaryShift = maxRefVerifyBoundaryShift;
    }

    public static SupplementaryConfig defaults()
    {
        // AnnotatedOnly=false: merge builds a better alignment from existing BWA records, so missing
        // annotation shouldn't block it. An annotated junction position is preferred when one is available.
        return new SupplementaryConfig(
                MIN_IMPLIED_INTRON_LENGTH, MAX_IMPLIED_INTRON_LENGTH,
                MAX_SUPP_MERGES,
                false, MAX_SUPP_READ_OVERLAP, MAX_REF_VERIFY_BOUNDARY_SHIFT);
    }
}
