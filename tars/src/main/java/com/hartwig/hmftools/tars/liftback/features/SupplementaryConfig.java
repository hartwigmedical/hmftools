package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_MERGES;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_READ_OVERLAP;

// Config for SupplementaryResolver; defaults() implements the annotated-junction policy.
public class SupplementaryConfig
{
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MaxSuppMerges;
    public final boolean AnnotatedOnly;
    public final int MaxSuppReadOverlap;

    public SupplementaryConfig(
            final int minIntronLength, final int maxIntronLength, final int maxSuppMerges,
            final boolean annotatedOnly, final int maxSuppReadOverlap)
    {
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MaxSuppMerges = maxSuppMerges;
        AnnotatedOnly = annotatedOnly;
        MaxSuppReadOverlap = maxSuppReadOverlap;
    }

    public static SupplementaryConfig defaults()
    {
        // AnnotatedOnly=false: a merge builds a better alignment from existing bwa records, so missing annotation must not block it.
        // An annotated junction position is still preferred where one exists.
        return new SupplementaryConfig(
                MIN_IMPLIED_INTRON_LENGTH, MAX_IMPLIED_INTRON_LENGTH,
                MAX_SUPP_MERGES,
                false, MAX_SUPP_READ_OVERLAP);
    }
}
