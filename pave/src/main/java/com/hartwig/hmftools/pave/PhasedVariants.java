package com.hartwig.hmftools.pave;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class PhasedVariants
{
    public final int LocalPhaseId;
    public final Set<VariantData> Variants;

    public PhasedVariants(final int localPhaseId)
    {
        LocalPhaseId = localPhaseId;
        Variants = Sets.newHashSet();
    }
}
