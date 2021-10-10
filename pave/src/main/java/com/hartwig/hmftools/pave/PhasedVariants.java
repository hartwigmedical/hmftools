package com.hartwig.hmftools.pave;

import java.util.List;

import com.google.common.collect.Lists;

public class PhasedVariants
{
    public final int LocalPhaseId;

    private final List<VariantData> mVariants;

    public PhasedVariants(final int localPhaseId)
    {
        LocalPhaseId = localPhaseId;
        mVariants = Lists.newArrayList();
    }

    public List<VariantData> variants() { return mVariants; }

    public void addVariant(final VariantData variant)
    {
        if(!mVariants.contains(variant))
            mVariants.add(variant);
    }
}
