package com.hartwig.hmftools.purple.copynumber.sv;

import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

public class StructuralVariantLegs
{
    private StructuralVariantLeg mStart;
    private StructuralVariantLeg mEnd;

    public StructuralVariantLegs(final StructuralVariantLeg start, final StructuralVariantLeg end)
    {
        mStart = start;
        mEnd = end;
    }

    public boolean hasEither() { return mStart != null || mEnd != null; }
    public boolean hasBoth() { return mStart != null && mEnd != null; }

    public StructuralVariantLeg start() { return mStart; }
    public void setStart(final StructuralVariantLeg leg) { mStart = leg; }

    public StructuralVariantLeg end() { return mEnd; }
    public void setEnd(final StructuralVariantLeg leg) { mEnd = leg; }
}
