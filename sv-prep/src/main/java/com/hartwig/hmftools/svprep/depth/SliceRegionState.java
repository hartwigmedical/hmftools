package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

public class SliceRegionState
{
    public int VariantIndexStart;
    public int VariantIndexEnd;

    public int PositionMax;
    public int PositionMin;

    public int VariantCount;

    public SliceRegionState()
    {
        reset();
    }

    public void reset()
    {
        VariantIndexStart = 0;
        VariantIndexEnd = 0;
        PositionMin = 0;
        PositionMax = 0;
        VariantCount = 0;
    }

    public void addVariant(final int variantIndex, final VariantInfo variant)
    {
        if(VariantCount == 0)
        {
            VariantIndexStart = variantIndex;
            PositionMin = variant.PositionMin;
        }
        else
        {
            PositionMin = min(variant.PositionMin, PositionMin);
        }

        VariantIndexEnd = variantIndex;
        PositionMax = max(variant.PositionMax, PositionMax);
        ++VariantCount;
    }

    public String toString()
    {
        return format("variants(%d) indices(%d:%d), range(%d - %d)",
                VariantCount, VariantIndexStart, VariantIndexEnd, PositionMin, PositionMax);
    }
}
