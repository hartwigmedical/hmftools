package com.hartwig.hmftools.esvee.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class SliceRegionState
{
    public int PositionMax;
    public int PositionMin;

    public final List<VariantInfo> Variants;
    public final List<VariantInfo> UncappedVariants;

    public int MinPositionIndex;

    public SliceRegionState()
    {
        Variants = Lists.newArrayList();
        UncappedVariants = Lists.newArrayList();
        MinPositionIndex = 0;
        reset();
    }

    public void reset()
    {
        Variants.clear();
        UncappedVariants.clear();
        PositionMin = 0;
        PositionMax = 0;
        MinPositionIndex = 0;
    }

    public void resetUncappedVariants()
    {
        UncappedVariants.clear();
        UncappedVariants.addAll(Variants);
        MinPositionIndex = 0;
    }

    public void addVariant(final VariantInfo variant)
    {
        if(Variants.isEmpty())
        {
            PositionMin = variant.PositionMin;
        }
        else
        {
            PositionMin = min(variant.PositionMin, PositionMin);
        }

        Variants.add(variant);
        PositionMax = max(variant.PositionMax, PositionMax);
    }

    public int variantCount() { return Variants.size(); }

    public String toString()
    {
        return format("variants(%d) uncapped(%d) range(%d - %d) minPosIndex(%d)",
                Variants.size(), UncappedVariants.size(), PositionMin, PositionMax, MinPositionIndex);
    }
}
