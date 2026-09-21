package com.hartwig.hmftools.virusdetect.integration;

import static java.util.Objects.requireNonNull;

import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.jetbrains.annotations.Nullable;

public record InsertRepeat(
        String repeatClass,
        String repeatType,
        // Portion of the insert the repeat covers.
        double coverage
)
{
    @Nullable
    public static InsertRepeat from(StructuralVariant variant)
    {
        String repeatClass = variant.insertSequenceRepeatClass();
        if(repeatClass == null)
        {
            return null;
        }

        return new InsertRepeat(
                repeatClass,
                requireNonNull(variant.insertSequenceRepeatType()),
                requireNonNull(variant.insertSequenceRepeatCoverage()));
    }
}
