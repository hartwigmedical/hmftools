package com.hartwig.hmftools.viridian.integration.extract;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

// An SV which might be a viral integration into the host genome.
public record HostVariantCandidate(
        StructuralVariantType type,
        String filter,
        HostBreakend startBreakend,
        // Null for an SGL.
        @Nullable HostBreakend endBreakend,
        String insertSequence,
        boolean lineSite,
        @Nullable InsertRepeat insertRepeat,
        // ESVEE's own alignments of the insert against the host genome. Empty when it placed nowhere in the host.
        String insertHostAlignments
)
{
    public HostVariantCandidate
    {
        if(insertSequence.isEmpty())
        {
            throw new IllegalArgumentException("insertSequence cannot be empty");
        }
    }
}
