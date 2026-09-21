package com.hartwig.hmftools.virusdetect.integration;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

// An SV which might be a viral integration into the host genome.
// ESVEE's own annotations ride along as data; nothing in this phase acts on them.
public record CandidateIntegration(
        String svId,
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
    public CandidateIntegration
    {
        if(svId.isEmpty())
        {
            throw new IllegalArgumentException("SV id is empty");
        }
        if(insertSequence.isEmpty())
        {
            throw new IllegalArgumentException("Insert sequence is empty");
        }
    }
}
