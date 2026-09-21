package com.hartwig.hmftools.virusdetect.integration;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

// An SV whose inserted sequence is long enough to be worth aligning to the viral reference.
// ESVEE's own annotations ride along as data; nothing in this phase acts on them.
public record IntegrationCandidate(
        String svId,
        StructuralVariantType type,
        String filter,
        HostBreakend startBreakend,
        // Null for a single breakend, which has no other host side.
        @Nullable HostBreakend endBreakend,
        String insertSequence,
        boolean lineSite,
        @Nullable String insertRepeatClass,
        @Nullable String insertRepeatType,
        @Nullable Byte insertRepeatOrientation,
        @Nullable Double insertRepeatCoverage,
        // ESVEE's own alignments of the insert against the host genome.
        @Nullable String insertHostAlignments
)
{
    public IntegrationCandidate
    {
        if(svId.isEmpty())
        {
            throw new IllegalArgumentException("SV id is empty");
        }
        if(insertSequence.isEmpty())
        {
            throw new IllegalArgumentException("insert sequence is empty");
        }
    }

    public int insertLength()
    {
        return insertSequence.length();
    }

    public boolean isSingleBreakend()
    {
        return endBreakend == null;
    }
}
