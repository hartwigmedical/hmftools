package com.hartwig.hmftools.viridian.integration;

import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSv;

import org.jetbrains.annotations.Nullable;

// Possible integration of a virus genome into the tumor host genome.
public record IntegrationCandidate(
        CandidateSv sv,
        @Nullable ViralInsertAlignment alignment
)
{
    public boolean isAligned()
    {
        return alignment != null;
    }

    public boolean isIntegration()
    {
        return alignment != null && alignment.passesFilters();
    }
}
