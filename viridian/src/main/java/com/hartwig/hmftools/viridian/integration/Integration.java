package com.hartwig.hmftools.viridian.integration;

import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.HostVariantCandidate;

import org.jetbrains.annotations.Nullable;

// Possible integration of a virus genome into the tumor host genome.
public record Integration(
        HostVariantCandidate hostVariant,
        @Nullable ViralInsertAlignment alignment
)
{
    public boolean isAligned()
    {
        return alignment != null;
    }

    public boolean isPlausible()
    {
        return alignment != null && alignment.passesFilters();
    }
}
