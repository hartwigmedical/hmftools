package com.hartwig.hmftools.viridian.reporting;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;

import org.jetbrains.annotations.Nullable;

// The reporting verdict for one oncology group, plus supporting data.
// Produced for every oncology group that is present in the reads or has integrations.
public record VirusReport(
        OncologyGroup oncologyGroup,
        // Null if not a reportable virus.
        @Nullable VirusType reportingType,
        int integrations,
        boolean hasMinIntegrations,
        @Nullable VirusProfile profile,
        // Null if not profiled or Purple data wasn't usable.
        @Nullable Double copiesPerTumorCell,
        ClonalityStatus clonalityStatus,
        boolean isReported)
{
    public VirusReport
    {
        if(copiesPerTumorCell != null && profile == null)
        {
            throw new IllegalArgumentException("Copies per tumor cell needs a profile: " + oncologyGroup);
        }
    }

    public boolean isProfiled()
    {
        return profile != null;
    }
}
