package com.hartwig.hmftools.viridian.reporting;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.jetbrains.annotations.Nullable;

// The reporting verdict for one virus group, plus supporting data.
// Produced for every virus group that is present in the reads or has integrations.
public record VirusReport(
        VirusGroup virusGroup,
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
            throw new IllegalArgumentException("Copies per tumor cell needs a profile: " + virusGroup);
        }
    }

    public boolean isProfiled()
    {
        return profile != null;
    }
}
