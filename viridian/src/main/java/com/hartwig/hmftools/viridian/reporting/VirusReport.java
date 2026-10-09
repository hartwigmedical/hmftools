package com.hartwig.hmftools.viridian.reporting;

import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;

import org.jetbrains.annotations.Nullable;

// The reporting verdict for one oncology group plus supporting data.
// Produced for every oncology group that is detected or has integrations.
public record VirusReport(
        OncologyGroup oncologyGroup,
        // Null if not a reportable virus.
        @Nullable VirusType reportingType,
        @Nullable VirusLikelihoodType driverLikelihood,
        int integrations,
        boolean isDetected,
        // Null unless the group is detected and resolved to a representative.
        @Nullable ContigStats representativeStats,
        // Null if Purple data wasn't usable.
        @Nullable Double copiesPerTumorCell,
        VirusReportStatus status)
{
    public VirusReport
    {
        if(representativeStats != null && !isDetected)
        {
            throw new IllegalArgumentException("Resolved group must be detected: " + oncologyGroup);
        }
        if(copiesPerTumorCell != null && representativeStats == null)
        {
            throw new IllegalArgumentException("Copies per tumor cell needs a representative: " + oncologyGroup);
        }
    }

    public boolean isReported()
    {
        return status.isReported();
    }
}
