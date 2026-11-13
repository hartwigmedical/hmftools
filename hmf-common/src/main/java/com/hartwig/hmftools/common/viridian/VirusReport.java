package com.hartwig.hmftools.common.viridian;

import org.jetbrains.annotations.Nullable;

// The reporting verdict for one virus group, plus supporting data.
// Produced for every virus group that is present in the reads or has integrations.
public record VirusReport(
        String virusName,
        String virusGroup,
        ReportQcStatus qcStatus,
        // Null if not a reportable virus.
        @Nullable VirusType reportingType,
        int integrations,
        double coveragePercentage,
        double meanDepth,
        // Null if Purple data wasn't usable.
        @Nullable Double copiesPerTumorCell,
        ClonalityStatus clonalityStatus,
        boolean isReported)
{
}
