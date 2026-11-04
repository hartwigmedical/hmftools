package com.hartwig.hmftools.viridian.reporting;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.reference.Virus;

import org.jetbrains.annotations.Nullable;

// The reporting verdict for one virus group, plus supporting data.
// Produced for every virus group that is present in the reads or has integrations.
public record VirusReport(
        Virus virus,
        ReportQcStatus qcStatus,
        // Null if not a reportable virus.
        @Nullable VirusType reportingType,
        int integrations,
        @Nullable VirusProfile profile,
        // Null if Purple data wasn't usable.
        @Nullable Double copiesPerTumorCell,
        ClonalityStatus clonalityStatus,
        boolean isReported)
{
}
