package com.hartwig.hmftools.common.virus;

import org.jetbrains.annotations.NotNull;

public enum VirusBreakendQCStatus {
    NO_ABNORMALITIES,
    LOW_VIRAL_COVERAGE,
    EXCESSIVE_VIRAL_COVERAGE,
    ASSEMBLY_DOWNSAMPLED,
    CHILD_TAXID_REFERENCE,
    UNCLEAR_TAXID_ASSIGNMENT;

    @NotNull
    public static VirusBreakendQCStatus convert(@NotNull String qcStatusString) {
        if (qcStatusString.isEmpty()) {
            return NO_ABNORMALITIES;
        }

        for (VirusBreakendQCStatus qcStatus : VirusBreakendQCStatus.values()) {
            if (qcStatus.toString().equals(qcStatusString)) {
                return qcStatus;
            }
        }

        throw new IllegalStateException("Cannot resolve QC Status of VIRUSBreakend: " + qcStatusString);
    }
}
