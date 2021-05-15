package com.hartwig.hmftools.common.virusbreakend;

import org.jetbrains.annotations.NotNull;

public enum VirusBreakendQCStatus {
    LOW_VIRAL_COVERAGE,
    EXCESSIVE_VIRAL_COVERAGE,
    ASSEMBLY_DOWNSAMPLED,
    CHILD_TAXID_REFERENCE,
    UNCLEAR_TAXID_ASSIGNMENT,
    UNKNOWN;

    @NotNull
    public static VirusBreakendQCStatus extractVirusBreakendQCStatus(@NotNull String qcStatusString) {
        if (qcStatusString.isEmpty()) {
            return UNKNOWN;
        }

        for (VirusBreakendQCStatus qcStatus : VirusBreakendQCStatus.values()) {
            if (qcStatus.toString().equals(qcStatusString)) {
                return qcStatus;
            }
        }

        throw new IllegalStateException("Cannot resolve QC Status of VIRUSBreakend: " + qcStatusString);
    }
}
