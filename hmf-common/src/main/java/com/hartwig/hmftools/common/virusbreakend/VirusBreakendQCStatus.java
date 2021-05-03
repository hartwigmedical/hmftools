package com.hartwig.hmftools.common.virusbreakend;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum VirusBreakendQCStatus {
    LOW_VIRAL_COVERAGE,
    EXCESSIVE_VIRAL_COVERAGE,
    ASSEMBLY_DOWNSAMPLED,
    CHILD_TAXID_REFERENCE,
    UNCLEAR_TAXID_ASSIGNMENT,
    UNKNOWN;

    @NotNull
    public static VirusBreakendQCStatus extractVirusBreakendQCStatus(@NotNull String QcStatus) {
        switch (QcStatus) {
            case Strings.EMPTY:
                return UNKNOWN;
            case "LOW_VIRAL_COVERAGE":
                return LOW_VIRAL_COVERAGE;
            case "EXCESSIVE_VIRAL_COVERAGE":
                return EXCESSIVE_VIRAL_COVERAGE;
            case "ASSEMBLY_DOWNSAMPLED":
                return ASSEMBLY_DOWNSAMPLED;
            case "CHILD_TAXID_REFERENCE":
                return CHILD_TAXID_REFERENCE;
            case "UNCLEAR_TAXID_ASSIGNMENT":
                return UNCLEAR_TAXID_ASSIGNMENT;
            default:
                throw new IllegalStateException(
                        "Cannot resolve QC Status of VIRUSBreakend " + QcStatus);
        }
    }
}
