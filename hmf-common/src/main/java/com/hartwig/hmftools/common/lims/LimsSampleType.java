package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum LimsSampleType {
    CORE("CORE"),
    WIDE("WIDE"),
    CPCT("CPCT"),
    DRUP("DRUP"),
    UNDEFINED(Strings.EMPTY);

    @NotNull
    private final String studyType;

    LimsSampleType(@NotNull final String studyType) {
        this.studyType = studyType;
    }

    @NotNull
   public String labelSample() {
        return studyType;
    }

    @NotNull
    public static LimsSampleType fromSample(@NotNull final String label) {

        switch (label) {
            case "CORE": return CORE;
            case "WIDE": return WIDE;
            case "CPCT": return CPCT;
            case "DRUP": return DRUP;
            default: return UNDEFINED;
        }
    }
}
