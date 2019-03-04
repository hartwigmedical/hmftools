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
        if (label.equals("CORE")) {
            return CORE;
        } else if (label.equals("WIDE")) {
            return WIDE;
        } else if (label.equals("CPCT")) {
            return CPCT;
        } else if (label.equals("DRUP")) {
            return DRUP;
        } else {
            return UNDEFINED;
        }
    }
}
