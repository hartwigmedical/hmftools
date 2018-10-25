package com.hartwig.hmftools.patientreporter.qcfail;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotAnalysableStudy {
    CPCT("CPCT", "CPCT-02"),
    DRUP("DRUP", "DRUP-01");

    private static final String CPCT_IDENTIFIER = "CPCT";
    private static final String DRUP_IDENTIFIER = "DRUP";

    @NotNull
    private final String studyName;
    @NotNull
    private final String studyCode;

    NotAnalysableStudy(@NotNull final String studyName, @NotNull final String studyCode) {
        this.studyName = studyName;
        this.studyCode = studyCode;
    }

    @Nullable
    public static NotAnalysableStudy fromSample(@NotNull final String sample) {
        if (sample.contains(CPCT_IDENTIFIER)) {
            return CPCT;
        } else if (sample.contains(DRUP_IDENTIFIER)) {
            return DRUP;
        }
        return null;
    }

    @NotNull
    public String studyName() {
        return studyName;
    }

    @NotNull
    public String studyCode() {
        return studyCode;
    }
}

