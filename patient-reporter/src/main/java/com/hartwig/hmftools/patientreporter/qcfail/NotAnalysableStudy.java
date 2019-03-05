package com.hartwig.hmftools.patientreporter.qcfail;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum NotAnalysableStudy {
    CPCT("CPCT", "CPCT-02"),
    DRUP("DRUP", "DRUP-01"),
    CORE_OLD("CORE", "CORE1"),
    CORE("CORE", "CORE-01"),
    WIDE("WIDE", "WIDE-01");

    private static final String CPCT_IDENTIFIER = "CPCT";
    private static final String DRUP_IDENTIFIER = "DRUP";
    private static final String CORE_IDENTIFIER_OLD = "CORE1";
    private static final String CORE_IDENTIFIER = "CORE";
    private static final String WIDE_IDENTIFIER = "WIDE";


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
        } else if (sample.contains(CORE_IDENTIFIER)) {
            return CORE;
        } else if (sample.contains(CORE_IDENTIFIER_OLD)) {
            return CORE_OLD;
        } else if (sample.contains(WIDE_IDENTIFIER)) {
            return WIDE;
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

