package com.hartwig.hmftools.patientreporter.qcfail;

import com.hartwig.hmftools.common.lims.LimsStudy;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailStudy {
    CPCT("CPCT", "CPCT-02"),
    DRUP("DRUP", "DRUP-01"),
    CORE("CORE", "CORE"),
    WIDE("WIDE", "WIDE");

    @NotNull
    private final String studyName;
    @NotNull
    private final String studyCode;

    QCFailStudy(@NotNull final String studyName, @NotNull final String studyCode) {
        this.studyName = studyName;
        this.studyCode = studyCode;
    }

    @Nullable
    public static QCFailStudy fromSample(@NotNull final String sample) {
        LimsStudy type = LimsStudy.fromSampleId(sample);

        switch (type) {
            case CPCT:
                return CPCT;
            case DRUP:
                return DRUP;
            case CORE:
                return CORE;
            case WIDE:
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

