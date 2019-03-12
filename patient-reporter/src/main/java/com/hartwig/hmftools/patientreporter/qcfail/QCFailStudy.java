package com.hartwig.hmftools.patientreporter.qcfail;

import com.hartwig.hmftools.common.lims.LimsSampleType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailStudy {
    CPCT("CPCT", "CPCT-02"),
    DRUP("DRUP", "DRUP-01"),
    CORE("CORE", "CORE"),
    WIDE("WIDE", "WIDE-01");

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
        LimsSampleType type = LimsSampleType.fromSampleId(sample);

        if (type == LimsSampleType.CPCT) {
            return CPCT;
        } else if (type == LimsSampleType.DRUP) {
            return DRUP;
        } else if (type == LimsSampleType.CORE) {
            return CORE;
        } else if (type == LimsSampleType.WIDE) {
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

