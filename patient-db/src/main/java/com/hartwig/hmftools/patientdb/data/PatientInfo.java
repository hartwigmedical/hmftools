package com.hartwig.hmftools.patientdb.data;

import org.jetbrains.annotations.Nullable;

public class PatientInfo {
    @Nullable
    private final String cpctId;
    @Nullable
    private final String drupId;
    @Nullable
    private final String sex;
    @Nullable
    private final Integer birthYear;
    @Nullable
    private final String hospital;
    @Nullable
    private final String ethnicity;

    public PatientInfo(@Nullable final String cpctId, @Nullable final String drupId, @Nullable final String sex,
            @Nullable final Integer birthYear, @Nullable final String hospital, @Nullable final String ethnicity) {
        this.cpctId = cpctId;
        this.drupId = drupId;
        this.sex = sex;
        this.birthYear = birthYear;
        this.hospital = hospital;
        this.ethnicity = ethnicity;
    }

    public String cpctId() {
        return cpctId;
    }

    String drupId() {
        return drupId;
    }

    public String sex() {
        return sex;
    }

    public Integer birthYear() {
        return birthYear;
    }

    public String hospital() {
        return hospital;
    }

    public String ethnicity() {
        return ethnicity;
    }
}
