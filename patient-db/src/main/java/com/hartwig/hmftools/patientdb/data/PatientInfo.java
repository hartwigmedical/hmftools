package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

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
    @Nullable
    private final LocalDate deathDate;

    public PatientInfo(@Nullable final String cpctId, @Nullable final String drupId, @Nullable final String sex,
            @Nullable final Integer birthYear, @Nullable final String hospital, @Nullable final String ethnicity,
            @Nullable final LocalDate deathDate) {
        this.cpctId = cpctId;
        this.drupId = drupId;
        this.sex = sex;
        this.birthYear = birthYear;
        this.hospital = hospital;
        this.ethnicity = ethnicity;
        this.deathDate = deathDate;
    }

    @Nullable
    public String cpctId() {
        return cpctId;
    }

    @Nullable
    String drupId() {
        return drupId;
    }

    @Nullable
    public String sex() {
        return sex;
    }

    @Nullable
    public Integer birthYear() {
        return birthYear;
    }

    @Nullable
    public String hospital() {
        return hospital;
    }

    @Nullable
    public String ethnicity() {
        return ethnicity;
    }

    @Nullable
    public LocalDate deathDate() {
        return deathDate;
    }
}
