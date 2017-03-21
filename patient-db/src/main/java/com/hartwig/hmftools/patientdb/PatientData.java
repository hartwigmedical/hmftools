package com.hartwig.hmftools.patientdb;

import org.jetbrains.annotations.Nullable;

public class PatientData {
    private final String cpctId;
    private final String drupId;
    private final String sex;
    private final Integer birthYear;
    private final String hospital;
    private final String ethnicity;

    PatientData(@Nullable String cpctId, @Nullable String drupId, @Nullable String sex, @Nullable Integer birthYear,
            @Nullable String hospital, @Nullable String ethnicity) {
        this.cpctId = cpctId;
        this.drupId = drupId;
        this.sex = sex;
        this.birthYear = birthYear;
        this.hospital = hospital;
        this.ethnicity = ethnicity;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(cpctId).append("/").append(drupId).append(": ").append(birthYear).append(" - ").append(sex).append(
                " - ").append(hospital).append(" - ").append(ethnicity).append("\n");
        return bf.toString();
    }
}
