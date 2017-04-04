package com.hartwig.hmftools.patientdb;

import org.jetbrains.annotations.Nullable;

class PatientInfo {
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

    PatientInfo(@Nullable final String cpctId, @Nullable final String drupId, @Nullable final String sex,
            @Nullable final Integer birthYear, @Nullable final String hospital, @Nullable final String ethnicity) {
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

    String cpctId() {
        return cpctId;
    }

    String drupId() {
        return drupId;
    }

    String sex() {
        return sex;
    }

    Integer birthYear() {
        return birthYear;
    }

    String hospital() {
        return hospital;
    }

    String ethnicity() {
        return ethnicity;
    }
}
