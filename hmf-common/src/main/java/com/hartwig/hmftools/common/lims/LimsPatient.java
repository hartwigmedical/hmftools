package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public class LimsPatient {

    @NotNull
    private final String patientId;

    public LimsPatient(@NotNull final String patientId) {
        this.patientId = patientId;
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

}
