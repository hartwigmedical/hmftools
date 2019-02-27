package com.hartwig.hmftools.common.lims;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class LimsPatient {

    @NotNull
    private final List<String> patients;

    public LimsPatient(@NotNull final List<String> patients) {
        this.patients = patients;
    }

    @NotNull
    public List<String> patients() {
        return patients;
    }

}
