package com.hartwig.hmftools.patientdb.dao;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class ClinicalTrialDAO {

    @NotNull
    private final DSLContext context;

    ClinicalTrialDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalTrial(@NotNull String sample) {

             
    }
}
