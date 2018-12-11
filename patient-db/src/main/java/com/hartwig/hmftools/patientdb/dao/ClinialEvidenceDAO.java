package com.hartwig.hmftools.patientdb.dao;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class ClinialEvidenceDAO {

    @NotNull
    private final DSLContext context;

    ClinialEvidenceDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalEvidence(@NotNull String sample) {

    }
}
