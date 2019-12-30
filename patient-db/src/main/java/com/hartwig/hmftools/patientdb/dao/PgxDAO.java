package com.hartwig.hmftools.patientdb.dao;

import com.hartwig.hmftools.common.pharmacogenetics.PGXCalls;
import com.hartwig.hmftools.common.pharmacogenetics.PGXGenotype;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class PgxDAO {

    @NotNull
    private final DSLContext context;

    PgxDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writePgx(@NotNull String sample, @NotNull PGXGenotype pgxGenotype, @NotNull PGXCalls pgxCalls) {
        deletePgxForSample(sample);
    }

    void deletePgxForSample(@NotNull String sample) {
       // context.delete().where(.eq(sample)).execute();
    }


}
