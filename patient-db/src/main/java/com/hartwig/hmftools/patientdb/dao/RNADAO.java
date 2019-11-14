package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.RNA;

import java.util.Set;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class RNADAO {

    @NotNull
    private final DSLContext context;

    RNADAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull Set<String> samples) {
        context.truncate(RNA).execute();

        for (String sample : samples) {
            context.insertInto(RNA, RNA.SAMPLEID).values(sample).execute();
        }
    }
}
