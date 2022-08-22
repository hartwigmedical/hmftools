package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SNPCHECK;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class SnpCheckDAO {

    @NotNull
    private final DSLContext context;

    SnpCheckDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull String isolationBarcode, boolean isPass) {
        deleteSnpCheckForSample(sample);
        context.insertInto(SNPCHECK, SNPCHECK.SAMPLEID, SNPCHECK.ISOLATIONBARCODE, SNPCHECK.ISPASS)
                .values(sample, isolationBarcode, isPass ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteSnpCheckForSample(@NotNull String sample) {
        context.delete(SNPCHECK).where(SNPCHECK.SAMPLEID.eq(sample)).execute();
    }
}
