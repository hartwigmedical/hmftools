package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SNPCHECK;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class SnpCheckDAO
{
    private final DSLContext context;

    SnpCheckDAO(final DSLContext context) {
        this.context = context;
    }

    void write(final String sample, boolean isPass)
    {
        deleteSnpCheckForSample(sample);
        context.insertInto(SNPCHECK, SNPCHECK.SAMPLEID, SNPCHECK.ISPASS).values(sample, isPass ? (byte) 1 : (byte) 0).execute();
    }

    void deleteSnpCheckForSample(final String sample) {
        context.delete(SNPCHECK).where(SNPCHECK.SAMPLEID.eq(sample)).execute();
    }
}
