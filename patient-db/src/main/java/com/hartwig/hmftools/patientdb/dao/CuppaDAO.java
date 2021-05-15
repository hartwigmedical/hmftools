package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.Date;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPARESULT;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class CuppaDAO {

    @NotNull
    private final DSLContext context;

    CuppaDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeCuppa(@NotNull String sample, @NotNull String cuppaResult) {
        deleteCuppaForSample(sample);
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.insertInto(CUPPARESULT,
                CUPPARESULT.MODIFIED,
                CUPPARESULT.SAMPLEID,
                CUPPARESULT.CUPPARESULT_)
                .values(timestamp, sample, cuppaResult)
                .execute();
    }

    void deleteCuppaForSample(@NotNull String sample) {
        context.delete(CUPPARESULT).where(CUPPARESULT.SAMPLEID.eq(sample)).execute();
    }
}
