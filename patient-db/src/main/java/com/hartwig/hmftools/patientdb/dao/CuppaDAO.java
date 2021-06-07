package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.Date;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPA;

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

        String tumorLocation = cuppaResult.split("\\(")[0];
        String prediction = cuppaResult.split("\\(")[1];
        prediction = prediction.substring(0, prediction.length()-1);
        context.insertInto(CUPPA,
                CUPPA.MODIFIED,
                CUPPA.SAMPLEID,
                CUPPA.CUPPATUMORLOCATION,
                CUPPA.CUPPAPREDICTION)
                .values(timestamp, sample, tumorLocation, prediction)
                .execute();
    }

    void deleteCuppaForSample(@NotNull String sample) {
        context.delete(CUPPA).where(CUPPA.SAMPLEID.eq(sample)).execute();
    }
}
