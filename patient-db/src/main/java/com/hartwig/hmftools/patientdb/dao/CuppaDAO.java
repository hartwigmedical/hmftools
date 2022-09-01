package com.hartwig.hmftools.patientdb.dao;

import java.time.LocalDateTime;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPA;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class CuppaDAO {

    @NotNull
    private final DSLContext context;

    CuppaDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeCuppa(@NotNull String sample, @NotNull String cancerType, double likelihood) {
        deleteCuppaForSample(sample);
        LocalDateTime timestamp = LocalDateTime.now();

        context.insertInto(CUPPA, CUPPA.MODIFIED, CUPPA.SAMPLEID, CUPPA.CUPPATUMORLOCATION, CUPPA.CUPPAPREDICTION)
                .values(timestamp, sample, cancerType, Double.toString(likelihood))
                .execute();
    }

    void deleteCuppaForSample(@NotNull String sample) {
        context.delete(CUPPA).where(CUPPA.SAMPLEID.eq(sample)).execute();
    }
}
