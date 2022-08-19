package com.hartwig.hmftools.patientdb.dao;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLEIDMAPPING;


public class SampleIdMappingDAO {

    @NotNull
    private final DSLContext context;

    SampleIdMappingDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull String barcode) {
        deleteSampleIdMappingForSample(sample);

        context.insertInto(SAMPLEIDMAPPING,
                        SAMPLEIDMAPPING.SAMPLEID,
                        SAMPLEIDMAPPING.ISOLATIONBARCODE)
                .values(sample, barcode)
                .execute();
    }

    void deleteSampleIdMappingForSample(@NotNull String sample) {
        context.delete(SAMPLEIDMAPPING).where(SAMPLEIDMAPPING.SAMPLEID.eq(sample)).execute();
    }
}
