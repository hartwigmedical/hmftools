package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.AMBER;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep5;

class AmberDAO {

    @NotNull
    private final DSLContext context;

    AmberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull List<AmberBAF> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteAmberRecordsForSample(sample);

        for (List<AmberBAF> splitRegions : Iterables.partition(variants, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep5 inserter =
                    context.insertInto(AMBER, AMBER.SAMPLEID, AMBER.CHROMOSOME, AMBER.POSITION, AMBER.HETEROZYGOUS, AMBER.MODIFIED);
            splitRegions.forEach(variant -> addRecord(timestamp, inserter, sample, variant));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep5 inserter, @NotNull String sample,
            @NotNull AmberBAF variant) {
        inserter.values(sample, variant.chromosome(), variant.position(), Doubles.greaterThan(variant.normalBAF(), 0), timestamp);
    }

    void deleteAmberRecordsForSample(@NotNull String sample) {
        context.delete(AMBER).where(AMBER.SAMPLEID.eq(sample)).execute();
    }
}