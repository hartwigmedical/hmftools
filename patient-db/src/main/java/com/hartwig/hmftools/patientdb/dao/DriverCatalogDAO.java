package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERCATALOG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep6;
import org.jooq.InsertValuesStep7;

class DriverCatalogDAO {

    @NotNull
    private final DSLContext context;

    DriverCatalogDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull final String sample, @NotNull List<DriverCatalog> driverCatalog) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteForSample(sample);

        for (List<DriverCatalog> splitRegions : Iterables.partition(driverCatalog, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep7 inserter = context.insertInto(DRIVERCATALOG,
                    DRIVERCATALOG.SAMPLEID,
                    DRIVERCATALOG.GENE,
                    DRIVERCATALOG.CATEGORY,
                    DRIVERCATALOG.DRIVER,
                    DRIVERCATALOG.DNDSLIKELIHOOD,
                    DRIVERCATALOG.DRIVERLIKELIHOOD,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep7 inserter, @NotNull String sample,
            @NotNull DriverCatalog entry) {
        inserter.values(sample, entry.gene(), entry.category(), entry.driver(), DatabaseUtil.decimal(entry.dndsLikelihood()),  DatabaseUtil.decimal(entry.driverLikelihood()), timestamp);
    }

    void deleteForSample(@NotNull String sample) {
        context.delete(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).execute();
    }
}