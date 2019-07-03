package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERCATALOG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.Record;
import org.jooq.Result;

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
            InsertValuesStep13 inserter = context.insertInto(DRIVERCATALOG,
                    DRIVERCATALOG.SAMPLEID,
                    DRIVERCATALOG.GENE,
                    DRIVERCATALOG.CATEGORY,
                    DRIVERCATALOG.DRIVER,
                    DRIVERCATALOG.DNDSLIKELIHOOD,
                    DRIVERCATALOG.DRIVERLIKELIHOOD,
                    DRIVERCATALOG.MISSENSE,
                    DRIVERCATALOG.NONSENSE,
                    DRIVERCATALOG.SPLICE,
                    DRIVERCATALOG.INFRAME,
                    DRIVERCATALOG.FRAMESHIFT,
                    DRIVERCATALOG.BIALLELIC,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep13 inserter, @NotNull String sample,
            @NotNull DriverCatalog entry) {
        //noinspection unchecked
        inserter.values(sample,
                entry.gene(),
                entry.category(),
                entry.driver(),
                DatabaseUtil.decimal(entry.dndsLikelihood()),
                DatabaseUtil.decimal(entry.driverLikelihood()),
                entry.missense(),
                entry.nonsense(),
                entry.splice(),
                entry.inframe(),
                entry.frameshift(),
                entry.biallelic(),
                timestamp);
    }

    void deleteForSample(@NotNull String sample) {
        context.delete(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).execute();
    }

    @NotNull
    List<DriverCatalog> readDriverData(@NotNull final String sample) {
        final List<DriverCatalog> dcList = Lists.newArrayList();

        final Result<Record> result = context.select().from(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            final DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                    .gene(record.getValue(DRIVERCATALOG.GENE))
                    .category(DriverCategory.valueOf(record.getValue(DRIVERCATALOG.CATEGORY)))
                    .driver(DriverType.valueOf(record.getValue(DRIVERCATALOG.DRIVER)))
                    .driverLikelihood(record.getValue(DRIVERCATALOG.DRIVERLIKELIHOOD))
                    .dndsLikelihood(record.getValue(DRIVERCATALOG.DNDSLIKELIHOOD))
                    .missense(record.getValue(DRIVERCATALOG.MISSENSE))
                    .nonsense(record.getValue(DRIVERCATALOG.NONSENSE))
                    .splice(record.getValue(DRIVERCATALOG.SPLICE))
                    .inframe(record.getValue(DRIVERCATALOG.INFRAME))
                    .frameshift(record.getValue(DRIVERCATALOG.FRAMESHIFT))
                    .biallelic(record.getValue(DRIVERCATALOG.BIALLELIC) != 0)
                    .build();

            dcList.add(driverCatalog);
        }

        return dcList;
    }
}