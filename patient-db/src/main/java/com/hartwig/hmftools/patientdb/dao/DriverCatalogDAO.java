package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_PURPLE_GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_PURPLE_SOMATIC;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERCATALOG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Collection;
import java.util.Date;
import java.util.EnumSet;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;
import org.jooq.Record;
import org.jooq.Result;

class DriverCatalogDAO
{
    private final DSLContext context;

    DriverCatalogDAO(final DSLContext context)
    {
        this.context = context;
    }

    void writePurpleDrivers(final String sample, @Nullable List<DriverCatalog> somaticCatalog, @Nullable List<DriverCatalog> germlineCatalog)
    {
        if(somaticCatalog != null)
            write(sample, somaticCatalog, Sets.newHashSet(DRIVERS_PURPLE_SOMATIC));

        if(germlineCatalog != null)
            write(sample, germlineCatalog, Sets.newHashSet(DRIVERS_PURPLE_GERMLINE));
    }

    void writeLinxDrivers(final String sample, final List<DriverCatalog> driverCatalog, final EnumSet<DriverType> driverTypes)
    {
        write(sample, driverCatalog, driverTypes);
    }

    void write(final String sample, final List<DriverCatalog> driverCatalog, final Collection<DriverType> types)
    {
        final List<DriverCatalog> filtered = driverCatalog.stream().filter(x -> types.contains(x.driver())).collect(Collectors.toList());
        deleteForSample(sample, types);
        insert(sample, filtered);
    }

    private void insert(final String sample, final List<DriverCatalog> driverCatalog)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        for(List<DriverCatalog> splitRegions : Iterables.partition(driverCatalog, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep19 inserter = context.insertInto(DRIVERCATALOG,
                    DRIVERCATALOG.SAMPLEID,
                    DRIVERCATALOG.CHROMOSOME,
                    DRIVERCATALOG.CHROMOSOMEBAND,
                    DRIVERCATALOG.GENE,
                    DRIVERCATALOG.TRANSCRIPTID,
                    DRIVERCATALOG.CANONICALTRANSCRIPT,
                    DRIVERCATALOG.DRIVER,
                    DRIVERCATALOG.CATEGORY,
                    DRIVERCATALOG.LIKELIHOODMETHOD,
                    DRIVERCATALOG.DRIVERLIKELIHOOD,
                    DRIVERCATALOG.MISSENSE,
                    DRIVERCATALOG.NONSENSE,
                    DRIVERCATALOG.SPLICE,
                    DRIVERCATALOG.INFRAME,
                    DRIVERCATALOG.FRAMESHIFT,
                    DRIVERCATALOG.BIALLELIC,
                    DRIVERCATALOG.MINCOPYNUMBER,
                    DRIVERCATALOG.MAXCOPYNUMBER,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addRecord(
            Timestamp timestamp, InsertValuesStep19 inserter, String sample, DriverCatalog entry)
    {
        inserter.values(
                sample,
                entry.chromosome(),
                entry.chromosomeBand(),
                entry.gene(),
                entry.transcript(),
                entry.isCanonical(),
                entry.driver(),
                entry.category(),
                entry.likelihoodMethod(),
                DatabaseUtil.decimal(entry.driverLikelihood()),
                entry.missense(),
                entry.nonsense(),
                entry.splice(),
                entry.inframe(),
                entry.frameshift(),
                entry.biallelic(),
                DatabaseUtil.decimal(entry.minCopyNumber()),
                DatabaseUtil.decimal(entry.maxCopyNumber()),
                timestamp);
    }

    void deleteForSample(String sample)
    {
        context.delete(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).execute();
    }

    void deleteForSample(String sample, Collection<DriverType> types)
    {
        final List<String> stringTypes = types.stream().map(Enum::toString).collect(Collectors.toList());
        context.delete(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).and(DRIVERCATALOG.DRIVER.in(stringTypes)).execute();
    }

    @NotNull
    List<DriverCatalog> readDriverData(String sample)
    {
        List<DriverCatalog> dcList = Lists.newArrayList();

        Result<Record> result = context.select().from(DRIVERCATALOG).where(DRIVERCATALOG.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                    .gene(record.getValue(DRIVERCATALOG.GENE))
                    .chromosome(record.getValue(DRIVERCATALOG.CHROMOSOME))
                    .chromosomeBand(record.getValue(DRIVERCATALOG.CHROMOSOMEBAND))
                    .driver(DriverType.checkConvertType(record.getValue(DRIVERCATALOG.DRIVER)))
                    .category(DriverCategory.valueOf(record.getValue(DRIVERCATALOG.CATEGORY)))
                    .transcript(record.getValue(DRIVERCATALOG.TRANSCRIPTID))
                    .isCanonical(record.getValue(DRIVERCATALOG.CANONICALTRANSCRIPT) != 0)
                    .likelihoodMethod(LikelihoodMethod.valueOf(record.getValue(DRIVERCATALOG.LIKELIHOODMETHOD)))
                    .driverLikelihood(record.getValue(DRIVERCATALOG.DRIVERLIKELIHOOD))
                    .missense(record.getValue(DRIVERCATALOG.MISSENSE))
                    .nonsense(record.getValue(DRIVERCATALOG.NONSENSE))
                    .splice(record.getValue(DRIVERCATALOG.SPLICE))
                    .inframe(record.getValue(DRIVERCATALOG.INFRAME))
                    .frameshift(record.getValue(DRIVERCATALOG.FRAMESHIFT))
                    .biallelic(record.getValue(DRIVERCATALOG.BIALLELIC) != 0)
                    .minCopyNumber(record.getValue(DRIVERCATALOG.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(DRIVERCATALOG.MAXCOPYNUMBER))
                    .build();

            dcList.add(driverCatalog);
        }

        return dcList;
    }
}