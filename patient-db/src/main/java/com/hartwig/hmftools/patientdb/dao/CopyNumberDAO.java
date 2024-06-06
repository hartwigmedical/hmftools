package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep18;
import org.jooq.Record;
import org.jooq.Result;

class CopyNumberDAO
{
    private final DSLContext context;

    CopyNumberDAO(final DSLContext context)
    {
        this.context = context;
    }

    public List<PurpleCopyNumber> read(final String sample)
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        Result<Record> result = context.select().from(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            copyNumbers.add(ImmutablePurpleCopyNumber.builder()
                    .chromosome(record.getValue(COPYNUMBER.CHROMOSOME))
                    .start(record.getValue(COPYNUMBER.START))
                    .end(record.getValue(COPYNUMBER.END))
                    .bafCount(record.getValue(COPYNUMBER.BAFCOUNT))
                    .method(CopyNumberMethod.valueOf(record.getValue(COPYNUMBER.COPYNUMBERMETHOD)))
                    .segmentStartSupport(SegmentSupport.valueOf(record.getValue(COPYNUMBER.SEGMENTSTARTSUPPORT)))
                    .segmentEndSupport(SegmentSupport.valueOf(record.getValue(COPYNUMBER.SEGMENTENDSUPPORT)))
                    .averageActualBAF(record.getValue(COPYNUMBER.BAF))
                    .averageObservedBAF(record.getValue(COPYNUMBER.OBSERVEDBAF))
                    .averageTumorCopyNumber(record.getValue(COPYNUMBER.COPYNUMBER_))
                    .depthWindowCount(record.getValue(COPYNUMBER.DEPTHWINDOWCOUNT))
                    .gcContent(record.getValue(COPYNUMBER.GCCONTENT))
                    .minStart(record.getValue(COPYNUMBER.MINSTART))
                    .maxStart(record.getValue(COPYNUMBER.MAXSTART))
                    .build());
        }

        Collections.sort(copyNumbers);
        return copyNumbers;
    }

    void writeCopyNumber(final String sample, final List<PurpleCopyNumber> copyNumbers)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        for(List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep18 inserter = context.insertInto(COPYNUMBER,
                    COPYNUMBER.SAMPLEID,
                    COPYNUMBER.CHROMOSOME,
                    COPYNUMBER.START,
                    COPYNUMBER.END,
                    COPYNUMBER.COPYNUMBERMETHOD,
                    COPYNUMBER.SEGMENTSTARTSUPPORT,
                    COPYNUMBER.SEGMENTENDSUPPORT,
                    COPYNUMBER.BAFCOUNT,
                    COPYNUMBER.OBSERVEDBAF,
                    COPYNUMBER.BAF,
                    COPYNUMBER.COPYNUMBER_,
                    COPYNUMBER.MINORALLELECOPYNUMBER,
                    COPYNUMBER.MAJORALLELECOPYNUMBER,
                    COPYNUMBER.DEPTHWINDOWCOUNT,
                    COPYNUMBER.GCCONTENT,
                    COPYNUMBER.MINSTART,
                    COPYNUMBER.MAXSTART,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberRecord(
            final Timestamp timestamp, final InsertValuesStep18 inserter, final String sample, final PurpleCopyNumber region)
    {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.method(),
                region.segmentStartSupport(),
                region.segmentEndSupport(),
                region.bafCount(),
                DatabaseUtil.decimal(region.averageObservedBAF()),
                DatabaseUtil.decimal(region.averageActualBAF()),
                DatabaseUtil.decimal(region.averageTumorCopyNumber()),
                DatabaseUtil.decimal(region.minorAlleleCopyNumber()),
                DatabaseUtil.decimal(region.majorAlleleCopyNumber()),
                region.depthWindowCount(),
                DatabaseUtil.decimal(region.gcContent()),
                region.minStart(),
                region.maxStart(),
                timestamp);
    }

    void deleteCopyNumberForSample(final String sample)
    {
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();
    }
}
