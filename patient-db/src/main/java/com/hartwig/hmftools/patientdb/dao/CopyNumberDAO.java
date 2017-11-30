package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.COPYNUMBERGERMLINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumberregion.COPYNUMBERREGION;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep12;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Result;

class CopyNumberDAO {

    @NotNull
    private final DSLContext context;

    CopyNumberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public List<PurpleCopyNumber> read(@NotNull final String sample) {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        Result<Record> result = context.select().from(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            copyNumbers.add(ImmutablePurpleCopyNumber.builder()
                    .chromosome(record.getValue(COPYNUMBER.CHROMOSOME))
                    .start(record.getValue(COPYNUMBER.START))
                    .end(record.getValue(COPYNUMBER.END))
                    .bafCount(record.getValue(COPYNUMBER.BAFCOUNT))
                    .method(CopyNumberMethod.valueOf(record.getValue(COPYNUMBER.COPYNUMBERMETHOD)))
                    .segmentStartSupport(SegmentSupport.valueOf(record.getValue(COPYNUMBER.SEGMENTSTARTSUPPORT)))
                    .segmentEndSupport(SegmentSupport.valueOf(record.getValue(COPYNUMBER.SEGMENTENDSUPPORT)))
                    .averageActualBAF(record.getValue(COPYNUMBER.ACTUALBAF))
                    .averageObservedBAF(record.getValue(COPYNUMBER.OBSERVEDBAF))
                    .averageTumorCopyNumber(record.getValue(COPYNUMBER.COPYNUMBER_))
                    .build());
        }

        Collections.sort(copyNumbers);
        return copyNumbers;
    }

    @NotNull
    public List<FittedRegion> readCopyNumberRegions(@NotNull final String sample) {
        List<FittedRegion> results = Lists.newArrayList();

        Result<Record> query = context.select().from(COPYNUMBERREGION).where(COPYNUMBERREGION.SAMPLEID.eq(sample)).fetch();

        for (Record record : query) {
            results.add(ImmutableFittedRegion.builder()
                    .chromosome(record.getValue(COPYNUMBERREGION.CHROMOSOME))
                    .start(record.getValue(COPYNUMBERREGION.START))
                    .end(record.getValue(COPYNUMBERREGION.END))
                    .status(GermlineStatus.valueOf(record.getValue(COPYNUMBERREGION.GERMLINESTATUS)))
                    .svCluster(false)
                    .ratioSupport(true)
                    .support(SegmentSupport.valueOf(record.getValue(COPYNUMBERREGION.SEGMENTSTARTSUPPORT)))
                    .bafCount(record.getValue(COPYNUMBERREGION.BAFCOUNT))
                    .observedBAF(record.getValue(COPYNUMBERREGION.OBSERVEDBAF))
                    .observedTumorRatio(record.getValue(COPYNUMBERREGION.OBSERVEDTUMORRATIO))
                    .observedNormalRatio(record.getValue(COPYNUMBERREGION.OBSERVEDNORMALRATIO))
                    .observedTumorRatioCount(record.getValue(COPYNUMBERREGION.OBSERVEDTUMORRATIOCOUNT))
                    .gcContent(record.getValue(COPYNUMBERREGION.GCCONTENT))
                    .modelPloidy(record.getValue(COPYNUMBERREGION.MODELPLOIDY))
                    .modelBAF(record.getValue(COPYNUMBERREGION.MODELBAF))
                    .modelTumorRatio(record.getValue(COPYNUMBERREGION.MODELTUMORRATIO))
                    .tumorBAF(record.getValue(COPYNUMBERREGION.ACTUALTUMORBAF))
                    .tumorCopyNumber(record.getValue(COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER))
                    .refNormalisedCopyNumber(record.getValue(COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER))
                    .cnvDeviation(record.getValue(COPYNUMBERREGION.CNVDEVIATION))
                    .deviation(record.getValue(COPYNUMBERREGION.TOTALDEVIATION))
                    .bafDeviation(record.getValue(COPYNUMBERREGION.BAFDEVIATION))
                    .segmentBAF(record.getValue(COPYNUMBERREGION.FITTEDBAF))
                    .segmentTumorCopyNumber(record.getValue(COPYNUMBERREGION.FITTEDCOPYNUMBER))
                    .build());
        }

        Collections.sort(results);
        return results;
    }

    void writeCopyNumber(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        for (List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, BATCH_INSERT_SIZE)) {
            InsertValuesStep12 inserter = context.insertInto(COPYNUMBER,
                    COPYNUMBER.SAMPLEID,
                    COPYNUMBER.CHROMOSOME,
                    COPYNUMBER.START,
                    COPYNUMBER.END,
                    COPYNUMBER.COPYNUMBERMETHOD,
                    COPYNUMBER.SEGMENTSTARTSUPPORT,
                    COPYNUMBER.SEGMENTENDSUPPORT,
                    COPYNUMBER.BAFCOUNT,
                    COPYNUMBER.OBSERVEDBAF,
                    COPYNUMBER.ACTUALBAF,
                    COPYNUMBER.COPYNUMBER_,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    void writeGermlineCopyNumber(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBERGERMLINE).where(COPYNUMBERGERMLINE.SAMPLEID.eq(sample)).execute();

        for (List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, BATCH_INSERT_SIZE)) {
            InsertValuesStep12 inserter = context.insertInto(COPYNUMBERGERMLINE,
                    COPYNUMBERGERMLINE.SAMPLEID,
                    COPYNUMBERGERMLINE.CHROMOSOME,
                    COPYNUMBERGERMLINE.START,
                    COPYNUMBERGERMLINE.END,
                    COPYNUMBERGERMLINE.COPYNUMBERMETHOD,
                    COPYNUMBERGERMLINE.SEGMENTSTARTSUPPORT,
                    COPYNUMBERGERMLINE.SEGMENTENDSUPPORT,
                    COPYNUMBERGERMLINE.BAFCOUNT,
                    COPYNUMBERGERMLINE.OBSERVEDBAF,
                    COPYNUMBERGERMLINE.ACTUALBAF,
                    COPYNUMBERGERMLINE.COPYNUMBER,
                    COPYNUMBERGERMLINE.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStep12 inserter, String sample, PurpleCopyNumber region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.method(),
                region.segmentStartSupport(),
                region.segmentEndSupport(),
                region.bafCount(),
                region.averageObservedBAF(),
                region.averageActualBAF(),
                region.averageTumorCopyNumber(),
                timestamp);
    }

    void writeCopyNumberRegions(@NotNull final String sample, @NotNull List<FittedRegion> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBERREGION).where(COPYNUMBERREGION.SAMPLEID.eq(sample)).execute();

        for (List<FittedRegion> splitRegions : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(COPYNUMBERREGION,
                    COPYNUMBERREGION.SAMPLEID,
                    COPYNUMBERREGION.CHROMOSOME,
                    COPYNUMBERREGION.START,
                    COPYNUMBERREGION.END,
                    COPYNUMBERREGION.GERMLINESTATUS,
                    COPYNUMBERREGION.SVCLUSTER,
                    COPYNUMBERREGION.RATIOSUPPORT,
                    COPYNUMBERREGION.SEGMENTSTARTSUPPORT,
                    COPYNUMBERREGION.BAFCOUNT,
                    COPYNUMBERREGION.OBSERVEDBAF,
                    COPYNUMBERREGION.OBSERVEDTUMORRATIO,
                    COPYNUMBERREGION.OBSERVEDNORMALRATIO,
                    COPYNUMBERREGION.OBSERVEDTUMORRATIOCOUNT,
                    COPYNUMBERREGION.GCCONTENT,
                    COPYNUMBERREGION.MODELPLOIDY,
                    COPYNUMBERREGION.MODELBAF,
                    COPYNUMBERREGION.MODELTUMORRATIO,
                    COPYNUMBERREGION.ACTUALTUMORBAF,
                    COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER,
                    COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER,
                    COPYNUMBERREGION.CNVDEVIATION,
                    COPYNUMBERREGION.BAFDEVIATION,
                    COPYNUMBERREGION.TOTALDEVIATION,
                    COPYNUMBERREGION.FITTEDBAF,
                    COPYNUMBERREGION.FITTEDCOPYNUMBER,
                    COPYNUMBERREGION.MODIFIED);
            splitRegions.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStepN inserter, String sample, FittedRegion region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.status(),
                region.svCluster(),
                region.ratioSupport(),
                region.support(),
                region.bafCount(),
                region.observedBAF(),
                region.observedTumorRatio(),
                region.observedNormalRatio(),
                region.observedTumorRatioCount(),
                region.gcContent(),
                region.modelPloidy(),
                region.modelBAF(),
                region.modelTumorRatio(),
                region.tumorBAF(),
                region.tumorCopyNumber(),
                region.refNormalisedCopyNumber(),
                region.cnvDeviation(),
                region.bafDeviation(),
                region.deviation(),
                region.segmentBAF(),
                region.segmentTumorCopyNumber(),
                timestamp);
    }
}
