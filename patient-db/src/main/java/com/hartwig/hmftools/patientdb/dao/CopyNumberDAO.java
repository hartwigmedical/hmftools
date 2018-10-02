package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
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
import org.jooq.InsertValuesStep14;
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
                    .depthWindowCount(record.getValue(COPYNUMBER.DEPTHWINDOWCOUNT))
                    .gcContent(record.getValue(COPYNUMBER.GCCONTENT))
                    .build());
        }

        Collections.sort(copyNumbers);
        return copyNumbers;
    }

    void writeCopyNumber(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        for (List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep14 inserter = context.insertInto(COPYNUMBER,
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
                    COPYNUMBER.DEPTHWINDOWCOUNT,
                    COPYNUMBER.GCCONTENT,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    void writeGermlineCopyNumber(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBERGERMLINE).where(COPYNUMBERGERMLINE.SAMPLEID.eq(sample)).execute();

        for (List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep14 inserter = context.insertInto(COPYNUMBERGERMLINE,
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
                    COPYNUMBERGERMLINE.DEPTHWINDOWCOUNT,
                    COPYNUMBERGERMLINE.GCCONTENT,
                    COPYNUMBERGERMLINE.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep14 inserter, @NotNull String sample,
            @NotNull PurpleCopyNumber region) {
        //noinspection unchecked
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
                region.depthWindowCount(),
                DatabaseUtil.decimal(region.gcContent()),
                timestamp);
    }

    void writeCopyNumberRegions(@NotNull final String sample, @NotNull List<FittedRegion> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBERREGION).where(COPYNUMBERREGION.SAMPLEID.eq(sample)).execute();

        for (List<FittedRegion> splitRegions : Iterables.partition(regions, DB_BATCH_INSERT_SIZE)) {
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
                    COPYNUMBERREGION.DEPTHWINDOWCOUNT,
                    COPYNUMBERREGION.GCCONTENT,
                    COPYNUMBERREGION.MINORALLELEPLOIDY,
                    COPYNUMBERREGION.MAJORALLELEPLOIDY,
                    COPYNUMBERREGION.ACTUALTUMORBAF,
                    COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER,
                    COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER,
                    COPYNUMBERREGION.MINORALLELEPLOIDYDEVIATION,
                    COPYNUMBERREGION.MAJORALLELEPLOIDYDEVIATION,
                    COPYNUMBERREGION.TOTALDEVIATION,
                    COPYNUMBERREGION.PLOIDYPENALTY,
                    COPYNUMBERREGION.FITTEDBAF,
                    COPYNUMBERREGION.FITTEDCOPYNUMBER,
                    COPYNUMBERREGION.DEPTHWINDOWCOUNT,
                    COPYNUMBERREGION.MINSTART,
                    COPYNUMBERREGION.MAXSTART,
                    COPYNUMBERREGION.MODIFIED);
            splitRegions.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull FittedRegion region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.status(),
                region.svCluster(),
                region.ratioSupport(),
                region.support(),
                region.bafCount(),
                DatabaseUtil.decimal(region.observedBAF()),
                DatabaseUtil.decimal(region.observedTumorRatio()),
                DatabaseUtil.decimal(region.observedNormalRatio()),
                DatabaseUtil.decimal(region.depthWindowCount()),
                DatabaseUtil.decimal(region.gcContent()),
                DatabaseUtil.decimal(region.minorAllelePloidy()),
                DatabaseUtil.decimal(region.majorAllelePloidy()),
                DatabaseUtil.decimal(region.tumorBAF()),
                DatabaseUtil.decimal(region.tumorCopyNumber()),
                DatabaseUtil.decimal(region.refNormalisedCopyNumber()),
                DatabaseUtil.decimal(region.minorAllelePloidyDeviation()),
                DatabaseUtil.decimal(region.majorAllelePloidyDeviation()),
                DatabaseUtil.decimal(region.deviation()),
                DatabaseUtil.decimal(region.ploidyPenalty()),
                DatabaseUtil.decimal(region.fittedBAF()),
                DatabaseUtil.decimal(region.fittedTumorCopyNumber()),
                region.depthWindowCount(),
                region.minStart(),
                region.maxStart(),
                timestamp);
    }

    @NotNull
    List<FittedRegion> readCopyNumberRegions(@NotNull final String sample) {
        List<FittedRegion> fittedRegions = Lists.newArrayList();

        Result<Record> result = context.select().from(COPYNUMBERREGION).where(COPYNUMBERREGION.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            fittedRegions.add(ImmutableFittedRegion.builder()
                    .chromosome(record.getValue(COPYNUMBERREGION.CHROMOSOME))
                    .start(record.getValue(COPYNUMBERREGION.START))
                    .end(record.getValue(COPYNUMBERREGION.END))
                    .status(GermlineStatus.valueOf(record.getValue(COPYNUMBERREGION.GERMLINESTATUS)))
                    .svCluster(record.getValue(COPYNUMBERREGION.SVCLUSTER) != 0)
                    .ratioSupport(record.getValue(COPYNUMBERREGION.RATIOSUPPORT) != 0)
                    .support(SegmentSupport.valueOf(record.getValue(COPYNUMBERREGION.SEGMENTSTARTSUPPORT)))
                    .bafCount(record.getValue(COPYNUMBERREGION.BAFCOUNT))
                    .observedBAF(record.getValue(COPYNUMBERREGION.OBSERVEDBAF))
                    .observedTumorRatio(record.getValue(COPYNUMBERREGION.OBSERVEDTUMORRATIO))
                    .observedNormalRatio(record.getValue(COPYNUMBERREGION.OBSERVEDNORMALRATIO))
                    .depthWindowCount(record.getValue(COPYNUMBERREGION.DEPTHWINDOWCOUNT))
                    .gcContent(record.getValue(COPYNUMBERREGION.GCCONTENT))
                    .minorAllelePloidy(record.getValue(COPYNUMBERREGION.MINORALLELEPLOIDY))
                    .majorAllelePloidy(record.getValue(COPYNUMBERREGION.MAJORALLELEPLOIDY))
                    .tumorBAF(record.getValue(COPYNUMBERREGION.ACTUALTUMORBAF))
                    .tumorCopyNumber(record.getValue(COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER))
                    .refNormalisedCopyNumber(record.getValue(COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER))
                    .minorAllelePloidyDeviation(record.getValue(COPYNUMBERREGION.MINORALLELEPLOIDYDEVIATION))
                    .majorAllelePloidyDeviation(record.getValue(COPYNUMBERREGION.MAJORALLELEPLOIDYDEVIATION))
                    .deviation(record.getValue(COPYNUMBERREGION.TOTALDEVIATION))
                    .fittedBAF(record.getValue(COPYNUMBERREGION.FITTEDBAF))
                    .fittedTumorCopyNumber(record.getValue(COPYNUMBERREGION.FITTEDCOPYNUMBER))
                    .ploidyPenalty(record.getValue(COPYNUMBERREGION.PLOIDYPENALTY))
                    .depthWindowCount(record.getValue(COPYNUMBERREGION.DEPTHWINDOWCOUNT))
                    .minStart(record.getValue(COPYNUMBERREGION.MINSTART))
                    .maxStart(record.getValue(COPYNUMBERREGION.MAXSTART))
                    .build());
        }

        Collections.sort(fittedRegions);
        return fittedRegions;
    }

    void deleteCopyNumberForSample(@NotNull String sample) {
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();
        context.delete(COPYNUMBERGERMLINE).where(COPYNUMBERGERMLINE.SAMPLEID.eq(sample)).execute();
        context.delete(COPYNUMBERREGION).where(COPYNUMBERREGION.SAMPLEID.eq(sample)).execute();
    }
}
