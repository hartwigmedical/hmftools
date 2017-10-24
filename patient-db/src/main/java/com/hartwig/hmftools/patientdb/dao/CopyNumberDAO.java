package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumberregion.COPYNUMBERREGION;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentStatus;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep11;
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
                    .ratioSupport(true)
                    .structuralVariantSupport(StructuralVariantSupport.valueOf(record.getValue(COPYNUMBER.STRUCTURALVARIANTSUPPORT)))
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
                    .status(SegmentStatus.valueOf(record.getValue(COPYNUMBERREGION.STATUS)))
                    .ratioSupport(true)
                    .structuralVariantSupport(StructuralVariantSupport.valueOf(record.getValue(COPYNUMBERREGION.STRUCTURALVARIANTSUPPORT)))
                    .bafCount(record.getValue(COPYNUMBERREGION.BAFCOUNT))
                    .observedBAF(record.getValue(COPYNUMBERREGION.OBSERVEDBAF))
                    .observedTumorRatio(record.getValue(COPYNUMBERREGION.OBSERVEDTUMORRATIO))
                    .observedNormalRatio(record.getValue(COPYNUMBERREGION.OBSERVEDNORMALRATIO))
                    .observedTumorRatioCount(record.getValue(COPYNUMBERREGION.OBSERVEDTUMORRATIOCOUNT))
                    .modelPloidy(record.getValue(COPYNUMBERREGION.MODELPLOIDY))
                    .modelBAF(record.getValue(COPYNUMBERREGION.MODELBAF))
                    .modelTumorRatio(record.getValue(COPYNUMBERREGION.MODELTUMORRATIO))
                    .tumorCopyNumber(record.getValue(COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER))
                    .refNormalisedCopyNumber(record.getValue(COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER))
                    .cnvDeviation(record.getValue(COPYNUMBERREGION.CNVDEVIATION))
                    .deviation(record.getValue(COPYNUMBERREGION.TOTALDEVIATION))
                    .bafDeviation(record.getValue(COPYNUMBERREGION.BAFDEVIATION))
                    .broadBAF(record.getValue(COPYNUMBERREGION.HIGHCONFIDENCEBAF))
                    .broadTumorCopyNumber(record.getValue(COPYNUMBERREGION.HIGHCONFIDENCECOPYNUMBER))
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
            InsertValuesStep11 inserter = context.insertInto(COPYNUMBER,
                    COPYNUMBER.SAMPLEID,
                    COPYNUMBER.CHROMOSOME,
                    COPYNUMBER.START,
                    COPYNUMBER.END,
                    COPYNUMBER.RATIOSUPPORT,
                    COPYNUMBER.STRUCTURALVARIANTSUPPORT,
                    COPYNUMBER.BAFCOUNT,
                    COPYNUMBER.OBSERVEDBAF,
                    COPYNUMBER.ACTUALBAF,
                    COPYNUMBER.COPYNUMBER_,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }

    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStep11 inserter, String sample, PurpleCopyNumber region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.ratioSupport(),
                region.structuralVariantSupport(),
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
                    COPYNUMBERREGION.STATUS,
                    COPYNUMBERREGION.RATIOSUPPORT,
                    COPYNUMBERREGION.STRUCTURALVARIANTSUPPORT,
                    COPYNUMBERREGION.BAFCOUNT,
                    COPYNUMBERREGION.OBSERVEDBAF,
                    COPYNUMBERREGION.OBSERVEDTUMORRATIO,
                    COPYNUMBERREGION.OBSERVEDNORMALRATIO,
                    COPYNUMBERREGION.OBSERVEDTUMORRATIOCOUNT,
                    COPYNUMBERREGION.MODELPLOIDY,
                    COPYNUMBERREGION.MODELBAF,
                    COPYNUMBERREGION.MODELTUMORRATIO,
                    COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER,
                    COPYNUMBERREGION.REFNORMALISEDTUMORCOPYNUMBER,
                    COPYNUMBERREGION.CNVDEVIATION,
                    COPYNUMBERREGION.BAFDEVIATION,
                    COPYNUMBERREGION.TOTALDEVIATION,
                    COPYNUMBERREGION.HIGHCONFIDENCEBAF,
                    COPYNUMBERREGION.HIGHCONFIDENCECOPYNUMBER,
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
                region.ratioSupport(),
                region.structuralVariantSupport(),
                region.bafCount(),
                region.observedBAF(),
                region.observedTumorRatio(),
                region.observedNormalRatio(),
                region.observedTumorRatioCount(),
                region.modelPloidy(),
                region.modelBAF(),
                region.modelTumorRatio(),
                region.tumorCopyNumber(),
                region.refNormalisedCopyNumber(),
                region.cnvDeviation(),
                region.bafDeviation(),
                region.deviation(),
                region.broadBAF(),
                region.broadTumorCopyNumber(),
                region.segmentBAF(),
                region.segmentTumorCopyNumber(),
                timestamp);
    }
}
