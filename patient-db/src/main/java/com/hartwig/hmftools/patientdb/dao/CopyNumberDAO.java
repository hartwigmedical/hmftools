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

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep20;
import org.jooq.InsertValuesStep9;
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
                    .averageActualBAF(record.getValue(COPYNUMBER.ACTUALBAF))
                    .averageObservedBAF(record.getValue(COPYNUMBER.OBSERVEDBAF))
                    .averageTumorCopyNumber(record.getValue(COPYNUMBER.COPYNUMBER_))
                    .build());
        }

        Collections.sort(copyNumbers);
        return copyNumbers;
    }

    void writeCopyNumber(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        for (List<PurpleCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(COPYNUMBER,
                    COPYNUMBER.SAMPLEID,
                    COPYNUMBER.CHROMOSOME,
                    COPYNUMBER.START,
                    COPYNUMBER.END,
                    COPYNUMBER.BAFCOUNT,
                    COPYNUMBER.OBSERVEDBAF,
                    COPYNUMBER.ACTUALBAF,
                    COPYNUMBER.COPYNUMBER_,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }

    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStep9 inserter, String sample, PurpleCopyNumber region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
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
            InsertValuesStep20 inserter = context.insertInto(COPYNUMBERREGION,
                    COPYNUMBERREGION.SAMPLEID,
                    COPYNUMBERREGION.CHROMOSOME,
                    COPYNUMBERREGION.START,
                    COPYNUMBERREGION.END,
                    COPYNUMBERREGION.SOURCE,
                    COPYNUMBERREGION.BAFCOUNT,
                    COPYNUMBERREGION.OBSERVEDBAF,
                    COPYNUMBERREGION.OBSERVEDTUMORRATIO,
                    COPYNUMBERREGION.OBSERVEDNORMALRATIO,
                    COPYNUMBERREGION.MODELPLOIDY,
                    COPYNUMBERREGION.MODELBAF,
                    COPYNUMBERREGION.MODELTUMORRATIO,
                    COPYNUMBERREGION.ACTUALTUMORCOPYNUMBER,
                    COPYNUMBERREGION.CNVDEVIATION,
                    COPYNUMBERREGION.BAFDEVIATION,
                    COPYNUMBERREGION.HIGHCONFIDENCEBAF,
                    COPYNUMBERREGION.HIGHCONFIDENCECOPYNUMBER,
                    COPYNUMBERREGION.FITTEDBAF,
                    COPYNUMBERREGION.FITTEDCOPYNUMBER,
                    COPYNUMBERREGION.MODIFIED);
            splitRegions.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStep20 inserter, String sample, FittedRegion region) {
        inserter.values(sample,
                region.chromosome(),
                region.start(),
                region.end(),
                region.source(),
                region.bafCount(),
                region.observedBAF(),
                region.observedTumorRatio(),
                region.observedNormalRatio(),
                region.fittedPloidy(),
                region.modelBAF(),
                region.modelTumorRatio(),
                region.tumorCopyNumber(),
                region.cnvDeviation(),
                region.bafDeviation(),
                region.broadBAF(),
                region.broadTumorCopyNumber(),
                region.segmentBAF(),
                region.segmentTumorCopyNumber(),
                timestamp);
    }
}
