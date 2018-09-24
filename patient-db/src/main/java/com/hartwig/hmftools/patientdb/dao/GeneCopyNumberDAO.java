package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Result;

class GeneCopyNumberDAO {

    @NotNull
    private final DSLContext context;

    GeneCopyNumberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public final List<GeneCopyNumber> read(@NotNull final String sample) {
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        final Result<Record> result = context.select().from(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            geneCopyNumbers.add(ImmutableGeneCopyNumber.builder()
                    .chromosome(String.valueOf(record.getValue(GENECOPYNUMBER.CHROMOSOME)))
                    .start(record.getValue(GENECOPYNUMBER.START))
                    .end(record.getValue(GENECOPYNUMBER.END))
                    .gene(record.getValue(GENECOPYNUMBER.GENE))
                    .minCopyNumber(record.getValue(GENECOPYNUMBER.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(GENECOPYNUMBER.MAXCOPYNUMBER))
                    .somaticRegions(record.getValue(GENECOPYNUMBER.SOMATICREGIONS))
                    .germlineHomRegions(record.getValue(GENECOPYNUMBER.GERMLINEHOMREGIONS))
                    .germlineHet2HomRegions(record.getValue(GENECOPYNUMBER.GERMLINEHETREGIONS))
                    .transcriptID(record.getValue(GENECOPYNUMBER.TRANSCRIPTID))
                    .transcriptVersion(record.getValue(GENECOPYNUMBER.TRANSCRIPTVERSION))
                    .chromosomeBand(record.getValue(GENECOPYNUMBER.CHROMOSOMEBAND))
                    .minRegions(record.getValue(GENECOPYNUMBER.MINREGIONS))
                    .minRegionStart(record.getValue(GENECOPYNUMBER.MINREGIONSTART))
                    .minRegionEnd(record.getValue(GENECOPYNUMBER.MINREGIONEND))
                    .minRegionStartSupport(SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONSTARTSUPPORT)))
                    .minRegionEndSupport(SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONENDSUPPORT)))
                    .minRegionMethod(CopyNumberMethod.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONMETHOD)))
                    .nonsenseBiallelicCount(record.getValue(GENECOPYNUMBER.NONSENSEBIALLELICVARIANTS))
                    .nonsenseNonBiallelicCount(record.getValue(GENECOPYNUMBER.NONSENSENONBIALLELICVARIANTS))
                    .nonsenseNonBiallelicPloidy(record.getValue(GENECOPYNUMBER.NONSENSENONBIALLELICPLOIDY))
                    .spliceBiallelicCount(record.getValue(GENECOPYNUMBER.SPLICEBIALLELICVARIANTS))
                    .spliceNonBiallelicCount(record.getValue(GENECOPYNUMBER.SPLICENONBIALLELICVARIANTS))
                    .spliceNonBiallelicPloidy(record.getValue(GENECOPYNUMBER.SPLICENONBIALLELICPLOIDY))
                    .missenseBiallelicCount(record.getValue(GENECOPYNUMBER.MISSENSEBIALLELICVARIANTS))
                    .missenseNonBiallelicCount(record.getValue(GENECOPYNUMBER.MISSENSENONBIALLELICVARIANTS))
                    .missenseNonBiallelicPloidy(record.getValue(GENECOPYNUMBER.MISSENSENONBIALLELICPLOIDY))
                    .minMinorAllelePloidy(record.getValue(GENECOPYNUMBER.MINMINORALLELEPLOIDY))
                    .build());
        }
        return geneCopyNumbers;
    }

    void writeCopyNumber(@NotNull final String sample, @NotNull List<GeneCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteGeneCopyNumberForSample(sample);

        for (List<GeneCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(GENECOPYNUMBER,
                    GENECOPYNUMBER.SAMPLEID,
                    GENECOPYNUMBER.CHROMOSOME,
                    GENECOPYNUMBER.START,
                    GENECOPYNUMBER.END,
                    GENECOPYNUMBER.GENE,
                    GENECOPYNUMBER.MINCOPYNUMBER,
                    GENECOPYNUMBER.MAXCOPYNUMBER,
                    GENECOPYNUMBER.SOMATICREGIONS,
                    GENECOPYNUMBER.GERMLINEHOMREGIONS,
                    GENECOPYNUMBER.GERMLINEHETREGIONS,
                    GENECOPYNUMBER.TRANSCRIPTID,
                    GENECOPYNUMBER.TRANSCRIPTVERSION,
                    GENECOPYNUMBER.CHROMOSOMEBAND,
                    GENECOPYNUMBER.MINREGIONS,
                    GENECOPYNUMBER.MINREGIONSTART,
                    GENECOPYNUMBER.MINREGIONEND,
                    GENECOPYNUMBER.MINREGIONSTARTSUPPORT,
                    GENECOPYNUMBER.MINREGIONENDSUPPORT,
                    GENECOPYNUMBER.MINREGIONMETHOD,
                    GENECOPYNUMBER.NONSENSEBIALLELICVARIANTS,
                    GENECOPYNUMBER.NONSENSENONBIALLELICVARIANTS,
                    GENECOPYNUMBER.NONSENSENONBIALLELICPLOIDY,
                    GENECOPYNUMBER.SPLICEBIALLELICVARIANTS,
                    GENECOPYNUMBER.SPLICENONBIALLELICVARIANTS,
                    GENECOPYNUMBER.SPLICENONBIALLELICPLOIDY,
                    GENECOPYNUMBER.MISSENSEBIALLELICVARIANTS,
                    GENECOPYNUMBER.MISSENSENONBIALLELICVARIANTS,
                    GENECOPYNUMBER.MISSENSENONBIALLELICPLOIDY,
                    GENECOPYNUMBER.MINMINORALLELEPLOIDY,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull GeneCopyNumber gene) {
        //noinspection unchecked
        inserter.values(sample,
                gene.chromosome(),
                gene.start(),
                gene.end(),
                gene.gene(),
                DatabaseUtil.decimal(gene.minCopyNumber()),
                DatabaseUtil.decimal(gene.maxCopyNumber()),
                gene.somaticRegions(),
                gene.germlineHomRegions(),
                gene.germlineHet2HomRegions(),
                gene.transcriptID(),
                gene.transcriptVersion(),
                gene.chromosomeBand(),
                gene.minRegions(),
                gene.minRegionStart(),
                gene.minRegionEnd(),
                gene.minRegionStartSupport(),
                gene.minRegionEndSupport(),
                gene.minRegionMethod(),
                gene.nonsenseBiallelicCount(),
                gene.nonsenseNonBiallelicCount(),
                DatabaseUtil.decimal(gene.nonsenseNonBiallelicPloidy()),
                gene.spliceBiallelicCount(),
                gene.spliceNonBiallelicCount(),
                DatabaseUtil.decimal(gene.spliceNonBiallelicPloidy()),
                gene.missenseBiallelicCount(),
                gene.missenseNonBiallelicCount(),
                DatabaseUtil.decimal(gene.missenseNonBiallelicPloidy()),
                DatabaseUtil.decimal(gene.minMinorAllelePloidy()),
                timestamp);
    }

    void deleteGeneCopyNumberForSample(@NotNull String sample) {
        context.delete(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).execute();
    }
}
