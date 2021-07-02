package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
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
import org.jooq.InsertValuesStep20;
import org.jooq.Record;
import org.jooq.Result;

class GeneCopyNumberDAO {

    @NotNull
    private final DSLContext context;

    GeneCopyNumberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    List<GeneCopyNumber> read(@NotNull String sample, @NotNull List<String> genes) {
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        Result<Record> result = genes.isEmpty()
                ? context.select().from(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).fetch()
                : context.select()
                        .from(GENECOPYNUMBER)
                        .where(GENECOPYNUMBER.SAMPLEID.eq(sample))
                        .and(GENECOPYNUMBER.GENE.in(genes))
                        .fetch();

        for (Record record : result) {
            geneCopyNumbers.add(ImmutableGeneCopyNumber.builder()
                    .chromosome(String.valueOf(record.getValue(GENECOPYNUMBER.CHROMOSOME)))
                    .start(record.getValue(GENECOPYNUMBER.START))
                    .end(record.getValue(GENECOPYNUMBER.END))
                    .gene(record.getValue(GENECOPYNUMBER.GENE))
                    .minCopyNumber(record.getValue(GENECOPYNUMBER.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(GENECOPYNUMBER.MAXCOPYNUMBER))
                    .somaticRegions(record.getValue(GENECOPYNUMBER.SOMATICREGIONS))
                    .germlineHomRegions(record.getValue(GENECOPYNUMBER.GERMLINEHOMDELETIONREGIONS))
                    .germlineHet2HomRegions(record.getValue(GENECOPYNUMBER.GERMLINEHETTOHOMDELETIONREGIONS))
                    .transcriptID(record.getValue(GENECOPYNUMBER.TRANSCRIPTID))
                    .chromosomeBand(record.getValue(GENECOPYNUMBER.CHROMOSOMEBAND))
                    .minRegions(record.getValue(GENECOPYNUMBER.MINREGIONS))
                    .minRegionStart(record.getValue(GENECOPYNUMBER.MINREGIONSTART))
                    .minRegionEnd(record.getValue(GENECOPYNUMBER.MINREGIONEND))
                    .minRegionStartSupport(SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONSTARTSUPPORT)))
                    .minRegionEndSupport(SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONENDSUPPORT)))
                    .minRegionMethod(CopyNumberMethod.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONMETHOD)))
                    .minMinorAlleleCopyNumber(record.getValue(GENECOPYNUMBER.MINMINORALLELECOPYNUMBER))
                    .build());
        }
        return geneCopyNumbers;
    }

    void writeCopyNumber(@NotNull String sample, @NotNull List<GeneCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteGeneCopyNumberForSample(sample);

        for (List<GeneCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep20 inserter = context.insertInto(GENECOPYNUMBER,
                    GENECOPYNUMBER.SAMPLEID,
                    GENECOPYNUMBER.CHROMOSOME,
                    GENECOPYNUMBER.START,
                    GENECOPYNUMBER.END,
                    GENECOPYNUMBER.GENE,
                    GENECOPYNUMBER.MINCOPYNUMBER,
                    GENECOPYNUMBER.MAXCOPYNUMBER,
                    GENECOPYNUMBER.SOMATICREGIONS,
                    GENECOPYNUMBER.GERMLINEHOMDELETIONREGIONS,
                    GENECOPYNUMBER.GERMLINEHETTOHOMDELETIONREGIONS,
                    GENECOPYNUMBER.TRANSCRIPTID,
                    GENECOPYNUMBER.CHROMOSOMEBAND,
                    GENECOPYNUMBER.MINREGIONS,
                    GENECOPYNUMBER.MINREGIONSTART,
                    GENECOPYNUMBER.MINREGIONEND,
                    GENECOPYNUMBER.MINREGIONSTARTSUPPORT,
                    GENECOPYNUMBER.MINREGIONENDSUPPORT,
                    GENECOPYNUMBER.MINREGIONMETHOD,
                    GENECOPYNUMBER.MINMINORALLELECOPYNUMBER,
                    COPYNUMBER.MODIFIED);
            splitCopyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep20 inserter, @NotNull String sample,
            @NotNull GeneCopyNumber gene) {
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
                gene.chromosomeBand(),
                gene.minRegions(),
                gene.minRegionStart(),
                gene.minRegionEnd(),
                gene.minRegionStartSupport(),
                gene.minRegionEndSupport(),
                gene.minRegionMethod(),
                DatabaseUtil.decimal(gene.minMinorAlleleCopyNumber()),
                timestamp);
    }

    void deleteGeneCopyNumberForSample(@NotNull String sample) {
        context.delete(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).execute();
    }
}
