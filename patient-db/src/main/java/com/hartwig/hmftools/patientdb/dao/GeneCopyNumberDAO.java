package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEDELETION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep17;
import org.jooq.InsertValuesStep19;
import org.jooq.Record;
import org.jooq.Result;

class GeneCopyNumberDAO
{
    @NotNull
    private final DSLContext context;

    GeneCopyNumberDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    @NotNull
    List<GeneCopyNumber> read(@NotNull String sample, @NotNull List<String> genes)
    {
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        Result<Record> result = genes.isEmpty()
                ? context.select().from(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).fetch()
                : context.select()
                        .from(GENECOPYNUMBER)
                        .where(GENECOPYNUMBER.SAMPLEID.eq(sample))
                        .and(GENECOPYNUMBER.GENE.in(genes))
                        .fetch();

        for(Record record : result)
        {
            geneCopyNumbers.add(ImmutableGeneCopyNumber.builder()
                    .chromosome(String.valueOf(record.getValue(GENECOPYNUMBER.CHROMOSOME)))
                    .start(record.getValue(GENECOPYNUMBER.START))
                    .end(record.getValue(GENECOPYNUMBER.END))
                    .geneName(record.getValue(GENECOPYNUMBER.GENE))
                    .minCopyNumber(record.getValue(GENECOPYNUMBER.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(GENECOPYNUMBER.MAXCOPYNUMBER))
                    .somaticRegions(record.getValue(GENECOPYNUMBER.SOMATICREGIONS))
                    .transName(record.getValue(GENECOPYNUMBER.TRANSCRIPTID))
                    .isCanonical(record.getValue(GENECOPYNUMBER.CANONICALTRANSCRIPT) == 1)
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

    public void writeCopyNumber(final String sample, final List<GeneCopyNumber> copyNumbers)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteGeneCopyNumberForSample(sample);

        for(List<GeneCopyNumber> splitCopyNumbers : Iterables.partition(copyNumbers, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep19 inserter = context.insertInto(GENECOPYNUMBER,
                    GENECOPYNUMBER.SAMPLEID,
                    GENECOPYNUMBER.CHROMOSOME,
                    GENECOPYNUMBER.START,
                    GENECOPYNUMBER.END,
                    GENECOPYNUMBER.GENE,
                    GENECOPYNUMBER.MINCOPYNUMBER,
                    GENECOPYNUMBER.MAXCOPYNUMBER,
                    GENECOPYNUMBER.SOMATICREGIONS,
                    GENECOPYNUMBER.TRANSCRIPTID,
                    GENECOPYNUMBER.CANONICALTRANSCRIPT,
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

    private static void addCopynumberRecord(
            final Timestamp timestamp, final InsertValuesStep19 inserter, final String sample, final GeneCopyNumber gene)
    {
        inserter.values(sample,
                gene.chromosome(),
                gene.start(),
                gene.end(),
                gene.geneName(),
                DatabaseUtil.decimal(gene.minCopyNumber()),
                DatabaseUtil.decimal(gene.maxCopyNumber()),
                gene.somaticRegions(),
                gene.transName(),
                gene.isCanonical(),
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

    void deleteGeneCopyNumberForSample(final String sample)
    {
        context.delete(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).execute();
    }

    public void writeGermlineDeletions(final String sample, final List<GermlineDeletion> deletions)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteGermlineDeletionsForSample(sample);

        InsertValuesStep17 inserter = context.insertInto(GERMLINEDELETION,
                GERMLINEDELETION.SAMPLEID,
                GERMLINEDELETION.GENE,
                GERMLINEDELETION.CHROMOSOME,
                GERMLINEDELETION.REGIONSTART,
                GERMLINEDELETION.REGIONEND,
                GERMLINEDELETION.DEPTHWINDOWCOUNT,
                GERMLINEDELETION.EXONSTART,
                GERMLINEDELETION.EXONEND,
                GERMLINEDELETION.DETECTIONMETHOD,
                GERMLINEDELETION.GERMLINESTATUS,
                GERMLINEDELETION.TUMORSTATUS,
                GERMLINEDELETION.GERMLINECOPYNUMBER,
                GERMLINEDELETION.TUMORCOPYNUMBER,
                GERMLINEDELETION.FILTER,
                GERMLINEDELETION.COHORTFREQUENCY,
                GERMLINEDELETION.REPORTED,
                COPYNUMBER.MODIFIED);

        for(GermlineDeletion deletion : deletions)
        {
            addDeletionRecord(timestamp, inserter, sample, deletion);
        }

        inserter.execute();
    }

    private static void addDeletionRecord(
            final Timestamp timestamp, final InsertValuesStep17 inserter, final String sample, final GermlineDeletion deletion)
    {
        inserter.values(
                sample,
                deletion.GeneName,
                deletion.Chromosome,
                deletion.RegionStart,
                deletion.RegionEnd,
                deletion.DepthWindowCount,
                deletion.ExonStart,
                deletion.ExonEnd,
                deletion.DetectionMethod.toString(),
                deletion.NormalStatus.toString(),
                deletion.TumorStatus.toString(),
                DatabaseUtil.decimal(deletion.GermlineCopyNumber),
                DatabaseUtil.decimal(deletion.TumorCopyNumber),
                deletion.Filter,
                deletion.CohortFrequency,
                deletion.Reported,
                timestamp);
    }

    private void deleteGermlineDeletionsForSample(final String sample)
    {
        context.delete(GERMLINEDELETION).where(GERMLINEDELETION.SAMPLEID.eq(sample)).execute();
    }

}
