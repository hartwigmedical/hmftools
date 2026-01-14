package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.checkStringLength;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEDELETION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jooq.DSLContext;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep19;
import org.jooq.Record;
import org.jooq.Record15;
import org.jooq.Result;

public class GeneCopyNumberDAO
{
    private final DSLContext context;

    GeneCopyNumberDAO(final DSLContext context)
    {
        this.context = context;
    }

    public List<GeneCopyNumber> readCopyNumbers(final String sample, final List<String> genes)
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
            GeneCopyNumber geneCopyNumber = new GeneCopyNumber(
                    String.valueOf(record.getValue(GENECOPYNUMBER.CHROMOSOME)),
                    record.getValue(GENECOPYNUMBER.START),
                    record.getValue(GENECOPYNUMBER.END),
                    record.getValue(GENECOPYNUMBER.GENE),
                    record.getValue(GENECOPYNUMBER.TRANSCRIPTID),
                    record.getValue(GENECOPYNUMBER.CANONICALTRANSCRIPT) == 1,
                    record.getValue(GENECOPYNUMBER.CHROMOSOMEBAND),
                    record.getValue(GENECOPYNUMBER.MAXCOPYNUMBER),
                    record.getValue(GENECOPYNUMBER.MINCOPYNUMBER),
                    record.getValue(GENECOPYNUMBER.MINMINORALLELECOPYNUMBER),
                    record.getValue(GENECOPYNUMBER.SOMATICREGIONS),
                    record.getValue(GENECOPYNUMBER.MINREGIONS),
                    record.getValue(GENECOPYNUMBER.MINREGIONSTART),
                    record.getValue(GENECOPYNUMBER.MINREGIONEND),
                    0, 0,
                    SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONSTARTSUPPORT)),
                    SegmentSupport.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONENDSUPPORT)),
                    CopyNumberMethod.valueOf(record.getValue(GENECOPYNUMBER.MINREGIONMETHOD)), 0);

            geneCopyNumbers.add(geneCopyNumber);
        }
        return geneCopyNumbers;
    }

    public List<GermlineAmpDel> readGermlineDeletions(final String sample)
    {
        List<GermlineAmpDel> germlineDeletions = Lists.newArrayList();

        Result<Record15<String,String,Integer,Integer,Integer,Integer,Integer,String,String,String,Double,Double,String,Integer,Byte>> result = context.select(
                GERMLINEDELETION.GENE, GERMLINEDELETION.CHROMOSOME, GERMLINEDELETION.REGIONSTART, GERMLINEDELETION.REGIONEND,
                GERMLINEDELETION.DEPTHWINDOWCOUNT, GERMLINEDELETION.EXONSTART, GERMLINEDELETION.EXONEND,
                GERMLINEDELETION.DETECTIONMETHOD, GERMLINEDELETION.GERMLINESTATUS, GERMLINEDELETION.TUMORSTATUS,
                GERMLINEDELETION.GERMLINECOPYNUMBER, GERMLINEDELETION.TUMORCOPYNUMBER, GERMLINEDELETION.FILTER,
                GERMLINEDELETION.COHORTFREQUENCY, GERMLINEDELETION.REPORTED)
                .from(GERMLINEDELETION).where(GERMLINEDELETION.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            germlineDeletions.add(new GermlineAmpDel(
                    record.getValue(GERMLINEDELETION.GENE),
                    record.getValue(GERMLINEDELETION.CHROMOSOME),
                    "", // record.getValue(GERMLINEDELETION.CHROMOSOMEBAND), // until 5.29 DB changes are applied to prod
                    record.getValue(GERMLINEDELETION.REGIONSTART),
                    record.getValue(GERMLINEDELETION.REGIONEND),
                    record.getValue(GERMLINEDELETION.DEPTHWINDOWCOUNT),
                    record.getValue(GERMLINEDELETION.EXONSTART),
                    record.getValue(GERMLINEDELETION.EXONEND),
                    GermlineDetectionMethod.valueOf(record.getValue(GERMLINEDELETION.DETECTIONMETHOD)),
                    GermlineStatus.valueOf(record.getValue(GERMLINEDELETION.GERMLINESTATUS)),
                    GermlineStatus.valueOf(record.getValue(GERMLINEDELETION.TUMORSTATUS)),
                    record.getValue(GERMLINEDELETION.GERMLINECOPYNUMBER),
                    record.getValue(GERMLINEDELETION.TUMORCOPYNUMBER),
                    record.getValue(GERMLINEDELETION.FILTER),
                    record.getValue(GERMLINEDELETION.COHORTFREQUENCY),
                    record.getValue(GERMLINEDELETION.REPORTED).intValue() == 1 ? ReportedStatus.REPORTED : ReportedStatus.NONE));
        }
        return germlineDeletions;
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
                gene.SomaticRegions,
                gene.TransName,
                gene.IsCanonical,
                gene.ChromosomeBand,
                gene.MinRegions,
                gene.MinRegionStart,
                gene.MinRegionEnd,
                gene.MinRegionStartSupport,
                gene.MinRegionEndSupport,
                gene.MinRegionMethod,
                DatabaseUtil.decimal(gene.MinMinorAlleleCopyNumber),
                timestamp);
    }

    public void deleteGeneCopyNumberForSample(final String sample)
    {
        context.delete(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).execute();
    }

    public void writeGermlineDeletions(final String sample, final List<GermlineAmpDel> deletions)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteGermlineDeletionsForSample(sample);

        InsertValuesStep18 inserter = context.insertInto(GERMLINEDELETION,
                GERMLINEDELETION.SAMPLEID,
                GERMLINEDELETION.GENE,
                GERMLINEDELETION.CHROMOSOME,
                GERMLINEDELETION.CHROMOSOMEBAND,
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

        for(GermlineAmpDel deletion : deletions)
        {
            addDeletionRecord(timestamp, inserter, sample, deletion);
        }

        inserter.execute();
    }

    private static void addDeletionRecord(
            final Timestamp timestamp, final InsertValuesStep18 inserter, final String sample, final GermlineAmpDel deletion)
    {
        inserter.values(
                sample,
                DatabaseUtil.checkStringLength(deletion.GeneName, GERMLINEDELETION.GENE),
                deletion.Chromosome,
                deletion.ChromosomeBand,
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
                checkStringLength(deletion.Filter, GERMLINEDELETION.FILTER),
                deletion.CohortFrequency,
                deletion.Reported == ReportedStatus.REPORTED ? 1 : 0,
                timestamp);
    }

    public void deleteGermlineDeletionsForSample(final String sample)
    {
        context.delete(GERMLINEDELETION).where(GERMLINEDELETION.SAMPLEID.eq(sample)).execute();
    }

}
