package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.ResultSet;
import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

class GeneCopyNumberDAO {

    @NotNull
    private final DSLContext context;

    GeneCopyNumberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeCopyNumber(@NotNull final String sample, @NotNull List<GeneCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(GENECOPYNUMBER).where(GENECOPYNUMBER.SAMPLEID.eq(sample)).execute();

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

    @NotNull
    Stream<PotentialActionableCNV> potentiallyActionableCNVs() {
        return context.select(GENECOPYNUMBER.SAMPLEID, BASELINE.PRIMARYTUMORLOCATION, GENECOPYNUMBER.GENE, GENECOPYNUMBER.MINCOPYNUMBER)
                .from(GENECOPYNUMBER.join(PURITY)
                        .on(GENECOPYNUMBER.SAMPLEID.eq(PURITY.SAMPLEID))
                        .join(SAMPLE)
                        .on(GENECOPYNUMBER.SAMPLEID.eq(SAMPLE.SAMPLEID))
                        .join(BASELINE)
                        .on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(GENECOPYNUMBER.MINCOPYNUMBER.ge(8.0)
                        .or(GENECOPYNUMBER.MINCOPYNUMBER.le(0.5))
                        .or(GENECOPYNUMBER.MINCOPYNUMBER.div(PURITY.PLOIDY).ge(3.0)))
                .fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream()
                .map(PotentialActionableCNV::of)
                .filter(cnv -> cnv.alteration() != CopyNumberAlteration.NEUTRAL);
    }

    private static void addCopynumberRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull GeneCopyNumber gene) {
        //noinspection unchecked
        inserter.values(sample,
                gene.chromosome(),
                gene.start(),
                gene.end(),
                gene.gene(),
                gene.minCopyNumber(),
                gene.maxCopyNumber(),
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
                gene.nonsenseNonBiallelicPloidy(),
                gene.spliceBiallelicCount(),
                gene.spliceNonBiallelicCount(),
                gene.spliceNonBiallelicPloidy(),
                gene.missenseBiallelicCount(),
                gene.missenseNonBiallelicCount(),
                gene.missenseNonBiallelicPloidy(),
                gene.minMinorAllelePloidy(),
                timestamp);
    }
}
