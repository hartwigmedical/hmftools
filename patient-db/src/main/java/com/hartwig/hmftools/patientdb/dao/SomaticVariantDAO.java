package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

class SomaticVariantDAO {
    private static final String PASS = "PASS";

    @NotNull
    private final DSLContext context;

    SomaticVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull final String sample, @NotNull List<EnrichedSomaticVariant> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<EnrichedSomaticVariant> splitRegions : Iterables.partition(variants, BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(SOMATICVARIANT,
                    SOMATICVARIANT.SAMPLEID,
                    SOMATICVARIANT.CHROMOSOME,
                    SOMATICVARIANT.POSITION,
                    SOMATICVARIANT.FILTER,
                    SOMATICVARIANT.TYPE,
                    SOMATICVARIANT.REF,
                    SOMATICVARIANT.ALT,
                    SOMATICVARIANT.GENE,
                    SOMATICVARIANT.COSMICID,
                    SOMATICVARIANT.DBSNPID,
                    SOMATICVARIANT.EFFECT,
                    SOMATICVARIANT.ALLELEREADCOUNT,
                    SOMATICVARIANT.TOTALREADCOUNT,
                    SOMATICVARIANT.ADJUSTEDCOPYNUMBER,
                    SOMATICVARIANT.ADJUSTEDVAF,
                    SOMATICVARIANT.HIGHCONFIDENCE,
                    SOMATICVARIANT.TRINUCLEOTIDECONTEXT,
                    SOMATICVARIANT.MICROHOMOLOGY,
                    SOMATICVARIANT.REPEATSEQUENCE,
                    SOMATICVARIANT.REPEATCOUNT,
                    SOMATICVARIANT.REFGENOMECONTEXT,
                    SOMATICVARIANT.CLONALITY,
                    SOMATICVARIANT.LOH,
                    SOMATICVARIANT.MAPPABILITY,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull EnrichedSomaticVariant region) {
        inserter.values(sample,
                region.chromosome(),
                region.position(),
                filter(region.filter()),
                region.type(),
                region.ref(),
                region.alt(),
                region.gene(),
                region.cosmicID() == null ? "" : region.cosmicID(),
                region.dbsnpID() == null ? "" : region.dbsnpID(),
                region.effect(),
                region.alleleReadCount(),
                region.totalReadCount(),
                region.adjustedCopyNumber(),
                region.adjustedVAF(),
                region.highConfidenceRegion(),
                region.trinucleotideContext(),
                region.microhomology(),
                region.repeatSequence(),
                region.repeatCount(),
                region.refGenomeContext(),
                region.clonality(),
                region.lossOfHeterozygosity(),
                region.mappability(),
                timestamp);
    }

    @NotNull
    private static String filter(@NotNull String filter) {
        return filter.equals(".") ? PASS : filter;
    }
}