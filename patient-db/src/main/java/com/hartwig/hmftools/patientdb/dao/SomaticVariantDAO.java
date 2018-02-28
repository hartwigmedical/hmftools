package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
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

        for (List<EnrichedSomaticVariant> splitRegions : Iterables.partition(variants, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(SOMATICVARIANT,
                    SOMATICVARIANT.SAMPLEID,
                    SOMATICVARIANT.CHROMOSOME,
                    SOMATICVARIANT.POSITION,
                    SOMATICVARIANT.FILTER,
                    SOMATICVARIANT.TYPE,
                    SOMATICVARIANT.REF,
                    SOMATICVARIANT.ALT,
                    SOMATICVARIANT.GENE,
                    SOMATICVARIANT.GENESEFFECTED,
                    SOMATICVARIANT.COSMICID,
                    SOMATICVARIANT.DBSNPID,
                    SOMATICVARIANT.WORSTEFFECT,
                    SOMATICVARIANT.WORSTCODINGEFFECT,
                    SOMATICVARIANT.WORSTEFFECTTRANSCRIPT,
                    SOMATICVARIANT.CANONICALEFFECT,
                    SOMATICVARIANT.CANONICALCODINGEFFECT,
                    SOMATICVARIANT.ALLELEREADCOUNT,
                    SOMATICVARIANT.TOTALREADCOUNT,
                    SOMATICVARIANT.ADJUSTEDCOPYNUMBER,
                    SOMATICVARIANT.ADJUSTEDVAF,
                    SOMATICVARIANT.HIGHCONFIDENCE,
                    SOMATICVARIANT.TRINUCLEOTIDECONTEXT,
                    SOMATICVARIANT.MICROHOMOLOGY,
                    SOMATICVARIANT.REPEATSEQUENCE,
                    SOMATICVARIANT.REPEATCOUNT,
                    SOMATICVARIANT.CLONALITY,
                    SOMATICVARIANT.BIALLELIC,
                    SOMATICVARIANT.HOTSPOT,
                    SOMATICVARIANT.MAPPABILITY,
                    SOMATICVARIANT.GERMLINESTATUS,
                    SOMATICVARIANT.MINORALLELEPLOIDY,
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
                region.genesEffected(),
                region.cosmicID() == null ? "" : region.cosmicID(),
                region.dbsnpID() == null ? "" : region.dbsnpID(),
                region.worstEffect(),
                region.worstCodingEffect(),
                region.worstEffectTranscript(),
                region.canonicalEffect(),
                region.canonicalCodingEffect(),
                region.alleleReadCount(),
                region.totalReadCount(),
                DatabaseUtil.decimal(region.adjustedCopyNumber()),
                DatabaseUtil.decimal(region.adjustedVAF()),
                region.highConfidenceRegion(),
                region.trinucleotideContext(),
                region.microhomology(),
                region.repeatSequence(),
                region.repeatCount(),
                region.clonality(),
                region.biallelic(),
                region.hotspot(),
                DatabaseUtil.decimal(region.mappability()),
                region.germlineStatus(),
                DatabaseUtil.decimal(region.minorAllelePloidy()),
                timestamp);
    }

    @NotNull
    private static String filter(@NotNull String filter) {
        return filter.equals(".") ? PASS : filter;
    }
}