package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

class SomaticVariantDAO {

    @NotNull
    private final DSLContext context;

    SomaticVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull final String sample, @NotNull List<EnrichedSomaticVariant> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteSomaticVariantForSample(sample);

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
            @NotNull EnrichedSomaticVariant variant) {
        inserter.values(sample,
                variant.chromosome(),
                variant.position(),
                variant.filter(),
                variant.type(),
                variant.ref(),
                variant.alt(),
                variant.gene(),
                variant.genesEffected(),
                variant.cosmicID() == null ? "" : variant.cosmicID(),
                variant.dbsnpID() == null ? "" : variant.dbsnpID(),
                variant.worstEffect(),
                variant.worstCodingEffect() == CodingEffect.UNDEFINED ? "" : variant.worstCodingEffect(),
                variant.worstEffectTranscript(),
                variant.canonicalEffect(),
                variant.canonicalCodingEffect() == CodingEffect.UNDEFINED ? "" : variant.canonicalCodingEffect(),
                variant.alleleReadCount(),
                variant.totalReadCount(),
                DatabaseUtil.decimal(variant.adjustedCopyNumber()),
                DatabaseUtil.decimal(variant.adjustedVAF()),
                variant.highConfidenceRegion(),
                variant.trinucleotideContext(),
                variant.microhomology(),
                variant.repeatSequence(),
                variant.repeatCount(),
                variant.clonality(),
                variant.biallelic(),
                variant.hotspot(),
                DatabaseUtil.decimal(variant.mappability()),
                variant.germlineStatus(),
                DatabaseUtil.decimal(variant.minorAllelePloidy()),
                timestamp);
    }

    void deleteSomaticVariantForSample(@NotNull String sample) {
        context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).execute();
    }
}