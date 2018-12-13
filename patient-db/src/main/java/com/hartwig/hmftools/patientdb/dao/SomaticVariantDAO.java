package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Result;

class SomaticVariantDAO {

    @NotNull
    private final DSLContext context;

    SomaticVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public final List<EnrichedSomaticVariant> read(@NotNull final String sample) {
        List<EnrichedSomaticVariant> variants = Lists.newArrayList();

        final Result<Record> result = context.select().from(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            variants.add(ImmutableEnrichedSomaticVariant.builder()
                    .chromosome(record.getValue(SOMATICVARIANT.CHROMOSOME))
                    .position(record.getValue(SOMATICVARIANT.POSITION))
                    .filter(record.getValue(SOMATICVARIANT.FILTER))
                    .type(VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE)))
                    .ref(record.getValue(SOMATICVARIANT.REF))
                    .alt(record.getValue(SOMATICVARIANT.ALT))
                    .gene(record.getValue(SOMATICVARIANT.GENE))
                    .genesEffected(record.getValue(SOMATICVARIANT.GENESEFFECTED))
                    .cosmicIDs(Lists.newArrayList())
                    .dbsnpID(record.getValue(SOMATICVARIANT.DBSNPID))
                    .worstEffect(record.getValue(SOMATICVARIANT.WORSTEFFECT))
                    .worstCodingEffect(CodingEffect.NONE)
                    .worstEffectTranscript(record.getValue(SOMATICVARIANT.WORSTEFFECTTRANSCRIPT))
                    .canonicalEffect(record.getValue(SOMATICVARIANT.CANONICALEFFECT))
                    .canonicalHgvsCodingImpact(record.getValue(SOMATICVARIANT.CANONICALHGVSCODINGIMPACT))
                    .canonicalHgvsProteinImpact(record.getValue(SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT))
                    .alleleReadCount(record.getValue(SOMATICVARIANT.ALLELEREADCOUNT))
                    .totalReadCount(record.getValue(SOMATICVARIANT.TOTALREADCOUNT))
                    .adjustedCopyNumber(record.getValue(SOMATICVARIANT.ADJUSTEDCOPYNUMBER))
                    .adjustedVAF(record.getValue(SOMATICVARIANT.ADJUSTEDVAF))
                    .highConfidenceRegion(byteToBoolean(record.getValue(SOMATICVARIANT.HIGHCONFIDENCE)))
                    .trinucleotideContext(record.getValue(SOMATICVARIANT.TRINUCLEOTIDECONTEXT))
                    .microhomology(record.getValue(SOMATICVARIANT.MICROHOMOLOGY))
                    .repeatSequence(record.getValue(SOMATICVARIANT.REPEATSEQUENCE))
                    .repeatCount(record.getValue(SOMATICVARIANT.REPEATCOUNT))
                    .clonality(Clonality.valueOf(record.getValue(SOMATICVARIANT.CLONALITY)))
                    .biallelic(byteToBoolean(record.getValue(SOMATICVARIANT.BIALLELIC)))
                    .hotspot(Hotspot.valueOf(record.getValue(SOMATICVARIANT.HOTSPOT)))
                    .mappability(record.getValue(SOMATICVARIANT.MAPPABILITY))
                    .germlineStatus(GermlineStatus.valueOf(record.getValue(SOMATICVARIANT.GERMLINESTATUS)))
                    .minorAllelePloidy(record.getValue(SOMATICVARIANT.MINORALLELEPLOIDY))
                    .ploidy(1)
                    .canonicalCodingEffect(CodingEffect.NONE)
                    .recovered(byteToBoolean(record.getValue(SOMATICVARIANT.RECOVERED)))
                    .build());
        }
        return variants;
    }

    @Nullable
    private static Boolean byteToBoolean(Byte b) {
        if (b == null) {
            return null;
        }
        return b != 0;
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
                    SOMATICVARIANT.CANONICALHGVSCODINGIMPACT,
                    SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT,
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
                    SOMATICVARIANT.RECOVERED,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(variant -> addRecord(timestamp, inserter, sample, variant));
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
                variant.canonicalCosmicID() != null ? variant.canonicalCosmicID() : Strings.EMPTY,
                variant.dbsnpID() != null ? variant.dbsnpID() : Strings.EMPTY,
                variant.worstEffect(),
                variant.worstCodingEffect() != CodingEffect.UNDEFINED ? variant.worstCodingEffect() : Strings.EMPTY,
                variant.worstEffectTranscript(),
                variant.canonicalEffect(),
                variant.canonicalCodingEffect() != CodingEffect.UNDEFINED ? variant.canonicalCodingEffect() : Strings.EMPTY,
                variant.canonicalHgvsCodingImpact(),
                variant.canonicalHgvsProteinImpact(),
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
                variant.recovered(),
                timestamp);
    }

    void deleteSomaticVariantForSample(@NotNull String sample) {
        context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).execute();
    }
}