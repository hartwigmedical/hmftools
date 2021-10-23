package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.genotype.GenotypeStatus.UNKNOWN;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import static org.jooq.impl.DSL.count;

import java.sql.Timestamp;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoad;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsVariant;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsMutationalLoad;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsVariant;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Record3;
import org.jooq.Result;

public class SomaticVariantDAO
{
    @NotNull
    private final DSLContext context;

    SomaticVariantDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    @NotNull
    public BufferedWriter<SomaticVariant> writer(String tumorSample)
    {
        BufferedWriterConsumer<SomaticVariant> consumer = new BufferedWriterConsumer<SomaticVariant>()
        {
            @Override
            public void initialise()
            {
                context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(tumorSample)).execute();
            }

            @Override
            public void accept(final Timestamp timestamp, final List<SomaticVariant> entries)
            {
                writeAll(timestamp, tumorSample, entries);
            }
        };

        return new BufferedWriter<>(consumer);
    }

    public DndsMutationalLoad readDndsLoad(@NotNull String sample)
    {
        Result<Record3<Byte, String, Integer>> result = context.select(SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.TYPE, count())
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.SAMPLEID.eq(sample))
                .and(SOMATICVARIANT.TYPE.in(VariantType.INDEL.toString(), VariantType.SNP.toString()))
                .and(SOMATICVARIANT.FILTER.eq("PASS"))
                .groupBy(SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.TYPE)
                .fetch();

        int snvBiallelic = 0;
        int snvNonBiallelic = 0;
        int indelBiallelic = 0;
        int indelNonBiallelic = 0;

        for(Record3<Byte, String, Integer> record : result)
        {
            boolean isBiallelic = byteToBoolean(record.value1());
            VariantType type = VariantType.valueOf(record.value2());
            int count = record.value3();
            if(isBiallelic && type == VariantType.INDEL)
            {
                indelBiallelic = count;
            }
            if(isBiallelic && type == VariantType.SNP)
            {
                snvBiallelic = count;
            }
            if(!isBiallelic && type == VariantType.INDEL)
            {
                indelNonBiallelic = count;
            }
            if(!isBiallelic && type == VariantType.SNP)
            {
                snvNonBiallelic = count;
            }
        }

        return ImmutableDndsMutationalLoad.builder()
                .sampleId(sample)
                .indelBiallelic(indelBiallelic)
                .indelNonBiallelic(indelNonBiallelic)
                .snvBiallelic(snvBiallelic)
                .snvNonBiallelic(snvNonBiallelic)
                .build();
    }

    @NotNull
    public List<DndsVariant> readDndsVariants(int maxRepeatCount, @NotNull String sample)
    {
        List<DndsVariant> variants = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.SAMPLEID.eq(sample))
                .and(SOMATICVARIANT.FILTER.eq("PASS"))
                .and(SOMATICVARIANT.GENE.ne(""))
                .and(SOMATICVARIANT.REPEATCOUNT.lessOrEqual(maxRepeatCount))
                .and(SOMATICVARIANT.TYPE.in(VariantType.INDEL.toString(), VariantType.SNP.toString()))
                .fetch();

        for(Record record : result)
        {
            variants.add(ImmutableDndsVariant.builder()
                    .sampleId(record.getValue(SOMATICVARIANT.SAMPLEID))
                    .chromosome(record.getValue(SOMATICVARIANT.CHROMOSOME))
                    .position(record.getValue(SOMATICVARIANT.POSITION))
                    .ref(record.getValue(SOMATICVARIANT.REF))
                    .alt(record.getValue(SOMATICVARIANT.ALT))
                    .gene(record.getValue(SOMATICVARIANT.GENE))
                    .worstCodingEffect(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT).isEmpty()
                            ? CodingEffect.UNDEFINED
                            : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT)))
                    .canonicalCodingEffect(record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT).isEmpty()
                            ? CodingEffect.UNDEFINED
                            : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT)))
                    .biallelic(byteToBoolean(record.getValue(SOMATICVARIANT.BIALLELIC)))
                    .repeatCount(record.getValue(SOMATICVARIANT.REPEATCOUNT))
                    .hotspot(Hotspot.valueOf(record.getValue(SOMATICVARIANT.HOTSPOT)).equals(Hotspot.HOTSPOT))
                    .build());
        }
        return variants;
    }

    @NotNull
    public List<SomaticVariant> read(@NotNull String sample, VariantType type)
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        Result<Record> result = type == VariantType.UNDEFINED
                ? context.select().from(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).fetch()
                : context.select()
                        .from(SOMATICVARIANT)
                        .where(SOMATICVARIANT.SAMPLEID.eq(sample))
                        .and(SOMATICVARIANT.TYPE.eq(type.toString()))
                        .fetch();

        for(Record record : result)
        {
            variants.add(buildFromRecord(record));
        }

        return variants;
    }

    public static SomaticVariant buildFromRecord(final Record record)
    {
        Integer referenceAlleleReadCount = record.getValue(SOMATICVARIANT.REFERENCEALLELEREADCOUNT);
        Integer referenceTotalCount = record.getValue(SOMATICVARIANT.REFERENCETOTALREADCOUNT);
        AllelicDepth referenceAllelicDepth = referenceAlleleReadCount != null && referenceTotalCount != null ? ImmutableAllelicDepthImpl
                .builder()
                .alleleReadCount(referenceAlleleReadCount)
                .totalReadCount(referenceTotalCount)
                .build() : null;

        Integer rnaAlleleReadCount = record.getValue(SOMATICVARIANT.RNAALLELEREADCOUNT);
        Integer rnaTotalCount = record.getValue(SOMATICVARIANT.RNATOTALREADCOUNT);
        AllelicDepth rnaAllelicDepth = rnaAlleleReadCount != null && rnaTotalCount != null ? ImmutableAllelicDepthImpl.builder()
                .alleleReadCount(rnaAlleleReadCount)
                .totalReadCount(rnaTotalCount)
                .build() : null;

        return ImmutableSomaticVariantImpl.builder()
                .chromosome(record.getValue(SOMATICVARIANT.CHROMOSOME))
                .position(record.getValue(SOMATICVARIANT.POSITION))
                .filter(record.getValue(SOMATICVARIANT.FILTER))
                .type(VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE)))
                .ref(record.getValue(SOMATICVARIANT.REF))
                .alt(record.getValue(SOMATICVARIANT.ALT))
                .gene(record.getValue(SOMATICVARIANT.GENE))
                .genesAffected(record.getValue(SOMATICVARIANT.GENESEFFECTED))
                .worstCodingEffect(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT).isEmpty()
                        ? CodingEffect.UNDEFINED
                        : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT)))
                .canonicalTranscript("")
                .canonicalEffect(record.getValue(SOMATICVARIANT.CANONICALEFFECT))
                .canonicalCodingEffect(record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT).isEmpty()
                        ? CodingEffect.UNDEFINED
                        : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT)))
                .canonicalHgvsCodingImpact(record.getValue(SOMATICVARIANT.CANONICALHGVSCODINGIMPACT))
                .canonicalHgvsProteinImpact(record.getValue(SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT))
                .spliceRegion(byteToBoolean(record.getValue(SOMATICVARIANT.SPLICEREGION)))
                .otherReportedEffects(DatabaseUtil.valueNotNull(record.getValue(SOMATICVARIANT.OTHERTRANSCRIPTEFFECTS)))
                .alleleReadCount(record.getValue(SOMATICVARIANT.ALLELEREADCOUNT))
                .totalReadCount(record.getValue(SOMATICVARIANT.TOTALREADCOUNT))
                .adjustedCopyNumber(record.getValue(SOMATICVARIANT.COPYNUMBER))
                .adjustedVAF(record.getValue(SOMATICVARIANT.ADJUSTEDVAF))
                .variantCopyNumber(record.getValue(SOMATICVARIANT.VARIANTCOPYNUMBER))
                .biallelic(byteToBoolean(record.getValue(SOMATICVARIANT.BIALLELIC)))
                .reported(byteToBoolean(record.getValue(SOMATICVARIANT.REPORTED)))
                .trinucleotideContext(record.getValue(SOMATICVARIANT.TRINUCLEOTIDECONTEXT))
                .microhomology(record.getValue(SOMATICVARIANT.MICROHOMOLOGY))
                .repeatSequence(record.getValue(SOMATICVARIANT.REPEATSEQUENCE))
                .repeatCount(record.getValue(SOMATICVARIANT.REPEATCOUNT))
                .subclonalLikelihood(record.getValue(SOMATICVARIANT.SUBCLONALLIKELIHOOD))
                .hotspot(Hotspot.valueOf(record.getValue(SOMATICVARIANT.HOTSPOT)))
                .mappability(record.getValue(SOMATICVARIANT.MAPPABILITY))
                .germlineStatus(GermlineStatus.valueOf(record.getValue(SOMATICVARIANT.GERMLINESTATUS)))
                .minorAlleleCopyNumber(record.getValue(SOMATICVARIANT.MINORALLELECOPYNUMBER))
                .recovered(byteToBoolean(record.getValue(SOMATICVARIANT.RECOVERED)))
                .kataegis(record.get(SOMATICVARIANT.KATAEGIS))
                .tier(VariantTier.fromString(record.get(SOMATICVARIANT.TIER)))
                .referenceDepth(referenceAllelicDepth)
                .rnaDepth(rnaAllelicDepth)
                .qual(record.get(SOMATICVARIANT.QUAL))
                .localPhaseSet(record.get(SOMATICVARIANT.LOCALPHASESET))
                .localRealignmentSet(record.get(SOMATICVARIANT.LOCALREALIGNMENTSET))
                .phasedInframeIndelIdentifier(record.get(SOMATICVARIANT.PHASEDINFRAMEINDEL))
                .genotypeStatus(UNKNOWN)
                .build();
    }

    private static boolean byteToBoolean(@Nullable Byte b)
    {
        if(b == null)
        {
            throw new IllegalStateException("NULL value present in database for non-null field");
        }
        return b != 0;
    }

    void writeAll(@NotNull final Timestamp timestamp, @NotNull String sample, @NotNull List<SomaticVariant> variants)
    {
        final InsertValuesStepN inserter = context.insertInto(SOMATICVARIANT,
                SOMATICVARIANT.SAMPLEID,
                SOMATICVARIANT.CHROMOSOME,
                SOMATICVARIANT.POSITION,
                SOMATICVARIANT.FILTER,
                SOMATICVARIANT.TYPE,
                SOMATICVARIANT.REF,
                SOMATICVARIANT.ALT,
                SOMATICVARIANT.GENE,
                SOMATICVARIANT.GENESEFFECTED,
                SOMATICVARIANT.REPORTED,
                SOMATICVARIANT.WORSTCODINGEFFECT,
                SOMATICVARIANT.CANONICALEFFECT,
                SOMATICVARIANT.CANONICALCODINGEFFECT,
                SOMATICVARIANT.CANONICALHGVSCODINGIMPACT,
                SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT,
                SOMATICVARIANT.SPLICEREGION,
                SOMATICVARIANT.OTHERTRANSCRIPTEFFECTS,
                SOMATICVARIANT.ALLELEREADCOUNT,
                SOMATICVARIANT.TOTALREADCOUNT,
                SOMATICVARIANT.COPYNUMBER,
                SOMATICVARIANT.ADJUSTEDVAF,
                SOMATICVARIANT.VARIANTCOPYNUMBER,
                SOMATICVARIANT.TRINUCLEOTIDECONTEXT,
                SOMATICVARIANT.MICROHOMOLOGY,
                SOMATICVARIANT.REPEATSEQUENCE,
                SOMATICVARIANT.REPEATCOUNT,
                SOMATICVARIANT.SUBCLONALLIKELIHOOD,
                SOMATICVARIANT.BIALLELIC,
                SOMATICVARIANT.HOTSPOT,
                SOMATICVARIANT.MAPPABILITY,
                SOMATICVARIANT.GERMLINESTATUS,
                SOMATICVARIANT.MINORALLELECOPYNUMBER,
                SOMATICVARIANT.RECOVERED,
                SOMATICVARIANT.KATAEGIS,
                SOMATICVARIANT.TIER,
                SOMATICVARIANT.REFERENCEALLELEREADCOUNT,
                SOMATICVARIANT.REFERENCETOTALREADCOUNT,
                SOMATICVARIANT.RNAALLELEREADCOUNT,
                SOMATICVARIANT.RNATOTALREADCOUNT,
                SOMATICVARIANT.QUAL,
                SOMATICVARIANT.LOCALPHASESET,
                SOMATICVARIANT.LOCALREALIGNMENTSET,
                SOMATICVARIANT.PHASEDINFRAMEINDEL,
                SOMATICVARIANT.MODIFIED);
        variants.forEach(variant -> addRecord(timestamp, inserter, sample, variant));
        inserter.execute();
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull SomaticVariant variant)
    {
        inserter.values(sample,
                variant.chromosome(),
                variant.position(),
                variant.filter(),
                variant.type(),
                variant.ref(),
                variant.alt(),
                variant.gene(),
                variant.genesAffected(),
                variant.reported(),
                variant.worstCodingEffect() != CodingEffect.UNDEFINED ? variant.worstCodingEffect() : Strings.EMPTY,
                variant.canonicalEffect(),
                variant.canonicalCodingEffect() != CodingEffect.UNDEFINED ? variant.canonicalCodingEffect() : Strings.EMPTY,
                variant.canonicalHgvsCodingImpact(),
                variant.canonicalHgvsProteinImpact(),
                variant.spliceRegion(),
                variant.otherReportedEffects(),
                variant.alleleReadCount(),
                variant.totalReadCount(),
                DatabaseUtil.decimal(variant.adjustedCopyNumber()),
                DatabaseUtil.decimal(variant.adjustedVAF()),
                DatabaseUtil.decimal(variant.variantCopyNumber()),
                variant.trinucleotideContext(),
                variant.microhomology(),
                variant.repeatSequence(),
                variant.repeatCount(),
                variant.subclonalLikelihood(),
                variant.biallelic(),
                variant.hotspot(),
                DatabaseUtil.decimal(variant.mappability()),
                variant.germlineStatus(),
                DatabaseUtil.decimal(variant.minorAlleleCopyNumber()),
                variant.recovered(),
                variant.kataegis(),
                variant.tier().toString(),
                Optional.ofNullable(variant.referenceDepth()).map(AllelicDepth::alleleReadCount).orElse(null),
                Optional.ofNullable(variant.referenceDepth()).map(AllelicDepth::totalReadCount).orElse(null),
                Optional.ofNullable(variant.rnaDepth()).map(AllelicDepth::alleleReadCount).orElse(null),
                Optional.ofNullable(variant.rnaDepth()).map(AllelicDepth::totalReadCount).orElse(null),
                variant.qual(),
                variant.localPhaseSet(),
                variant.localRealignmentSet(),
                variant.phasedInframeIndelIdentifier(),
                timestamp);
    }

    void deleteSomaticVariantForSample(@NotNull String sample)
    {
        context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).execute();
    }

    @NotNull
    List<String> getSamplesList()
    {
        Result<Record1<String>> result =
                context.select(SOMATICVARIANT.SAMPLEID).from(SOMATICVARIANT).groupBy(SOMATICVARIANT.SAMPLEID).fetch();

        List<String> samplesList = Lists.newArrayList();

        for(Record record : result)
        {
            samplesList.add(record.getValue(SOMATICVARIANT.SAMPLEID));
        }

        return samplesList;
    }
}