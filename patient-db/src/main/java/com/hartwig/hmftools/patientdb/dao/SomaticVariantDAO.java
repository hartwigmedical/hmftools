package com.hartwig.hmftools.patientdb.dao;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genotype.GenotypeStatus.UNKNOWN;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.byteToBoolean;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.checkStringLength;
import static com.hartwig.hmftools.patientdb.dao.GermlineVariantDAO.checkTrimHgsvString;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Result;

public class SomaticVariantDAO
{
    private final DSLContext context;

    private static final int DB_BATCH_INSERT_SIZE = 10000;

    SomaticVariantDAO(final DSLContext context)
    {
        this.context = context;
    }

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

        return new BufferedWriter<>(consumer, DB_BATCH_INSERT_SIZE);
    }

    public List<SomaticVariant> read(final String sample, VariantType type)
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

        AllelicDepth referenceAllelicDepth = referenceAlleleReadCount != null && referenceTotalCount != null ?
                new AllelicDepth(referenceTotalCount, referenceAlleleReadCount) : null;

        Integer rnaAlleleReadCount = record.getValue(SOMATICVARIANT.RNAALLELEREADCOUNT);
        Integer rnaTotalCount = record.getValue(SOMATICVARIANT.RNATOTALREADCOUNT);
        AllelicDepth rnaAllelicDepth = rnaAlleleReadCount != null && rnaTotalCount != null ?
                new AllelicDepth(rnaTotalCount, rnaAlleleReadCount) : null;

        return ImmutableSomaticVariantImpl.builder()
                .chromosome(record.getValue(SOMATICVARIANT.CHROMOSOME))
                .position(record.getValue(SOMATICVARIANT.POSITION))
                .filter(record.getValue(SOMATICVARIANT.FILTER))
                .type(VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE)))
                .ref(record.getValue(SOMATICVARIANT.REF))
                .alt(record.getValue(SOMATICVARIANT.ALT))
                .gene(record.getValue(SOMATICVARIANT.GENE))
                .genesAffected(record.getValue(SOMATICVARIANT.GENESAFFECTED))
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
                .allelicDepth(new AllelicDepth(record.getValue(SOMATICVARIANT.TOTALREADCOUNT), record.getValue(SOMATICVARIANT.ALLELEREADCOUNT)))
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
                .localPhaseSets(SomaticVariantFactory.localPhaseSetsStringToList(record.get(SOMATICVARIANT.LOCALPHASESET)))
                .genotypeStatus(UNKNOWN)
                .build();
    }


    void writeAll(final Timestamp timestamp, final String sample, final List<SomaticVariant> variants)
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
                SOMATICVARIANT.GENESAFFECTED,
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
                SOMATICVARIANT.CLINVARINFO,
                SOMATICVARIANT.GNOMADFREQUENCY,
                SOMATICVARIANT.SOMATICLIKELIHOOD,
                SOMATICVARIANT.MODIFIED);
        variants.forEach(variant -> addRecord(timestamp, inserter, sample, variant));
        inserter.execute();
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStepN inserter, String sample, SomaticVariant variant)
    {
        // append reportable status for each transcript where non-canonical may be reportable
        String otherReportedEffects = variant.otherReportedEffects();

        if(variant.reportableTranscripts() != null)
        {
            boolean hasCanonical = false;

            for(String reportedTrans : variant.reportableTranscripts())
            {
                if(reportedTrans.equals(variant.canonicalTranscript()))
                    hasCanonical = true;
                else
                    otherReportedEffects = otherReportedEffects.replaceFirst(reportedTrans, format("%s|REPORTED", reportedTrans));
            }

            if(!hasCanonical)
                otherReportedEffects = otherReportedEffects + ";CANONICAL_NOT_REPORTED";
        }

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
                checkTrimHgsvString(variant.canonicalHgvsCodingImpact(), SOMATICVARIANT.CANONICALHGVSCODINGIMPACT),
                checkTrimHgsvString(variant.canonicalHgvsProteinImpact(), SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT),
                variant.spliceRegion(),
                otherReportedEffects,
                variant.allelicDepth().AlleleReadCount,
                variant.allelicDepth().TotalReadCount,
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
                Optional.ofNullable(variant.referenceDepth() != null ? variant.referenceDepth().AlleleReadCount : null),
                Optional.ofNullable(variant.referenceDepth() != null ? variant.referenceDepth().TotalReadCount : null),
                Optional.ofNullable(variant.rnaDepth() != null ? variant.rnaDepth().AlleleReadCount : null),
                Optional.ofNullable(variant.rnaDepth() != null ? variant.rnaDepth().TotalReadCount : null),
                variant.qual(),
                variant.localPhaseSets() != null ? checkStringLength(variant.localPhaseSetsStr(), SOMATICVARIANT.LOCALPHASESET) : null,
                variant.clinvarInfo(),
                variant.gnomadFrequency(),
                variant.somaticLikelihood() == SomaticLikelihood.UNKNOWN ? Strings.EMPTY : variant.somaticLikelihood().toString(),
                timestamp);
    }

    void deleteSomaticVariantForSample(String sample)
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