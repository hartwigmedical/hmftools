package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Germlinevariant.GERMLINEVARIANT;

import java.sql.Timestamp;
import java.util.List;

import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariantDAO {

    @NotNull
    private final DSLContext context;

    public GermlineVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public BufferedWriter<VariantContext> writer(String tumorSample, String referenceSample, String rnaSample) {
        BufferedWriterConsumer<VariantContext> consumer = new BufferedWriterConsumer<VariantContext>() {
            @Override
            public void initialise() {
                deleteGermlineVariantsForSample(tumorSample);
            }

            @Override
            public void accept(final Timestamp timestamp, final List<VariantContext> entries) {
                writeAll(timestamp, tumorSample, referenceSample, rnaSample, entries);
            }
        };

        return new BufferedWriter<>(consumer);
    }

    private void writeAll(@NotNull final Timestamp timestamp, String tumorSample, String referenceSample, String rnaSample,
            @NotNull List<VariantContext> variants) {
        final InsertValuesStepN inserter = createInserter();
        variants.forEach(variant -> addRecord(timestamp, inserter, tumorSample, referenceSample, rnaSample, variant));
        inserter.execute();
    }

    void deleteGermlineVariantsForSample(@NotNull String sampleId) {
        context.delete(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    private InsertValuesStepN createInserter() {
        return context.insertInto(GERMLINEVARIANT,
                GERMLINEVARIANT.MODIFIED,
                GERMLINEVARIANT.SAMPLEID,
                GERMLINEVARIANT.CHROMOSOME,
                GERMLINEVARIANT.POSITION,
                GERMLINEVARIANT.FILTER,
                GERMLINEVARIANT.TYPE,
                GERMLINEVARIANT.REF,
                GERMLINEVARIANT.ALT,
                GERMLINEVARIANT.QUAL,
                GERMLINEVARIANT.TIER,
                GERMLINEVARIANT.GERMLINEGENOTYPE,
                GERMLINEVARIANT.GERMLINEALLELEREADCOUNT,
                GERMLINEVARIANT.GERMLINETOTALREADCOUNT,
                GERMLINEVARIANT.RNAALLELEREADCOUNT,
                GERMLINEVARIANT.RNATOTALREADCOUNT,
                GERMLINEVARIANT.TUMORALLELEREADCOUNT,
                GERMLINEVARIANT.TUMORTOTALREADCOUNT,
                GERMLINEVARIANT.LOCALPHASESET,
                GERMLINEVARIANT.ADJUSTEDVAF,
                GERMLINEVARIANT.VARIANTCOPYNUMBER,
                GERMLINEVARIANT.COPYNUMBER,
                GERMLINEVARIANT.BIALLELIC,
                GERMLINEVARIANT.MINORALLELECOPYNUMBER,
                GERMLINEVARIANT.CLINVARINFO,
                GERMLINEVARIANT.PATHOGENICITY,
                GERMLINEVARIANT.PATHOGENIC,
                GERMLINEVARIANT.GENE,
                GERMLINEVARIANT.GENESEFFECTED,
                GERMLINEVARIANT.WORSTEFFECT,
                GERMLINEVARIANT.WORSTCODINGEFFECT,
                GERMLINEVARIANT.WORSTEFFECTTRANSCRIPT,
                GERMLINEVARIANT.CANONICALEFFECT,
                GERMLINEVARIANT.CANONICALCODINGEFFECT,
                GERMLINEVARIANT.CANONICALHGVSCODINGIMPACT,
                GERMLINEVARIANT.CANONICALHGVSPROTEINIMPACT,
                GERMLINEVARIANT.MICROHOMOLOGY,
                GERMLINEVARIANT.REPEATSEQUENCE,
                GERMLINEVARIANT.REPEATCOUNT,
                GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                GERMLINEVARIANT.HOTSPOT,
                GERMLINEVARIANT.MAPPABILITY,
                GERMLINEVARIANT.REPORTED);
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStepN inserter, String tumorSample, String referenceSample,
            String rnaSample, VariantContext variantContext) {
        final VariantContextDecorator decorator = new VariantContextDecorator(variantContext);
        final AllelicDepth tumorDepth = decorator.allelicDepth(tumorSample);
        final AllelicDepth referenceDepth = decorator.allelicDepth(referenceSample);
        final AllelicDepth rnaDepth = decorator.allelicDepth(rnaSample);
        final VariantImpact variantImpact = decorator.variantImpact();
        final PathogenicSummary pathogenicSummary = decorator.pathogenicSummary();

        inserter.values(timestamp,
                tumorSample,
                decorator.chromosome(),
                decorator.position(),
                decorator.filter(),
                decorator.type(),
                decorator.ref(),
                decorator.alt(),
                decorator.qual(),
                decorator.tier(),
                decorator.genotypeStatus(referenceSample).simplifiedDisplay(),
                referenceDepth.alleleReadCount(),
                referenceDepth.totalReadCount(),
                rnaDepth.alleleReadCount(),
                rnaDepth.totalReadCount(),
                tumorDepth.alleleReadCount(),
                tumorDepth.totalReadCount(),
                decorator.localPhaseSet(),
                decorator.adjustedVaf(),
                decorator.variantCopyNumber(),
                decorator.adjustedCopyNumber(),
                decorator.biallelic(),
                decorator.minorAlleleCopyNumber(),
                pathogenicSummary.clinvarInfo(),
                pathogenicSummary.pathogenicity().toString(),
                decorator.isPathogenic(),
                variantImpact.gene(),
                variantImpact.GenesAffected,
                variantImpact.WorstEffect,
                variantImpact.WorstCodingEffect != CodingEffect.UNDEFINED ? variantImpact.WorstCodingEffect : Strings.EMPTY,
                variantImpact.WorstTranscript,
                variantImpact.CanonicalEffect,
                variantImpact.CanonicalCodingEffect != CodingEffect.UNDEFINED ? variantImpact.CanonicalCodingEffect : Strings.EMPTY,
                variantImpact.CanonicalHgvsCodingImpact,
                variantImpact.CanonicalHgvsProteinImpact,
                decorator.microhomology(),
                decorator.repeatSequence(),
                decorator.repeatCount(),
                decorator.trinucleotideContext(),
                decorator.hotspot().toString(),
                decorator.mappability(),
                decorator.reported()
        );
    }
}
