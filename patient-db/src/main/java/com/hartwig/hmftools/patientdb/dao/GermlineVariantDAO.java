package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT2;

import java.sql.Timestamp;
import java.util.List;

import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;

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
        context.delete(GERMLINEVARIANT2).where(GERMLINEVARIANT2.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    private InsertValuesStepN createInserter() {
        return context.insertInto(GERMLINEVARIANT2,
                GERMLINEVARIANT2.MODIFIED,
                GERMLINEVARIANT2.SAMPLEID,
                GERMLINEVARIANT2.CHROMOSOME,
                GERMLINEVARIANT2.POSITION,
                GERMLINEVARIANT2.FILTER,
                GERMLINEVARIANT2.TYPE,
                GERMLINEVARIANT2.REF,
                GERMLINEVARIANT2.ALT,
                GERMLINEVARIANT2.QUAL,
                GERMLINEVARIANT2.TIER,
                GERMLINEVARIANT2.GERMLINEGENOTYPE,
                GERMLINEVARIANT2.GERMLINEALLELEREADCOUNT,
                GERMLINEVARIANT2.GERMLINETOTALREADCOUNT,
                GERMLINEVARIANT2.RNAALLELEREADCOUNT,
                GERMLINEVARIANT2.RNATOTALREADCOUNT,
                GERMLINEVARIANT2.TUMORALLELEREADCOUNT,
                GERMLINEVARIANT2.TUMORTOTALREADCOUNT,
                GERMLINEVARIANT2.LOCALPHASESET,
                GERMLINEVARIANT2.ADJUSTEDVAF,
                GERMLINEVARIANT2.VARIANTCOPYNUMBER,
                GERMLINEVARIANT2.COPYNUMBER,
                GERMLINEVARIANT2.BIALLELIC,
                GERMLINEVARIANT2.MINORALLELECOPYNUMBER,
                GERMLINEVARIANT2.CLINVARINFO,
                GERMLINEVARIANT2.PATHOGENICITY,
                GERMLINEVARIANT2.PATHOGENIC,
                GERMLINEVARIANT2.GENE,
                GERMLINEVARIANT2.GENESEFFECTED,
                GERMLINEVARIANT2.WORSTEFFECT,
                GERMLINEVARIANT2.WORSTCODINGEFFECT,
                GERMLINEVARIANT2.WORSTEFFECTTRANSCRIPT,
                GERMLINEVARIANT2.CANONICALEFFECT,
                GERMLINEVARIANT2.CANONICALCODINGEFFECT,
                GERMLINEVARIANT2.CANONICALHGVSCODINGIMPACT,
                GERMLINEVARIANT2.CANONICALHGVSPROTEINIMPACT,
                GERMLINEVARIANT2.MICROHOMOLOGY,
                GERMLINEVARIANT2.REPEATSEQUENCE,
                GERMLINEVARIANT2.REPEATCOUNT,
                GERMLINEVARIANT2.TRINUCLEOTIDECONTEXT,
                GERMLINEVARIANT2.HOTSPOT,
                GERMLINEVARIANT2.MAPPABILITY,
                GERMLINEVARIANT2.REPORTED);
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStepN inserter, String tumorSample, String referenceSample,
            String rnaSample, VariantContext variantContext) {
        final VariantContextDecorator decorator = new VariantContextDecorator(variantContext);
        final AllelicDepth tumorDepth = decorator.allelicDepth(tumorSample);
        final AllelicDepth referenceDepth = decorator.allelicDepth(referenceSample);
        final AllelicDepth rnaDepth = decorator.allelicDepth(rnaSample);
        final SnpEffSummary snpEffSummary = decorator.snpEffSummary();
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
                snpEffSummary.gene(),
                snpEffSummary.genesAffected(),
                snpEffSummary.worstEffect(),
                snpEffSummary.worstCodingEffect() != CodingEffect.UNDEFINED ? snpEffSummary.worstCodingEffect() : Strings.EMPTY,
                snpEffSummary.worstTranscript(),
                snpEffSummary.canonicalEffect(),
                snpEffSummary.canonicalCodingEffect() != CodingEffect.UNDEFINED ? snpEffSummary.canonicalCodingEffect() : Strings.EMPTY,
                snpEffSummary.canonicalHgvsCodingImpact(),
                snpEffSummary.canonicalHgvsProteinImpact(),
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
