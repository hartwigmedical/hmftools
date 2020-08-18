package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.enrich.SnpEffEnrichment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SnpEffSummaryFactory {

    private final CanonicalAnnotation canonicalAnnotationFactory;

    public SnpEffSummaryFactory(@NotNull final DriverGenePanel genePanel, @NotNull final List<CanonicalTranscript> transcripts) {
        this.canonicalAnnotationFactory = new CanonicalAnnotation(genePanel, transcripts);
    }

    @NotNull
    public static SnpEffSummary fromSage(@NotNull final VariantContext context) {
        final List<String> worst = context.getAttributeAsStringList(SnpEffEnrichment.SNPEFF_WORST, Strings.EMPTY);
        final List<String> canonical = context.getAttributeAsStringList(SnpEffEnrichment.SNPEFF_CANONICAL, Strings.EMPTY);
        return SnpEffSummarySerialiser.fromDetails(worst, canonical);
    }

    @NotNull
    public SnpEffSummary fromAnnotations(@NotNull final VariantContext context) {
        final boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        return fromAnnotations(phasedInframeIndel, allAnnotations);
    }

    @NotNull
    public SnpEffSummary fromAnnotations(boolean phasedInframeIndel, @NotNull final List<SnpEffAnnotation> allAnnotations) {
        final ImmutableSnpEffSummary.Builder builder = SnpEffSummarySerialiser.createBuilder();

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        if (!transcriptAnnotations.isEmpty()) {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            builder.worstGene(worstAnnotation.gene());
            builder.worstEffect(worstAnnotation.consequenceString());
            builder.worstCodingEffect(codingEffect(phasedInframeIndel, worstAnnotation.gene(), worstAnnotation.consequences()));
            builder.worstTranscript(worstAnnotation.transcript());
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = canonicalAnnotationFactory.canonicalSnpEffAnnotation(transcriptAnnotations);
        if (canonicalAnnotation.isPresent()) {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            builder.canonicalGene(annotation.gene());
            builder.canonicalEffect(annotation.consequenceString());
            builder.canonicalCodingEffect(codingEffect(phasedInframeIndel, annotation.gene(), annotation.consequences()));
            builder.canonicalHgvsCodingImpact(annotation.hgvsCoding());
            builder.canonicalHgvsProteinImpact(annotation.hgvsProtein());
            builder.canonicalTranscript(annotation.transcript());
        }

        builder.genesAffected((int) transcriptAnnotations.stream()
                .map(SnpEffAnnotation::gene)
                .filter(x -> !x.isEmpty())
                .distinct()
                .count());

        return builder.build();
    }

    private CodingEffect codingEffect(boolean phasedInframeIndel, @NotNull final String gene,
            @NotNull final List<VariantConsequence> consequences) {
        return phasedInframeIndel ? CodingEffect.MISSENSE : CodingEffect.effect(gene, consequences);
    }

}
