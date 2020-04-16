package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SnpEffSummaryFactory {

    private final CanonicalAnnotation canonicalAnnotationFactory;

    public SnpEffSummaryFactory(@NotNull final CanonicalAnnotation canonicalAnnotationFactory) {
        this.canonicalAnnotationFactory = canonicalAnnotationFactory;
    }

    public SnpEffSummaryFactory(@NotNull final List<CanonicalTranscript> transcripts) {
        this.canonicalAnnotationFactory = new CanonicalAnnotation(transcripts);
    }

    @NotNull
    public SnpEffSummary fromAnnotations(@NotNull final List<SnpEffAnnotation> allAnnotations) {
        final ImmutableSnpEffSummary.Builder builder = ImmutableSnpEffSummary.builder()
                .genesAffected(0)
                .worstGene(Strings.EMPTY)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .worstTranscript(Strings.EMPTY)
                .canonicalGene(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY);

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        if (!transcriptAnnotations.isEmpty()) {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            builder.worstGene(worstAnnotation.gene());
            builder.worstEffect(worstAnnotation.consequenceString());
            builder.worstCodingEffect(CodingEffect.effect(worstAnnotation.gene(), worstAnnotation.consequences()));
            builder.worstTranscript(worstAnnotation.transcript());
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = canonicalAnnotationFactory.canonicalSnpEffAnnotation(transcriptAnnotations);
        if (canonicalAnnotation.isPresent()) {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            builder.canonicalGene(annotation.gene());
            builder.canonicalEffect(annotation.consequenceString());
            builder.canonicalCodingEffect(CodingEffect.effect(annotation.gene(), annotation.consequences()));
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

}
