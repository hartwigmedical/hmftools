package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.CodingEffectFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SnpEffSummaryFactory
{

    private final CanonicalAnnotation canonicalAnnotationFactory;
    private final CodingEffectFactory codingEffectFactory;

    public SnpEffSummaryFactory(@NotNull final Set<String> driverGenes, @NotNull final List<HmfTranscriptRegion> transcripts)
    {
        this.canonicalAnnotationFactory = new CanonicalAnnotation(driverGenes, transcripts);
        this.codingEffectFactory = new CodingEffectFactory(transcripts);
    }

    @NotNull
    // now only used in unit test
    public static SnpEffSummary fromSnpEffEnrichment(@NotNull final VariantContext context)
    {
        final List<String> worst = context.getAttributeAsStringList(SnpEffEnrichment.SNPEFF_WORST, Strings.EMPTY);
        final List<String> canonical = context.getAttributeAsStringList(SnpEffEnrichment.SNPEFF_CANONICAL, Strings.EMPTY);
        return SnpEffSummarySerialiser.fromDetails(worst, canonical);
    }

    @NotNull
    public SnpEffSummary fromSnpEffAnnotations(@NotNull final VariantContext context)
    {
        final boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        return fromSnpEffAnnotations(context, phasedInframeIndel, allAnnotations);
    }

    @NotNull
    private SnpEffSummary fromSnpEffAnnotations(@NotNull final VariantContext context, boolean phasedInframeIndel,
            @NotNull final List<SnpEffAnnotation> allAnnotations)
    {
        final ImmutableSnpEffSummary.Builder builder = SnpEffSummarySerialiser.createBuilder();

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        if(!transcriptAnnotations.isEmpty())
        {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            builder.worstGene(worstAnnotation.gene());
            builder.worstEffect(worstAnnotation.consequenceString());
            builder.worstCodingEffect(codingEffect(context, phasedInframeIndel, worstAnnotation));
            builder.worstTranscript(worstAnnotation.transcript());
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = canonicalAnnotationFactory.canonicalSnpEffAnnotation(transcriptAnnotations);
        if(canonicalAnnotation.isPresent())
        {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            builder.canonicalGene(annotation.gene());
            builder.canonicalEffect(annotation.consequenceString());
            builder.canonicalCodingEffect(codingEffect(context, phasedInframeIndel, annotation));
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

    @NotNull
    private CodingEffect codingEffect(@NotNull final VariantContext context, boolean phasedInframeIndel,
            @NotNull final SnpEffAnnotation annotation)
    {
        CodingEffect effect = codingEffectFactory.effect(context, annotation.gene(), annotation.consequences());
        return phasedInframeIndel && effect.equals(CodingEffect.NONSENSE_OR_FRAMESHIFT) ? CodingEffect.MISSENSE : effect;
    }
}
