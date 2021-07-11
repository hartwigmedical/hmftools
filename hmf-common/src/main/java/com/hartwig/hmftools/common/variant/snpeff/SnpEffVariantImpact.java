package com.hartwig.hmftools.common.variant.snpeff;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.CodingEffectFactory;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SnpEffVariantImpact
{
    private final CanonicalAnnotation mCanonicalAnnotation;
    private final CodingEffectFactory mCodingEffectFactory;

    public SnpEffVariantImpact(@NotNull final Set<String> driverGenes, @NotNull final List<HmfTranscriptRegion> transcripts)
    {
        mCanonicalAnnotation = new CanonicalAnnotation(driverGenes, transcripts);
        mCodingEffectFactory = new CodingEffectFactory(transcripts);
    }

    @NotNull
    public VariantImpact fromSnpEffAnnotations(@NotNull final VariantContext context)
    {
        boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        return fromSnpEffAnnotations(context, phasedInframeIndel, allAnnotations);
    }

    private VariantImpact fromSnpEffAnnotations(
            final VariantContext context, boolean phasedInframeIndel, final List<SnpEffAnnotation> allAnnotations)
    {
        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());

        int genesAffected = (int)transcriptAnnotations.stream()
                .map(SnpEffAnnotation::gene)
                .filter(x -> !x.isEmpty())
                .distinct()
                .count();

        String canonicalGene = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        String worstGene = "";
        String worstEffect = "";
        String worstTranscript = "";
        CodingEffect worstCodingEffect = UNDEFINED;

        if(!transcriptAnnotations.isEmpty())
        {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            worstGene = worstAnnotation.gene();
            worstEffect = worstAnnotation.consequenceString();
            worstCodingEffect = codingEffect(context, phasedInframeIndel, worstAnnotation);
            worstTranscript = worstAnnotation.transcript();
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = mCanonicalAnnotation.canonicalSnpEffAnnotation(transcriptAnnotations);
        if(canonicalAnnotation.isPresent())
        {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            canonicalGene = annotation.gene();
            canonicalEffect = annotation.consequenceString();
            canonicalCodingEffect = codingEffect(context, phasedInframeIndel, annotation);
            canonicalHgvsCodingImpact = annotation.hgvsCoding();
            canonicalHgvsProteinImpact = annotation.hgvsProtein();
            canonicalTranscript = annotation.transcript();
        }

        return new VariantImpact(
                genesAffected, canonicalGene, canonicalEffect, canonicalTranscript, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, worstGene, worstEffect, worstTranscript, worstCodingEffect);
    }

    private CodingEffect codingEffect(final VariantContext context, boolean phasedInframeIndel, final SnpEffAnnotation annotation)
    {
        CodingEffect effect = mCodingEffectFactory.effect(context, annotation.gene(), annotation.consequences());
        return phasedInframeIndel && effect.equals(CodingEffect.NONSENSE_OR_FRAMESHIFT) ? CodingEffect.MISSENSE : effect;
    }
}
