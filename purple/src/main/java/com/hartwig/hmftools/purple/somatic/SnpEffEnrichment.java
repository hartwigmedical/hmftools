package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_CANONICAL;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_WORST;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.CodingEffectFactory;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SnpEffEnrichment implements VariantContextEnrichment
{
    private final Consumer<VariantContext> mConsumer;

    private final CanonicalAnnotation mCanonicalAnnotation;
    private final CodingEffectFactory mCodingEffectFactory;

    private final Set<String> mMissingGenes;

    public SnpEffEnrichment(
            final Set<String> driverGenes, final List<HmfTranscriptRegion> transcripts, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;
        mCanonicalAnnotation = new CanonicalAnnotation(driverGenes, transcripts);
        mCodingEffectFactory = new CodingEffectFactory(transcripts);
        mMissingGenes = Sets.newHashSet();
    }

    public Set<String> missingGenes() { return mMissingGenes; }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        final VariantImpact variantImpact = formVariantImpact(context);
        VariantImpactSerialiser.writeImpactDetails(context, variantImpact, SNPEFF_CANONICAL, SNPEFF_WORST);
        mConsumer.accept(context);
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader header)
    {
        return VariantImpactSerialiser.writeHeader(header, SNPEFF_CANONICAL, SNPEFF_WORST);
    }

    private VariantImpact formVariantImpact(@NotNull final VariantContext context)
    {
        boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;

        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationParser.fromContext(context);

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());

        int genesAffected = (int)transcriptAnnotations.stream()
                .map(SnpEffAnnotation::gene)
                .filter(x -> !x.isEmpty())
                .distinct()
                .count();

        String canonicalGeneName = "";
        String canonicalGeneId = "";
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
            canonicalGeneName = annotation.gene();
            canonicalEffect = annotation.consequenceString();
            canonicalCodingEffect = codingEffect(context, phasedInframeIndel, annotation);
            canonicalHgvsCodingImpact = annotation.hgvsCoding();
            canonicalHgvsProteinImpact = annotation.hgvsProtein();
            canonicalTranscript = annotation.transcript();
        }

        return new VariantImpact(
                genesAffected, canonicalGeneId, canonicalGeneName, canonicalEffect, canonicalTranscript, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, worstGene, worstEffect, worstTranscript, worstCodingEffect);
    }

    private CodingEffect codingEffect(final VariantContext context, boolean phasedInframeIndel, final SnpEffAnnotation annotation)
    {
        /*
        if(!mMissingGenes.contains(annotation.gene()) && !mCodingEffectFactory.hasGene(annotation.gene()))
        {
            PPL_LOGGER.warn("gene({}: {}) not in global list", annotation.geneID(), annotation.gene());
            mMissingGenes.add(annotation.gene());
        }
        */

        CodingEffect effect = mCodingEffectFactory.effect(context, annotation.gene(), annotation.featureID(), annotation.consequences());
        return phasedInframeIndel && effect.equals(CodingEffect.NONSENSE_OR_FRAMESHIFT) ? CodingEffect.MISSENSE : effect;
    }

    @Override
    public void flush() {}

}
