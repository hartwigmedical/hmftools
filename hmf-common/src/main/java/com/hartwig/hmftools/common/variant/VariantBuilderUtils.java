package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import java.util.Optional;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBuilderUtils
{
    public static ImmutableVariantImpl.Builder createVariantBuilder(
            final AllelicDepth allelicDepth, final VariantContext context,
            @Nullable final String reference,
            @Nullable String rna)
    {
        final VariantContextDecorator decorator = new VariantContextDecorator(context);
        final VariantImpact variantImpact = decorator.variantImpact();

        final AllelicDepth rnaDepth = extractRnaDepth(context, rna);

        // to protect against changes in types, which aren't expected to impact downstream processing
        GermlineStatus germlineStatus = GermlineStatus.UNKNOWN;

        try
        {
            germlineStatus = GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, germlineStatus.toString()));
        }
        catch(Exception e) {}

        final GenotypeStatus genotypeStatus = reference != null ? decorator.genotypeStatus(reference) : GenotypeStatus.UNKNOWN;

        ImmutableVariantImpl.Builder builder = ImmutableVariantImpl.builder()
                .chromosome(decorator.chromosome())
                .position(decorator.position())
                .type(decorator.type())
                .gene(variantImpact.GeneName)
                .ref(decorator.ref())
                .alt(decorator.alt())
                .allelicDepth(allelicDepth)
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .qual(decorator.qual())
                .mappability(decorator.mappability())
                .filter(decorator.filter())
                .genesAffected(variantImpact.GenesAffected)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .otherReportedEffects(variantImpact.OtherReportableEffects)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .hotspot(decorator.hotspot())
                .adjustedCopyNumber(decorator.adjustedCopyNumber())
                .adjustedVAF(decorator.adjustedVaf())
                .minorAlleleCopyNumber(decorator.minorAlleleCopyNumber())
                .variantCopyNumber(decorator.variantCopyNumber())
                .biallelic(decorator.biallelic())
                .reported(decorator.reported())
                .genotypeStatus(genotypeStatus)
                .germlineStatus(germlineStatus)
                .trinucleotideContext(decorator.trinucleotideContext())
                .microhomology(decorator.microhomology())
                .repeatSequence(decorator.repeatSequence())
                .repeatCount(decorator.repeatCount())
                .tier(decorator.tier())
                .rnaDepth(rnaDepth)
                .reportableTranscripts(decorator.reportableTranscripts())
                .clinvarInfo(context.getAttributeAsString(CLNSIG, ""));

        if(context.hasAttribute(LOCAL_PHASE_SET))
        {
            builder.localPhaseSets(context.getAttributeAsIntList(LOCAL_PHASE_SET, 0));
        }

        return builder;
    }

    @Nullable
    private static AllelicDepth extractRnaDepth(VariantContext context, @Nullable String rna)
    {
        return Optional.ofNullable(context.getGenotype(rna))
                .filter(AllelicDepth::containsAllelicDepth)
                .map(AllelicDepth::fromGenotype)
                .orElse(null);
    }
}
