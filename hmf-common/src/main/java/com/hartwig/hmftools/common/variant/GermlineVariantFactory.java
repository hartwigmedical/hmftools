package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class GermlineVariantFactory
{
    public static List<GermlineVariant> fromVCFFile(final String tumor, final String vcfFile) throws Exception
    {
        List<GermlineVariant> variants = Lists.newArrayList();

        final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

        for(VariantContext variantContext : reader.iterator())
        {
            GermlineVariant variant = createVariant(tumor, "", variantContext);
            variants.add(variant);
        }

        return variants;
    }

    public static GermlineVariant createVariant(final String sample, final String reference, final VariantContext context)
    {
        final VariantContextDecorator decorator = new VariantContextDecorator(context);
        final GenotypeStatus genotypeStatus = reference != null ? decorator.genotypeStatus(reference) : null;

        final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

        final VariantImpact variantImpact = decorator.variantImpact();

        final PathogenicSummary pathogenicSummary = decorator.pathogenicSummary();

        ImmutableGermlineVariantImpl.Builder builder = ImmutableGermlineVariantImpl.builder()
                .qual(decorator.qual())
                .type(decorator.type())
                .filter(decorator.filter())
                .chromosome(decorator.chromosome())
                .position(decorator.position())
                .ref(decorator.ref())
                .alt(decorator.alt())
                .allelicDepth(tumorDepth)
                .hotspot(decorator.hotspot())
                .minorAlleleCopyNumber(decorator.minorAlleleCopyNumber())
                .adjustedCopyNumber(decorator.adjustedCopyNumber())
                .adjustedVAF(decorator.adjustedVaf())
                .variantCopyNumber(decorator.variantCopyNumber())
                .mappability(decorator.mappability())
                .tier(decorator.tier())
                .trinucleotideContext(decorator.trinucleotideContext())
                .microhomology(decorator.microhomology())
                .repeatCount(decorator.repeatCount())
                .repeatSequence(decorator.repeatSequence())
                .reported(decorator.reported())
                .biallelic(decorator.biallelic())
                .gene(variantImpact.GeneName)
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .otherReportedEffects(variantImpact.OtherReportableEffects)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .genesAffected(variantImpact.GenesAffected)
                .germlineStatus(GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, "UNKNOWN")))
                .clinvarInfo(pathogenicSummary.ClinvarInfo)
                .pathogenic(pathogenicSummary.Status.isPathogenic())
                .pathogenicity(pathogenicSummary.Status.toString());

        builder.genotypeStatus(genotypeStatus != null ? genotypeStatus : GenotypeStatus.UNKNOWN);

        return builder.build();
    }

}
