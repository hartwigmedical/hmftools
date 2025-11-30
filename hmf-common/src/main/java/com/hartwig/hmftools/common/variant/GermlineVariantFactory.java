package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;

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

        final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

        final PathogenicSummary pathogenicSummary = decorator.pathogenicSummary();

        return ImmutableGermlineVariantImpl.builder()
                .variant(VariantBuilderUtils.createVariantBuilder(tumorDepth, context, reference,null).build())
                .pathogenic(pathogenicSummary.Status.isPathogenic())
                .pathogenicity(pathogenicSummary.Status.toString())
                .build();
    }

}
