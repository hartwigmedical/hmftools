package com.hartwig.hmftools.common.sv;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sv.gridss.GridssSvFactory;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;

public final class StructuralVariantFileLoader
{
    public static List<StructuralVariant> fromFile(final String vcfFileLocation, final VariantContextFilter filter) throws IOException
    {
        final StructuralVariantFactory factory = StructuralVariantFactory.build(filter);

        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(), false))
        {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }

    public static List<StructuralVariant> fromGridssFile(final String vcfFileLocation, final VariantContextFilter filter) throws IOException
    {
        final GridssSvFactory factory = GridssSvFactory.build(filter);

        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(), false))
        {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }
}
