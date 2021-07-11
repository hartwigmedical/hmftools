package com.hartwig.hmftools.common.sv;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;

public final class StructuralVariantFileLoader {

    private StructuralVariantFileLoader() {
    }

    @NotNull
    public static List<StructuralVariant> fromFile(@NotNull String vcfFileLocation, @NotNull VariantContextFilter filter) throws IOException {
        final StructuralVariantFactory factory = new StructuralVariantFactory(filter);

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(), false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }
}
