package com.hartwig.hmftools.common.variant.structural;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class StructuralVariantFileLoader {

    @NotNull
    public static List<StructuralVariant> fromFile(@NotNull final String vcfFileLocation, final boolean filterOnPasses) throws IOException {
        final StructuralVariantFactory factory = new StructuralVariantFactory(filterOnPasses);

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(), false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }
}
