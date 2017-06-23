package com.hartwig.hmftools.purple.structural;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class StructuralVariantFileLoader {

    @NotNull
    public static List<StructuralVariant> fromFile(final String vcfFileLocation) throws IOException {
        final StructuralVariantFactory factory = new StructuralVariantFactory();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(),
                false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }

}
