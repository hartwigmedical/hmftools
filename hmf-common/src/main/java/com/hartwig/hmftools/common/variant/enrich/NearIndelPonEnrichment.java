package com.hartwig.hmftools.common.variant.enrich;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class NearIndelPonEnrichment implements SomaticEnrichment {

    private static final int DISTANCE = 10;
    private static final String FILTER = "NEAR_INDEL_PON";

    private final AbstractFeatureReader<VariantContext, LineIterator> reader;

    public NearIndelPonEnrichment(@NotNull final String vcfFile) {
        reader = getFeatureReader(vcfFile, new VCFCodec(), true);
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) throws IOException {
        if (context.getType() == VariantContext.Type.INDEL && context.isNotFiltered()) {

            int variantStart = context.getStart();
            int variantEnd = context.getStart() + context.getReference().length() - 1 + DISTANCE;

            try (CloseableTribbleIterator<VariantContext> iterator = reader.query(context.getContig(),
                    Math.max(0, variantStart - 5 * DISTANCE),
                    variantEnd)) {
                for (VariantContext pon : iterator) {
                    if (overlaps(pon, context)) {
                        builder.filter(FILTER);
                        return builder;
                    }
                }
            }

        }
        return builder;
    }

    @VisibleForTesting
    static boolean overlaps(@NotNull final VariantContext pon, @NotNull final VariantContext variant) {
        int variantStart = variant.getStart();
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

        int ponStart = pon.getStart();
        int ponEnd = pon.getStart() + pon.getReference().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }
}