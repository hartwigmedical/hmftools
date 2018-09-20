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

    public static NearIndelPonEnrichment germlinePon(@NotNull final String vcfFile) {
        return new NearIndelPonEnrichment(5, "GERMLINE_PON_COUNT", vcfFile);
    }

    public static NearIndelPonEnrichment somaticPon(@NotNull final String vcfFile) {
        return new NearIndelPonEnrichment(3, "SOMATIC_PON_COUNT", vcfFile);
    }

    private static final int DISTANCE = 10;
    private static final String FILTER = "NEAR_INDEL_PON";

    private final String flag;
    private final int minCount;
    private final AbstractFeatureReader<VariantContext, LineIterator> reader;

    private NearIndelPonEnrichment(int minCount, @NotNull final String ponFlag, @NotNull final String vcfFile) {
        this.reader = getFeatureReader(vcfFile, new VCFCodec(), true);
        this.flag = ponFlag;
        this.minCount = minCount;
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext variant) throws IOException {
        if (variant.isIndel() && variant.isNotFiltered()) {

            int variantStart = variant.getStart();
            int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

            try (CloseableTribbleIterator<VariantContext> iterator = reader.query(variant.getContig(),
                    Math.max(0, variantStart - 5 * DISTANCE),
                    variantEnd)) {
                for (VariantContext pon : iterator) {
                    if (isValidPonEntry(minCount, flag, pon) && overlapsPon(pon, variant)) {
                        builder.filter(FILTER);
                        return builder;
                    }
                }
            }

        }
        return builder;
    }

    @VisibleForTesting
    static boolean isValidPonEntry(int minCount, @NotNull final String ponFlag, @NotNull final VariantContext pon) {
        return pon.isIndel() && pon.getAttributeAsInt(ponFlag, 0) > minCount;
    }

    @VisibleForTesting
    static boolean overlapsPon(@NotNull final VariantContext pon, @NotNull final VariantContext variant) {
        if (!pon.isIndel() || !variant.isIndel()) {
            return false;
        }

        int variantStart = variant.getStart();
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

        int ponStart = pon.getStart();
        int ponEnd = pon.getStart() + pon.getReference().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }
}