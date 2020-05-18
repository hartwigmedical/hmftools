package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.enrich.StructuralRefContextEnrichment;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public final class EnrichedStructuralVariantFactory {

    private static final String PURPLE_AF = "PURPLE_AF";
    private static final String PURPLE_CN = "PURPLE_CN";
    private static final String PURPLE_CN_CHANGE = "PURPLE_CN_CHANGE";
    private static final String PURPLE_PLOIDY = "PURPLE_PLOIDY";

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull final List<StructuralVariant> variants) {
        final List<EnrichedStructuralVariant> result = Lists.newArrayList();

        for (StructuralVariant variant : variants) {
            final VariantContext startContext = variant.startContext();
            if (startContext != null) {

                final Double purplePloidy =
                        startContext.hasAttribute(PURPLE_PLOIDY) ? startContext.getAttributeAsDouble(PURPLE_PLOIDY, 0) : null;

                final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder()
                        .from(variant)
                        .ploidy(purplePloidy)
                        .start(createBuilder(startContext, variant.start()));

                @Nullable
                final StructuralVariantLeg endLeg = variant.end();
                @Nullable
                final VariantContext endContext = variant.endContext();
                if (endLeg != null && endContext != null) {
                    builder.end(createBuilder(endContext, endLeg));
                }

                result.add(builder.build());
            }
        }

        return result;
    }

    @NotNull
    private ImmutableEnrichedStructuralVariantLeg createBuilder(@NotNull final VariantContext context,
            @NotNull final StructuralVariantLeg leg) {
        final List<Double> purpleAF =
                context.hasAttribute(PURPLE_AF) ? context.getAttributeAsDoubleList(PURPLE_AF, 0.0) : Collections.emptyList();

        final List<Double> purpleCN =
                context.hasAttribute(PURPLE_CN) ? context.getAttributeAsDoubleList(PURPLE_CN, 0.0) : Collections.emptyList();

        final List<Double> purpleCNChange =
                context.hasAttribute(PURPLE_CN_CHANGE) ? context.getAttributeAsDoubleList(PURPLE_CN_CHANGE, 0.0) : Collections.emptyList();

        final ImmutableEnrichedStructuralVariantLeg.Builder builder = ImmutableEnrichedStructuralVariantLeg.builder()
                .from(leg)
                .refGenomeContext(context.getAttributeAsString(StructuralRefContextEnrichment.REF_CONTEXT_FLAG, null))
                .adjustedAlleleFrequency(!purpleAF.isEmpty() ? purpleAF.get(0) : null)
                .adjustedCopyNumber(!purpleCN.isEmpty() ? purpleCN.get(0) : null)
                .adjustedCopyNumberChange(!purpleCNChange.isEmpty() ? purpleCNChange.get(0) : null);
        return builder.build();
    }
}
