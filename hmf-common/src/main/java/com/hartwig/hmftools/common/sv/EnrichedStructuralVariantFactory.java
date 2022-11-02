package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_CONTEXT_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_INFO;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public final class EnrichedStructuralVariantFactory {

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull final List<StructuralVariant> variants) {
        final List<EnrichedStructuralVariant> result = Lists.newArrayList();

        for (StructuralVariant variant : variants) {
            final VariantContext startContext = variant.startContext();
            if (startContext != null) {

                final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder()
                        .from(variant)
                        .junctionCopyNumber(junctionCopyNumber(startContext))
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

    private static Double junctionCopyNumber(@NotNull final VariantContext startContext) {
        if (startContext.hasAttribute(StructuralVariantHeader.PURPLE_JUNCTION_COPY_NUMBER_INFO)) {
            return startContext.getAttributeAsDouble(StructuralVariantHeader.PURPLE_JUNCTION_COPY_NUMBER_INFO, 0);
        }

        if (startContext.hasAttribute(StructuralVariantHeader.PURPLE_PLOIDY_INFO)) {
            return startContext.getAttributeAsDouble(StructuralVariantHeader.PURPLE_PLOIDY_INFO, 0);
        }

        return null;
    }

    @NotNull
    private ImmutableEnrichedStructuralVariantLeg createBuilder(@NotNull final VariantContext context,
            @NotNull final StructuralVariantLeg leg) {
        final List<Double> purpleAF =
                context.hasAttribute(PURPLE_AF_INFO) ? context.getAttributeAsDoubleList(PURPLE_AF_INFO, 0.0) : Collections.emptyList();

        final List<Double> purpleCN =
                context.hasAttribute(PURPLE_CN_INFO) ? context.getAttributeAsDoubleList(PURPLE_CN_INFO, 0.0) : Collections.emptyList();

        final List<Double> purpleCNChange =
                context.hasAttribute(StructuralVariantHeader.PURPLE_CN_CHANGE_INFO) ? context.getAttributeAsDoubleList(StructuralVariantHeader.PURPLE_CN_CHANGE_INFO, 0.0) : Collections.emptyList();

        final ImmutableEnrichedStructuralVariantLeg.Builder builder = ImmutableEnrichedStructuralVariantLeg.builder()
                .from(leg)
                .refGenomeContext(context.getAttributeAsString(REF_CONTEXT_FLAG, null))
                .adjustedAlleleFrequency(!purpleAF.isEmpty() ? purpleAF.get(0) : null)
                .adjustedCopyNumber(!purpleCN.isEmpty() ? purpleCN.get(0) : null)
                .adjustedCopyNumberChange(!purpleCNChange.isEmpty() ? purpleCNChange.get(0) : null);
        return builder.build();
    }
}
