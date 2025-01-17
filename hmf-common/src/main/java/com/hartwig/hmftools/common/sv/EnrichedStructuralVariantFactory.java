package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_PLOIDY_INFO;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public final class EnrichedStructuralVariantFactory
{
    public List<EnrichedStructuralVariant> enrich(final List<StructuralVariant> variants)
    {
        final List<EnrichedStructuralVariant> result = Lists.newArrayList();

        for(StructuralVariant variant : variants)
        {
            final VariantContext startContext = variant.startContext();
            if(startContext != null)
            {
                final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder()
                        .from(variant)
                        .junctionCopyNumber(junctionCopyNumber(startContext))
                        .start(createBuilder(startContext, variant.start()));

                @Nullable
                final StructuralVariantLeg endLeg = variant.end();
                @Nullable
                final VariantContext endContext = variant.endContext();
                if(endLeg != null && endContext != null)
                {
                    builder.end(createBuilder(endContext, endLeg));
                }

                result.add(builder.build());
            }
        }

        return result;
    }

    private static Double junctionCopyNumber(final VariantContext startContext)
    {
        if(startContext.hasAttribute(PURPLE_JUNCTION_COPY_NUMBER))
        {
            return startContext.getAttributeAsDouble(PURPLE_JUNCTION_COPY_NUMBER, 0);
        }

        if(startContext.hasAttribute(PURPLE_PLOIDY_INFO))
        {
            return startContext.getAttributeAsDouble(PURPLE_PLOIDY_INFO, 0);
        }

        return null;
    }

    private ImmutableEnrichedStructuralVariantLeg createBuilder(final VariantContext context, final StructuralVariantLeg leg)
    {
        final List<Double> purpleAF = context.hasAttribute(PURPLE_AF) ?
                context.getAttributeAsDoubleList(PURPLE_AF, 0.0) : Collections.emptyList();

        final List<Double> purpleCN = context.hasAttribute(PURPLE_CN) ?
                context.getAttributeAsDoubleList(PURPLE_CN, 0.0) : Collections.emptyList();

        final List<Double> purpleCNChange = context.hasAttribute(PURPLE_CN_CHANGE) ?
                context.getAttributeAsDoubleList(PURPLE_CN_CHANGE, 0.0) : Collections.emptyList();

        final ImmutableEnrichedStructuralVariantLeg.Builder builder = ImmutableEnrichedStructuralVariantLeg.builder()
                .from(leg)
                .adjustedAlleleFrequency(!purpleAF.isEmpty() ? purpleAF.get(0) : null)
                .adjustedCopyNumber(!purpleCN.isEmpty() ? purpleCN.get(0) : null)
                .adjustedCopyNumberChange(!purpleCNChange.isEmpty() ? purpleCNChange.get(0) : null);

        return builder.build();
    }
}
