package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class SnpEffUtils
{
    public static final String SNPEFF_WORST = "SEW";
    public static final String SNPEFF_CANONICAL = "SEC";

    // now only used in unit test
    @NotNull
    public static VariantImpact fromSnpEffEnrichedVariant(@NotNull final VariantContext context)
    {
        final List<String> worst = context.getAttributeAsStringList(SNPEFF_WORST, Strings.EMPTY);
        final List<String> canonical = context.getAttributeAsStringList(SNPEFF_CANONICAL, Strings.EMPTY);

        // return SnpEffSummarySerialiser.fromDetails(worst, canonical);
        return VariantImpactSerialiser.fromVcfAnnotation(worst, canonical);
    }
}
