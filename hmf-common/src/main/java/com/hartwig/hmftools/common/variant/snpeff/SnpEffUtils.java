package com.hartwig.hmftools.common.variant.snpeff;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;

import java.util.List;

import com.hartwig.hmftools.common.variant.CodingEffect;
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

        String canonicalGeneName = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        boolean canonicalSpliceRegion = false;
        String otherReportableEffects = "";
        CodingEffect worstCodingEffect = UNDEFINED;
        int genesAffected = 0;

        if(worst.size() == 5)
        {
            worstCodingEffect = CodingEffect.valueOf(worst.get(3));
            genesAffected = Integer.parseInt(worst.get(4));
        }

        if(canonical.size() == 6)
        {
            canonicalGeneName = canonical.get(0);
            canonicalTranscript = canonical.get(1);
            canonicalEffect = canonical.get(2);
            canonicalCodingEffect = CodingEffect.valueOf(canonical.get(3));
            canonicalHgvsCodingImpact = canonical.get(4);
            canonicalHgvsProteinImpact = canonical.get(5);

            canonicalSpliceRegion = canonicalEffect.contains("splice");
        }

        return new VariantImpact(
                canonicalGeneName, canonicalTranscript, canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, canonicalSpliceRegion, otherReportableEffects, worstCodingEffect, genesAffected);
    }
}
