package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.SnvSigUtils.variantContext;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.util.List;
import java.util.Map;

public final class TrinucleotideCounts
{
    public static double[] extractTrinucleotideCounts(final List<SomaticVariant> variants, final Map<String,Integer> bucketNameMap)
    {
        final double[] counts = new double[bucketNameMap.size()];

        for(final SomaticVariant variant : variants)
        {
            if(variant.Type != SNP)
                continue;

            if(variant.Alt.length() != 1)
                continue;

            if(variant.TrinucleotideContext.contains("N"))
                continue;

            String bucketName = variantContext(variant.Ref, variant.Alt, variant.TrinucleotideContext);
            Integer bucketIndex = bucketNameMap.get(bucketName);

            if(bucketIndex == null)
            {
                CUP_LOGGER.error("invalid bucketName({}) from var({}>{}) context={})",
                        bucketName, variant.Ref, variant.Alt, variant.TrinucleotideContext);
                continue;
            }

            ++counts[bucketIndex];
        }

        return counts;
    }
}
