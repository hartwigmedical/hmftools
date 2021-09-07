package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

public class SomaticVariantQC
{
    private final int mVariantCount;
    private final double mVariantAlleleCount;

    public SomaticVariantQC(int variantCount, double variantAlleleCount)
    {
        mVariantCount = variantCount;
        mVariantAlleleCount = variantAlleleCount;
    }

    public List<String> header()
    {
        return Lists.newArrayList("SomaticVariantsMatched", "SomaticVariantsUnmatched");
    }

    public List<String> body()
    {
        return Lists.newArrayList(String.valueOf(matchedVariants()), String.valueOf(unmatchedVariants()));
    }

    public int matchedVariants() { return (int)round(min(mVariantAlleleCount, mVariantCount)); }

    public int unmatchedVariants() { return mVariantCount - matchedVariants(); }

    public static SomaticVariantQC create(int variantCount, final List<SomaticCodingCount> codingCount)
    {
        double totalCount = codingCount.stream().mapToDouble(x -> x.total()).sum();

        SomaticVariantQC result = new SomaticVariantQC(variantCount, totalCount);

        if(result.unmatchedVariants() > 0)
        {
            LL_LOGGER.warn("  UNASSIGNED_VARIANT - {} variants found but {} assigned", variantCount, totalCount);
        }

        return result;
    }
}
