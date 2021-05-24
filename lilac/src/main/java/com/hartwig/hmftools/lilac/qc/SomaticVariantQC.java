package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.abs;

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

    public final List<String> header()
    {
        return Lists.newArrayList("variantCount", "variantAlleleCount");
    }

    public final List<String> body()
    {
        return Lists.newArrayList(String.valueOf(mVariantCount), String.valueOf(mVariantAlleleCount));
    }

    public final boolean unmatchedVariants()
    {
        return abs((double)mVariantCount - mVariantAlleleCount) > 0.01;
    }

    public static SomaticVariantQC create(int variantCount, final List<SomaticCodingCount> codingCount)
    {
        double totalCount = codingCount.stream().mapToDouble(x -> x.total()).sum();

        SomaticVariantQC result = new SomaticVariantQC(variantCount, totalCount);

        if(result.unmatchedVariants())
        {
            LL_LOGGER.warn("  UNASSIGNED_VARIANT - {} variants found but {} assigned", variantCount, totalCount);
        }

        return result;
    }
}
