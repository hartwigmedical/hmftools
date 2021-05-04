package com.hartwig.hmftools.lilac.qc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

public class SomaticVariantQC
{
    private final int mVariantCount;
    private final double mVariantAlleleCount;

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
        int n = mVariantCount - mVariantCount;
        return (double) Math.abs(n) > 0.01;
    }

    public final int getVariantCount()
    {
        return mVariantCount;
    }

    public final double getVariantAlleleCount()
    {
        return mVariantAlleleCount;
    }

    public SomaticVariantQC(int variantCount, double variantAlleleCount)
    {
        mVariantCount = variantCount;
        mVariantAlleleCount = variantAlleleCount;
    }

    public static SomaticVariantQC create(int variantCount, final List<SomaticCodingCount> codingCount)
    {
        return null;

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull(codingCount, (String) "codingCount");
        Iterable iterable = $receiver$iv = (Iterable) codingCount;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            com.hartwig.hmftools.lilac.variant.SomaticCodingCount somaticCodingCount = (SomaticCodingCount) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            Double d = it.getTotal();
            collection.add(d);
        }
        double totalCount = CollectionsKt.sumOfDouble((Iterable) ((List) destination$iv$iv));
        SomaticVariantQC result = new SomaticVariantQC(variantCount, totalCount);
        if(result.unmatchedVariants())
        {
            getLogger().warn("    UNASSIGNED_VARIANT - " + variantCount + " variants found but " + totalCount + " assigned");
        }
        return result;

         */
    }
}
