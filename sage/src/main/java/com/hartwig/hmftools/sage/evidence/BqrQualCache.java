package com.hartwig.hmftools.sage.evidence;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.sage.bqr.BqrReadType;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

public class BqrQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<String,Double>[] mQualMapByIndex;

    public BqrQualCache(final int variantPosition, final String alt)
    {
        mVariantPosition = variantPosition;
        mVariantAlt = alt;
        mQualMapByIndex = new HashMap[mVariantAlt.length()];

        for(int i = 0; i < mVariantAlt.length(); ++i)
        {
            mQualMapByIndex[i] = new HashMap<>();
        }
    }

    public double getQual(final byte baseQual, final BqrReadType readType, final int refIndex, final QualityCalculator qualityCalculator)
    {
        String key = String.valueOf(baseQual) + "_" + readType.ordinal();
        Double bqrQual = mQualMapByIndex[refIndex].get(key);

        if(bqrQual != null)
            return bqrQual;

        byte[] trinucleotideContext = qualityCalculator.getTrinucleotideContext(mVariantPosition + refIndex);

        double bqrValue = qualityCalculator.lookupRecalibrateQuality(
                trinucleotideContext, (byte)mVariantAlt.charAt(refIndex), baseQual, readType);

        mQualMapByIndex[refIndex].put(key, bqrValue);

        return bqrValue;
    }
}
