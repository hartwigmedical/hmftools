package com.hartwig.hmftools.sage.evidence;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.sage.quality.QualityCalculator;

public class BqrQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<Byte,Double>[] mQualMapByIndex;

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

    public double getQual(final byte baseQual, final int refIndex, final QualityCalculator qualityCalculator)
    {
        Double bqrQual = mQualMapByIndex[refIndex].get(baseQual);

        if(bqrQual != null)
            return bqrQual;

        byte[] trinucleotideContext = qualityCalculator.getTrinucleotideContext(mVariantPosition + refIndex);
        double bqrValue = qualityCalculator.lookupRecalibrateQuality(trinucleotideContext, (byte)mVariantAlt.charAt(refIndex), baseQual);
        mQualMapByIndex[refIndex].put(baseQual, bqrValue);

        return bqrValue;
    }
}
