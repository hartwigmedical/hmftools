// TODO: REVIEW
package com.hartwig.hmftools.sage.quality;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class ReadContextQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<String,Double>[] mQualMapByIndex;
    private final QualityCalculator mQualityCalculator;
    private final MsiJitterQualCache mMsiJitterQualCache;

    public ReadContextQualCache(final VariantReadContext readContext, final QualityCalculator qualityCalculator, final String sampleId)
    {
        mVariantPosition = readContext.variant().Position;
        mVariantAlt = readContext.variant().alt();

        mQualityCalculator = qualityCalculator;
        mMsiJitterQualCache = new MsiJitterQualCache(readContext, qualityCalculator, sampleId);
        mQualMapByIndex = new HashMap[mVariantAlt.length()];
        for(int i = 0; i < mVariantAlt.length(); ++i)
        {
            mQualMapByIndex[i] = new HashMap<>();
        }
    }

    public double msiIndelErrorQual() { return mMsiJitterQualCache.msiIndelErrorQual(); }
    public boolean usesMsiIndelErrorQual() { return mMsiJitterQualCache.usesMsiIndelErrorQual(); }
    public boolean isMsiSampleAndVariant() { return mMsiJitterQualCache.isMsiSampleAndVariant(); }

    public double getQual(final byte baseQual, final BqrReadType readType, final int refIndex)
    {
        String key = String.valueOf(baseQual) + "_" + readType.ordinal();
        Double bqrQual = mQualMapByIndex[refIndex].get(key);

        if(bqrQual != null)
            return bqrQual;

        byte[] trinucleotideContext = mQualityCalculator.getTrinucleotideContext(mVariantPosition + refIndex);

        double bqrValue = mQualityCalculator.lookupRecalibrateQuality(
                trinucleotideContext, (byte)mVariantAlt.charAt(refIndex), baseQual, readType);

        mQualMapByIndex[refIndex].put(key, bqrValue);

        return bqrValue;
    }
}
