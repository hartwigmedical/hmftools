package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.quality.MsiJitterCalcs.getVariantRepeatInfo;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

public class ReadContextQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<String,Double>[] mQualMapByIndex;
    private final QualityCalculator mQualityCalculator;
    private final double mMsiIndelErrorRate;

    public ReadContextQualCache(final VariantReadContext readContext, final QualityCalculator qualityCalculator, final String sampleId)
    {
        mVariantPosition = readContext.variant().Position;
        mVariantAlt = readContext.variant().alt();

        mQualityCalculator = qualityCalculator;

        RepeatInfo repeatInfo = getVariantRepeatInfo(readContext);
        mMsiIndelErrorRate = qualityCalculator.msiJitterCalcs().calcErrorRate(readContext, sampleId);

        mQualMapByIndex = new HashMap[mVariantAlt.length()];

        for(int i = 0; i < mVariantAlt.length(); ++i)
        {
            mQualMapByIndex[i] = new HashMap<>();
        }
    }

    public double msiIndelErrorRate() { return mMsiIndelErrorRate; }

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
