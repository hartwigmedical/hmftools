package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class ReadContextQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<String,Double>[] mQualMapByRefIndexPosOrient;
    private final Map<String,Double>[] mQualMapByRefIndexNegOrient;
    private final QualityCalculator mQualityCalculator;

    private final RepeatInfo mMsiIndelRepeat;
    private final double mMsiIndelErrorQual;
    private final boolean mIsMsiSampleAndVariant;

    public ReadContextQualCache(final VariantReadContext readContext, final QualityCalculator qualityCalculator, final String sampleId)
    {
        mVariantPosition = readContext.variant().Position;
        mVariantAlt = readContext.variant().alt();

        mQualityCalculator = qualityCalculator;

        mMsiIndelRepeat = qualityCalculator.msiJitterCalcs().findRepeat(readContext);

        if(mMsiIndelRepeat != null)
        {
            double errorRate = qualityCalculator.msiJitterCalcs().calcErrorRate(readContext.variant(), sampleId, mMsiIndelRepeat);

            mMsiIndelErrorQual = errorRate > 0 ? BaseQualAdjustment.probabilityToPhredQual(errorRate) : INVALID_BASE_QUAL;
            mIsMsiSampleAndVariant = usesMsiIndelErrorQual() && qualityCalculator.msiJitterCalcs().getProbableMsiStatus(sampleId);
        }
        else
        {
            mMsiIndelErrorQual = INVALID_BASE_QUAL;
            mIsMsiSampleAndVariant = false;
        }

        mQualMapByRefIndexPosOrient = new HashMap[mVariantAlt.length()];
        mQualMapByRefIndexNegOrient = new HashMap[mVariantAlt.length()];

        for(int i = 0; i < mVariantAlt.length(); ++i)
        {
            mQualMapByRefIndexPosOrient[i] = new HashMap<>();
            mQualMapByRefIndexNegOrient[i] = new HashMap<>();
        }
    }

    public double msiIndelErrorQual() { return mMsiIndelErrorQual; }
    public RepeatInfo msiIndelRepeat() { return mMsiIndelRepeat; }
    public boolean usesMsiIndelErrorQual() { return mMsiIndelErrorQual != INVALID_BASE_QUAL; }
    public boolean isMsiSampleAndVariant() { return mIsMsiSampleAndVariant; }

    public double getQual(final byte baseQual, final ConsensusType readType, final int refIndex, final boolean posOrientation)
    {
        String key = String.valueOf(baseQual) + "_" + readType.ordinal();

        Map<String,Double>[] qualMapByRefIndex = posOrientation ? mQualMapByRefIndexPosOrient : mQualMapByRefIndexNegOrient;
        Double bqrQual = qualMapByRefIndex[refIndex].get(key);

        if(bqrQual != null)
            return bqrQual;

        byte[] trinucleotideContext = mQualityCalculator.getTrinucleotideContext(mVariantPosition + refIndex);
        byte altBase = (byte)mVariantAlt.charAt(refIndex);

        if(!posOrientation)
        {
            trinucleotideContext = Nucleotides.reverseComplementBases(trinucleotideContext);
            altBase = Nucleotides.swapDnaBase(altBase);
        }

        double bqrValue = mQualityCalculator.lookupRecalibrateQuality(trinucleotideContext, altBase, baseQual, readType);

        qualMapByRefIndex[refIndex].put(key, bqrValue);

        return bqrValue;
    }
}
