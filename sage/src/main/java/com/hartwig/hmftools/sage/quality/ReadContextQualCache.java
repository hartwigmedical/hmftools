package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BaseQualAdjustment;
import com.hartwig.hmftools.common.qual.BqrReadStrand;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class ReadContextQualCache
{
    private final int mVariantPosition;
    private final String mVariantAlt;
    private final Map<String,Double>[] mQualMapByIndex;
    private final QualityCalculator mQualityCalculator;
    private final Map<BqrReadType, Double> mMsiIndelErrorQual;
    private final boolean mIsMsiSampleAndVariant;

    public ReadContextQualCache(final VariantReadContext readContext, final QualityCalculator qualityCalculator, final String sampleId)
    {
        mVariantPosition = readContext.variant().Position;
        mVariantAlt = readContext.variant().alt();

        mQualityCalculator = qualityCalculator;

        mMsiIndelErrorQual = Maps.newHashMap();
        double errorRate;
        Set<BqrReadType> assignedReadTypes;
        List<MsiModelParams> sampleParams = qualityCalculator.msiJitterCalcs().getSampleParams(sampleId);
        if(sampleParams == null)
            assignedReadTypes = Arrays.stream(BqrReadType.values()).collect(Collectors.toSet());
        else
            assignedReadTypes = sampleParams.stream().map(x->x.params().ConsensusType).collect(Collectors.toSet());
        for(BqrReadType readType : assignedReadTypes)
        {
            errorRate = qualityCalculator.msiJitterCalcs().calcErrorRate(readContext, sampleId, readType);
            mMsiIndelErrorQual.put(readType, errorRate > 0 ? BaseQualAdjustment.probabilityToPhredQual(errorRate) : INVALID_BASE_QUAL);
        }
        mIsMsiSampleAndVariant = usesMsiIndelErrorQual() && qualityCalculator.msiJitterCalcs().getProbableMsiStatus(sampleId);

        mQualMapByIndex = new HashMap[mVariantAlt.length()];

        for(int i = 0; i < mVariantAlt.length(); ++i)
        {
            mQualMapByIndex[i] = new HashMap<>();
        }
    }

    public double msiIndelErrorQual(final BqrReadType readType) { return mMsiIndelErrorQual.get(readType); }
    public boolean usesMsiIndelErrorQual()
    {
        for(Map.Entry<BqrReadType, Double> entry : mMsiIndelErrorQual.entrySet())
        {
            if (entry.getValue() != INVALID_BASE_QUAL)
                return true;
        }
        return false;
    }
    public boolean isMsiSampleAndVariant() { return mIsMsiSampleAndVariant; }

    public double getQual(final byte baseQual, final BqrReadType readType, final BqrReadStrand readStrand, final int refIndex)
    {
        String key = String.valueOf(baseQual) + "_" + readType.ordinal() + "_" + readStrand.ordinal();
        Double bqrQual = mQualMapByIndex[refIndex].get(key);

        if(bqrQual != null)
            return bqrQual;

        byte[] trinucleotideContext = mQualityCalculator.getTrinucleotideContext(mVariantPosition + refIndex);

        double bqrValue = mQualityCalculator.lookupRecalibrateQuality(
                trinucleotideContext, (byte)mVariantAlt.charAt(refIndex), baseQual, readType, readStrand);

        mQualMapByIndex[refIndex].put(key, bqrValue);

        return bqrValue;
    }
}
