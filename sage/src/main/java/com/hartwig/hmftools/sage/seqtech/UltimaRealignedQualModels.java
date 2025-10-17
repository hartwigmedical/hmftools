package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.minQual;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.OTHER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BASE_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_DEFAULT_ERROR_RATE;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.BQR_CACHE;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class UltimaRealignedQualModels
{
    private final UltimaQualModel mOriginalQualModel;
    private final List<UltimaRealignedQualModel> mRealignedQualModels;

    public UltimaRealignedQualModels(final UltimaQualModel qualModel, final List<UltimaRealignedQualModel> realignedQualModels)
    {
        mOriginalQualModel = qualModel;
        mRealignedQualModels = realignedQualModels;
    }

    public UltimaRealignedQualModels(final UltimaQualModel qualModel)
    {
        mOriginalQualModel = qualModel;
        mRealignedQualModels = Collections.emptyList();
    }

    public double calculateQual(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        Byte modelQual = null;

        if(mOriginalQualModel.type() == OTHER && readContextCounter.qualCache().usesMsiIndelErrorQual())
        {
            modelQual = readContextCounter.qualCache().msiIndelErrorQual();
        }
        else if(mOriginalQualModel.canCompute())
        {
            modelQual = mOriginalQualModel.calculateQual(record, readIndex);
        }

        if(modelQual < 0)
            return INVALID_BASE_QUAL;

        // take the minimum qual across the models
        for(UltimaRealignedQualModel realignedQualModel : mRealignedQualModels)
        {
            if(!realignedQualModel.qualModel().canCompute())
                continue;

            byte realignedModelQual = realignedQualModel.calculateQual(record, readIndex);

            if(realignedModelQual < 0)
                return INVALID_BASE_QUAL;

            if(modelQual == null)
                modelQual = realignedModelQual;
            else
                modelQual = minQual(realignedModelQual, modelQual);
        }

        if(modelQual == null)
            modelQual = (byte)(BQR_CACHE.maxRawQual() + DEFAULT_BASE_QUAL_FIXED_PENALTY); // pick a sensible default, i.e. the max raw qual after deducting fixed penalty

        // floor the seq tech base qual due to realignment issues from nearby MSI sites dominating indel errors at quals above 40
        return readContextCounter.isIndel() ? Math.min(modelQual, probabilityToPhredQual(MSI_JITTER_DEFAULT_ERROR_RATE)) : modelQual;
    }
}
