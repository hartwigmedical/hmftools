package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.OTHER;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.quality.MsiJitterQualCache;

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
        double modelQual;

        if(mOriginalQualModel.type() == OTHER && readContextCounter.qualCache().usesMsiIndelErrorQual())
        {
            modelQual = readContextCounter.qualCache().msiIndelErrorQual();
        }
        else
        {
            modelQual = mOriginalQualModel.calculateQual(record, readIndex);
        }

        if(modelQual < 0)
            return INVALID_BASE_QUAL;

        if(mRealignedQualModels.isEmpty())
            return modelQual;

        // take the minimum qual across the models
        for(UltimaRealignedQualModel realignedQualModel : mRealignedQualModels)
        {
            if(!realignedQualModel.qualModel().canCompute())
                continue;

            double realignedModelQual = realignedQualModel.calculateQual(record, readIndex);

            if(realignedModelQual < 0)
                return INVALID_BASE_QUAL;

            modelQual = min(realignedModelQual, modelQual);
        }

        return modelQual;
    }
}
