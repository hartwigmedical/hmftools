package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.MICROSATELLITE;

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

    public UltimaRealignedQualModels(
            final VariantReadContext readContext, final UltimaQualModelBuilder ultimaQualModelBuilder,
            final List<UltimaRealignedQualModel> realignedQualModels)
    {
        mOriginalQualModel = originalQualModel(readContext, ultimaQualModelBuilder);
        mRealignedQualModels = realignedQualModels;
    }

    public UltimaRealignedQualModels(final UltimaQualModel qualModel)
    {
        mOriginalQualModel = qualModel;
        mRealignedQualModels = Collections.emptyList();
    }

    public UltimaRealignedQualModels(final VariantReadContext readContext, final UltimaQualModelBuilder ultimaQualModelBuilder)
    {
        this(readContext, ultimaQualModelBuilder, null);
    }

    private static UltimaQualModel originalQualModel(final VariantReadContext readContext, final UltimaQualModelBuilder ultimaQualModelBuilder)
    {
        byte[] triNucBases = Arrays.subsetArray(readContext.ReadBases, readContext.VarIndex - 1, readContext.VarIndex + 1);
        return ultimaQualModelBuilder.buildContext(readContext.variant(), triNucBases);
    }

    public double calculateQual(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        double modelQual;

        if(mOriginalQualModel.type() == MICROSATELLITE && readContextCounter.qualCache().usesMsiIndelErrorQual())
        {
            modelQual = readContextCounter.qualCache().msiIndelErrorQual();
        }
        else
        {
            modelQual = mOriginalQualModel.calculateQual(record, readIndex);
        }

        if(modelQual < 0)
            return INVALID_BASE_QUAL;

        if(mRealignedQualModels == null)
            return modelQual;

        // take the minimum qual across the models
        for(UltimaRealignedQualModel realignedUltimaQualModel : mRealignedQualModels)
        {
            // CHECK: do all QMs have the same MSI jitter cache
            MsiJitterQualCache qualCache = realignedUltimaQualModel.qualCache(
                    readContextCounter.readContext().RefBases,
                    readContextCounter.readContext().ReadBases,
                    readContextCounter.qualityCalculator(),
                    readContextCounter.sampleId());

            double realignedModelQual;

            if(realignedUltimaQualModel.type() == MICROSATELLITE && qualCache.usesMsiIndelErrorQual())
            {
                realignedModelQual = qualCache.msiIndelErrorQual();
            }
            else
            {
                realignedModelQual = realignedUltimaQualModel.calculateQual(record, readIndex);
            }

            if(realignedModelQual < 0)
                return INVALID_BASE_QUAL;

            modelQual = min(realignedModelQual, modelQual);
        }

        return modelQual;
    }
}
