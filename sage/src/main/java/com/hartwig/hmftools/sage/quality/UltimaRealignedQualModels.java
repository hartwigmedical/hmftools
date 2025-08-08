package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.ULTIMA_MAX_QUAL_T0;
import static com.hartwig.hmftools.sage.SageConstants.ULTIMA_MAX_QUAL_TP;
import static com.hartwig.hmftools.sage.SageConstants.ULTIMA_TP_0_BOOST;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.MICROSAT_ADJUSTMENT;

import java.util.List;

import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class UltimaRealignedQualModels
{
    private final UltimaQualModel mOriginalQualModel;
    private final List<UltimaRealignedQualModel> mRealignedQualModels;

    public UltimaRealignedQualModels(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator,
            final List<UltimaRealignedQualModel> realignedQualModels)
    {
        mOriginalQualModel = originalQualModel(readContext, ultimaQualCalculator);
        mRealignedQualModels = realignedQualModels;
    }

    public UltimaRealignedQualModels(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator)
    {
        this(readContext, ultimaQualCalculator, null);
    }

    private static UltimaQualModel originalQualModel(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator)
    {
        byte[] coreBases = Arrays.subsetArray(readContext.ReadBases, readContext.VarIndex - 1, readContext.VarIndex + 1);
        return ultimaQualCalculator.buildContext(readContext.variant(), coreBases);
    }

    public double calculateQual(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        double ultimaQual = Math.max(ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST, ULTIMA_MAX_QUAL_T0);
        if(mOriginalQualModel != null)
        {
            double modelQual;
            if(mOriginalQualModel.type() == MICROSAT_ADJUSTMENT && readContextCounter.qualCache().usesMsiIndelErrorQual())
            {
                modelQual = readContextCounter.qualCache().msiIndelErrorQual();
            }
            else
            {
                modelQual = mOriginalQualModel.calculateQual(record, readIndex);
            }

            if(modelQual < 0)
            {
                return INVALID_BASE_QUAL;
            }

            ultimaQual = min(ultimaQual, modelQual);
        }

        if(mRealignedQualModels == null)
        {
            return ultimaQual;
        }

        for(UltimaRealignedQualModel realignedUltimaQualModel : mRealignedQualModels)
        {
            MsiJitterQualCache qualCache = realignedUltimaQualModel.qualCache(
                    readContextCounter.readContext().RefBases,
                    readContextCounter.readContext().ReadBases,
                    readContextCounter.qualityCalculator(),
                    readContextCounter.sampleId());

            double modelQual;
            if(realignedUltimaQualModel.type() == MICROSAT_ADJUSTMENT && qualCache.usesMsiIndelErrorQual())
            {
                modelQual = qualCache.msiIndelErrorQual();
            }
            else
            {
                modelQual = realignedUltimaQualModel.calculateQual(record, readIndex);
            }

            if(modelQual < 0)
                return INVALID_BASE_QUAL;

            ultimaQual = min(ultimaQual, modelQual);
        }

        return ultimaQual;
    }
}
