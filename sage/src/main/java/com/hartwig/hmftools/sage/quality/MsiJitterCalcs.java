package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParamsFile;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

public class MsiJitterCalcs
{
    private final Map<String,List<JitterModelParams>> mSampleParams;

    public MsiJitterCalcs()
    {
        mSampleParams = Maps.newHashMap();
    }

    public static MsiJitterCalcs build(final List<String> sampleIds, @Nullable final String jitterParamsDir)
    {
        MsiJitterCalcs msiJitterCalcs = new MsiJitterCalcs();

        if(jitterParamsDir != null)
        {
            msiJitterCalcs.loadSampleJitterParams(sampleIds, jitterParamsDir);
        }

        return msiJitterCalcs;
    }

    public boolean loadSampleJitterParams(final List<String> sampleIds, final String jitterParamsDir)
    {
        try
        {
            for(String sampleId : sampleIds)
            {
                String jitterParamFile = JitterModelParamsFile.generateFilename(jitterParamsDir, sampleId);

                if(!Files.exists(Paths.get(jitterParamFile)))
                    return false;

                List<JitterModelParams> jitterParams = JitterModelParamsFile.read(jitterParamFile);
                mSampleParams.put(sampleId, jitterParams);
            }

            SG_LOGGER.debug("loaded {} fitter param files", sampleIds.size());
        }
        catch(Exception e)
        {
            return false;
        }

        return true;
    }

    public static RepeatInfo getVariantRepeatInfo(final VariantReadContext readContext)
    {
        if(!readContext.variant().isIndel())
            return null;

        int repeatIndexStart = readContext.variantRefIndex() + 1;

        return RepeatInfo.findMaxRepeat(
                readContext.RefBases, repeatIndexStart, repeatIndexStart, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1,
                false, repeatIndexStart);
    }

    public double calcErrorRate(final VariantReadContext readContext, final String sampleId)
    {
        if(!readContext.variant().isIndel())
            return 0;

        int repeatIndexStart = readContext.variantRefIndex() + 1;

        RepeatInfo refRepeat = RepeatInfo.findMaxRepeat(
                readContext.RefBases, repeatIndexStart, repeatIndexStart, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1,
                false, repeatIndexStart);

        if(refRepeat == null)
            return 0;

        // check if the alt adjusts the repeat by +/- one unit
        String altBases = readContext.variant().isInsert() ?
                readContext.variant().Alt.substring(1) : readContext.variant().Ref.substring(1);

        if(!altBases.equals(refRepeat.Bases))
            return 0;

        List<JitterModelParams> allParams = mSampleParams.get(sampleId);

        if(allParams == null)
            return 0;

        JitterModelParams varParams = findApplicableParams(allParams, refRepeat);

        if(varParams == null)
            return 0;

        if(refRepeat.Count == 4)
            return varParams.OptimalScaleRepeat4;
        else if(refRepeat.Count == 5)
            return varParams.OptimalScaleRepeat5;
        else if(refRepeat.Count == 6)
            return varParams.OptimalScaleRepeat6;

        // apply skew model

        return 0;
    }

    private static final int MAX_SPECIFIC_REPEAT_LENGTH = 2;

    private static JitterModelParams findApplicableParams(final List<JitterModelParams> allParams, final RepeatInfo refRepeat)
    {
        for(JitterModelParams params : allParams)
        {
            if(refRepeat.repeatLength() <= MAX_SPECIFIC_REPEAT_LENGTH)
            {
                if(params.repeatUnitLength() == refRepeat.repeatLength())
                {
                    if(params.RepeatUnit.contains(refRepeat.Bases))
                        return params;
                }
            }
            else
            {
                if(params.RepeatUnit.length() > MAX_SPECIFIC_REPEAT_LENGTH)
                    return params;
            }
        }

        return null;
    }
}
