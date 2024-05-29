package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.basequal.jitter.JitterModelParams.MAX_SPECIFIC_LENGTH_UNIT;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_MAX_REPEAT_CHANGE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParamsFile;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

public class MsiJitterCalcs
{
    private final Map<String,List<MsiModelParams>> mSampleParams;

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

                List<JitterModelParams> rawParams = JitterModelParamsFile.read(jitterParamFile);
                List<MsiModelParams> modelParams = rawParams.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList());
                mSampleParams.put(sampleId, modelParams);
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

        int impliedAltChange = altBases.length() / refRepeat.repeatLength();

        if(impliedAltChange > MSI_JITTER_MAX_REPEAT_CHANGE)
            return 0;

        List<MsiModelParams> allParams = mSampleParams.get(sampleId);

        if(allParams == null)
            return 0;

        MsiModelParams varParams = findApplicableParams(allParams, refRepeat);

        if(varParams == null)
            return 0;

        if(refRepeat.Count == 4)
            return varParams.params().OptimalScaleRepeat4;
        else if(refRepeat.Count == 5)
            return varParams.params().OptimalScaleRepeat5;
        else if(refRepeat.Count == 6)
            return varParams.params().OptimalScaleRepeat6;

        return varParams.calcSkew(refRepeat.Count, impliedAltChange);
    }

    private static MsiModelParams findApplicableParams(final List<MsiModelParams> allParams, final RepeatInfo refRepeat)
    {
        for(MsiModelParams params : allParams)
        {
            if(refRepeat.repeatLength() <= MAX_SPECIFIC_LENGTH_UNIT)
            {
                if(params.params().repeatUnitLength() == refRepeat.repeatLength())
                {
                    if(params.params().RepeatUnit.contains(refRepeat.Bases))
                        return params;
                }
            }
            else if(params.params().aboveSpecificLength())
            {
                return params;
            }
        }

        return null;
    }

    @VisibleForTesting
    public void setSampleParams(final String sampleId, final List<JitterModelParams> params)
    {
        mSampleParams.put(sampleId, params.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList()));
    }
}
