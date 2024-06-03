package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.basequal.jitter.JitterModelParams.MAX_SPECIFIC_LENGTH_UNIT;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_JITTER_PARAMS;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_DEFAULT_ERROR_RATE;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_MAX_REPEAT_CHANGE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
        else
        {
            for(String sampleId : sampleIds)
            {
                msiJitterCalcs.setSampleParams(sampleId, DEFAULT_JITTER_PARAMS);
            }
        }

        return msiJitterCalcs;
    }

    public List<MsiModelParams> getSampleParams(final String sampleId) { return mSampleParams.get(sampleId); }

    public void setSampleParams(final String sampleId, final List<JitterModelParams> params)
    {
        mSampleParams.put(sampleId, params.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList()));
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

        if(readContext.variant().isDelete())
            impliedAltChange *= -1;

        if(impliedAltChange > MSI_JITTER_MAX_REPEAT_CHANGE)
            return 0;

        List<MsiModelParams> allParams = mSampleParams.get(sampleId);

        if(allParams == null)
            return 0;

        MsiModelParams varParams = findApplicableParams(allParams, refRepeat.Bases);

        if(varParams == null)
            return 0;

        Double fixedScale = getScaleParam(varParams.params(), refRepeat.Count);

        return varParams.calcErrorRate(refRepeat.Count, impliedAltChange, fixedScale);
    }

    private Double getScaleParam(final JitterModelParams params, int repeatCount)
    {
        if(repeatCount == 4)
            return params.OptimalScaleRepeat4;
        else if(repeatCount == 5)
            return params.OptimalScaleRepeat5;
        else if(repeatCount == 6)
            return params.OptimalScaleRepeat6;
        else
            return null;
    }

    public double getErrorRate(final List<MsiModelParams> allParams, final String repeatBases, int repeatCount, int repeatChange)
    {
        // if variant measures shortened count for say a repeat of 5, then implies ref was 4 so get the error rate for 4 going to 5
        double errorRate = MSI_JITTER_DEFAULT_ERROR_RATE;

        if(repeatCount < MIN_REPEAT_COUNT)
            return errorRate;

        MsiModelParams modelParams = findApplicableParams(allParams, repeatBases);

        if(modelParams == null)
            return errorRate;

        Double fixedScale = getScaleParam(modelParams.params(), repeatCount);

        return modelParams.calcErrorRate(repeatCount, repeatChange, fixedScale);
    }

    private static MsiModelParams findApplicableParams(final List<MsiModelParams> allParams, final String repeatBases)
    {
        int repeatLength = repeatBases.length();

        for(MsiModelParams params : allParams)
        {
            if(repeatLength <= MAX_SPECIFIC_LENGTH_UNIT)
            {
                if(params.params().repeatUnitLength() == repeatLength)
                {
                    if(params.params().RepeatUnit.contains(repeatBases))
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
}
