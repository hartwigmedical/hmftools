// TODO: REVIEW
package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.basequal.jitter.JitterModelParams.MAX_SPECIFIC_LENGTH_UNIT;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HD_JITTER_PARAMS;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_JITTER_PARAMS;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_DEFAULT_ERROR_RATE;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_MAX_REPEAT_CHANGE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.basequal.jitter.JitterCountsTable;
import com.hartwig.hmftools.common.basequal.jitter.JitterCountsTableFile;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParamsFile;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

public class MsiJitterCalcs
{
    private final Map<String,List<MsiModelParams>> mSampleParams;
    private final Map<String,Boolean> mProbableMsiSample;

    public MsiJitterCalcs()
    {
        mSampleParams = Maps.newHashMap();
        mProbableMsiSample = Maps.newHashMap();
    }

    public static MsiJitterCalcs build(final List<String> sampleIds, @Nullable final String jitterParamsDir, final boolean highDepthMode)
    {
        MsiJitterCalcs msiJitterCalcs = new MsiJitterCalcs();

        List<JitterModelParams> jitterDefaults = highDepthMode ? DEFAULT_HD_JITTER_PARAMS : DEFAULT_JITTER_PARAMS;

        if(jitterParamsDir != null)
        {
            if(msiJitterCalcs.loadSampleJitterParams(sampleIds, jitterParamsDir, jitterDefaults))
                return msiJitterCalcs;
        }

        for(String sampleId : sampleIds)
        {
            msiJitterCalcs.setSampleParams(sampleId, jitterDefaults);
        }

        return msiJitterCalcs;
    }

    public List<MsiModelParams> getSampleParams(final String sampleId) { return mSampleParams.get(sampleId); }
    public boolean getProbableMsiStatus(final String sampleId) { return mProbableMsiSample.get(sampleId); }

    public void setSampleParams(final String sampleId, final List<JitterModelParams> params)
    {
        mSampleParams.put(sampleId, params.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList()));
        mProbableMsiSample.put(sampleId, false);
    }

    public boolean loadSampleJitterParams(final List<String> sampleIds, final String jitterParamsDir, final List<JitterModelParams> defaultParams)
    {
        try
        {
            for(String sampleId : sampleIds)
            {
                String jitterParamFile = JitterModelParamsFile.generateFilename(jitterParamsDir, sampleId);
                String jitterCountFile = JitterCountsTableFile.generateFilename(jitterParamsDir, sampleId);

                if(!Files.exists(Paths.get(jitterParamFile)) || !Files.exists(Paths.get(jitterCountFile)))
                    return false;

                List<JitterModelParams> rawParams = JitterModelParamsFile.read(jitterParamFile);
                List<MsiModelParams> msiParams = rawParams.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList());
                List<MsiModelParams> defaultMsiParams = defaultParams.stream().map(x -> new MsiModelParams(x)).collect(Collectors.toList());
                Collection<JitterCountsTable> jitterCounts = JitterCountsTableFile.read(jitterCountFile);

                PerSampleJitterParams sampleJitterParams = shouldRevertToDefaults(msiParams, defaultMsiParams, jitterCounts);

                mSampleParams.put(sampleId, sampleJitterParams.UseDefaults ? sampleJitterParams.ParamList : msiParams);
                mProbableMsiSample.put(sampleId, sampleJitterParams.UseDefaults);
            }

            SG_LOGGER.debug("loaded {} jitter param files", sampleIds.size());
        }
        catch(Exception e)
        {
            SG_LOGGER.error("missing jitter param file: {}", e.toString());
            return false;
        }

        return true;
    }

    private class PerSampleJitterParams
    {
        public final List<MsiModelParams> ParamList;
        public final boolean UseDefaults;

        public PerSampleJitterParams(final List<MsiModelParams> paramList, final boolean useDefaults)
        {
            ParamList = paramList;
            UseDefaults = useDefaults;
        }
    }

    private PerSampleJitterParams shouldRevertToDefaults(
            final List<MsiModelParams> msiParams, final List<MsiModelParams> defaultParams, final Collection<JitterCountsTable> jitterCounts)
    {
        double comparisonScore = 0;

        List<MsiModelParams> sampleParamList = Lists.newArrayListWithCapacity(msiParams.size());

        for(JitterCountsTable unitParams : jitterCounts)
        {
            String repeatUnit = unitParams.RepeatUnit.split("/")[0];
            MsiModelParams relevantMsiParams = findApplicableParams(msiParams, repeatUnit);
            MsiModelParams relevantDefaultParams = findApplicableParams(defaultParams, repeatUnit);

            double relevantMsiSkew = relevantMsiParams.params().MicrosatelliteSkew;

            JitterModelParams sampleJitterParams = new JitterModelParams(
                    relevantDefaultParams.params().RepeatUnit, relevantDefaultParams.params().OptimalScaleRepeat4,
                    relevantDefaultParams.params().OptimalScaleRepeat5, relevantDefaultParams.params().OptimalScaleRepeat6,
                    relevantDefaultParams.params().ScaleFitGradient, relevantDefaultParams.params().ScaleFitIntercept, relevantMsiSkew);

            MsiModelParams sampleModelParams = new MsiModelParams(sampleJitterParams);

            sampleParamList.add(sampleModelParams);

            for(JitterCountsTable.Row perRepeatData : unitParams.getRows())
            {
                int refLength = perRepeatData.refNumUnits;
                for(Map.Entry<Integer, Integer> entry : perRepeatData.jitterCounts.entrySet())
                {
                    int jitterLength = entry.getKey();

                    if(Math.abs(jitterLength) > MSI_JITTER_MAX_REPEAT_CHANGE || jitterLength == 0)
                        continue;

                    Double rawScale = getScaleParam(relevantMsiParams.params(), refLength);
                    double rawErrorRate = relevantMsiParams.calcErrorRate(refLength, jitterLength, rawScale);

                    double rawPhredScore = probabilityToPhredQual(rawErrorRate);
                    Double defaultScale = getScaleParam(sampleModelParams.params(), refLength);
                    double defaultErrorRate = sampleModelParams.calcErrorRate(refLength, jitterLength, defaultScale);
                    double defaultPhredScore = probabilityToPhredQual(defaultErrorRate);
                    comparisonScore += (rawPhredScore - defaultPhredScore) * entry.getValue();
                }
            }
        }
        // TODO: need separate jitter defaults for Ultima, rather than just ignoring the existing defaults
        return new PerSampleJitterParams(sampleParamList, false);
    }

    @VisibleForTesting
    public double calcErrorRate(final VariantReadContext readContext, final String sampleId)
    {
        return calcErrorRate(readContext.variant(), readContext.VarIndex, readContext.variantRefIndex(), readContext.RefBases, readContext.ReadBases, sampleId);
    }

    public double calcErrorRate(final SimpleVariant variant, int varIndex, int variantRefIndex, final byte[] refBases,
            final byte[] readBases, final String sampleId)
    {
        if(!variant.isIndel())
            return 0;

        int repeatIndexStart = variantRefIndex + 1;
        int readRepeatIndexStart = varIndex + 1;

        RepeatInfo refRepeat = RepeatInfo.findMaxRepeat(
                refBases, repeatIndexStart, repeatIndexStart, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1,
                false, repeatIndexStart);

        RepeatInfo inferredRefRepeat = RepeatInfo.findMaxRepeat(
                readBases, readRepeatIndexStart, readRepeatIndexStart, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1,
                false, readRepeatIndexStart);

        String altBases = variant.isInsert() ? variant.Alt.substring(1) : variant.Ref.substring(1);

        RepeatInfo repeatToUse;
        if(inferredRefRepeat == null || !altBases.startsWith(inferredRefRepeat.Bases))
        {
            repeatToUse = refRepeat;
        }
        else
        {
            int refRepeatCount = refRepeat == null ? 0 : refRepeat.Count;
            int inferredRefRepeatCount = inferredRefRepeat.Count - getImpliedAltChange(variant.isDelete(), altBases, inferredRefRepeat);
            repeatToUse = new RepeatInfo(inferredRefRepeat.Index, inferredRefRepeat.Bases, Math.max(refRepeatCount, inferredRefRepeatCount));
        }

        if(repeatToUse == null)
            return 0;

        int impliedAltChange = getImpliedAltChange(variant.isDelete(), altBases, repeatToUse);

        if(impliedAltChange > MSI_JITTER_MAX_REPEAT_CHANGE || impliedAltChange == 0)
            return 0;

        List<MsiModelParams> allParams = mSampleParams.get(sampleId);

        if(allParams == null)
            return 0;

        MsiModelParams varParams = findApplicableParams(allParams, repeatToUse.Bases);

        if(varParams == null)
            return 0;

        Double fixedScale = getScaleParam(varParams.params(), repeatToUse.Count);

        return varParams.calcErrorRate(repeatToUse.Count, impliedAltChange, fixedScale);
    }

    private static int getImpliedAltChange(boolean isDelete, final String altBases, final RepeatInfo repeat)
    {
        int impliedAltChange = altBases.length() / repeat.repeatLength();
        if(isDelete)
            impliedAltChange *= -1;

        return impliedAltChange;
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
