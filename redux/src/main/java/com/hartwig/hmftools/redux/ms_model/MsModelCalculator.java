package com.hartwig.hmftools.redux.ms_model;

import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.MULTI_BASE_REPEAT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.REPEAT_3_5_GROUP;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_UNIT_MIN_READ_COUNT;
import static com.hartwig.hmftools.redux.ms_model.RepeatUnitData.repeatUnitAndCountKey;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.math.Quantiles;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterTableRow;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class MsModelCalculator
{
    private final MsModelParams mModelParams;

    private final Map<String,double[]> mRepeatUnitCountCoefficients;
    private final Map<String,Double> mRepeatUnitCountErrorRates;

    // training data
    private final Map<String,PurplePurity> mSamplePurities;
    private final Map<String,List<RepeatUnitData>> mSampleRepeatUnits;

    public MsModelCalculator(final MsModelParams modelParams)
    {
        mModelParams = modelParams;

        mSampleRepeatUnits = Maps.newHashMap();
        mSamplePurities = Maps.newHashMap();
        mRepeatUnitCountErrorRates = Maps.newHashMap();
        mRepeatUnitCountCoefficients = Maps.newHashMap();
    }

    // constructor for external use
    public MsModelCalculator(final MsModelParams modelParams, final String coefficientsFile, final String errorRatesFile)
    {
        mModelParams = modelParams;

        mSampleRepeatUnits = Maps.newHashMap();
        mSamplePurities = Maps.newHashMap();
        mRepeatUnitCountErrorRates = Maps.newHashMap();
        mRepeatUnitCountCoefficients = Maps.newHashMap();

        loadCoefficients(coefficientsFile);
        loadErrorRates(errorRatesFile);
    }

    // compute method for external use
    public double calcPredictedMsiRate(final Collection<JitterCountsTable> jitterCountsTable)
    {
        List<RepeatUnitData> repeatUnitDataList = buildRepeatUnitData(jitterCountsTable);

        List<Double> repeatUnitPredictedValues = calcMsIndelPerMbValues(repeatUnitDataList);

        double predictedMsIndelsPerMb = calcMsIndelPerMb(repeatUnitPredictedValues);

        return predictedMsIndelsPerMb;
    }

    protected void loadCoefficients(final String coefficientsFile)
    {
        List<MsModelCoefficients> modelCoefficients = MsModelCoefficients.read(coefficientsFile);
        for(MsModelCoefficients coefficients : modelCoefficients)
        {
            String rucKey = repeatUnitAndCountKey(coefficients.RepeatUnit, coefficients.RepeatCount);
            mRepeatUnitCountCoefficients.put(rucKey, new double[] { coefficients.Coefficient1, coefficients.Coefficient2} );
        }
    }

    public void loadErrorRates(final String errorRatesFile)
    {
        List<MsModelErrorRates> modelErrorRates = MsModelErrorRates.read(errorRatesFile);

        for(MsModelErrorRates errorRates : modelErrorRates)
        {
            String rucKey = repeatUnitAndCountKey(errorRates.RepeatUnit, errorRates.RepeatCount);
            mRepeatUnitCountErrorRates.put(rucKey, errorRates.ErrorRate);
        }
    }

    public Map<String,double[]> repeatUnitCountCoefficients() { return mRepeatUnitCountCoefficients; }
    public Map<String,Double> repeatUnitCountErrorRates() { return mRepeatUnitCountErrorRates; }

    public void calculateErrorRates(final Map<String,Collection<JitterCountsTable>> sampleJitterCounts)
    {
        RD_LOGGER.info("calculating error rates");

        for(Map.Entry<String,Collection<JitterCountsTable>> entry : sampleJitterCounts.entrySet())
        {
            addSampleData(entry.getKey(), entry.getValue());
        }

        calculateErrorRates();
    }

    private void addSampleData(final String sampleId, final Collection<JitterCountsTable> jitterCounts)
    {
        List<RepeatUnitData> repeatUnitDataList = buildRepeatUnitData(jitterCounts);
        mSampleRepeatUnits.put(sampleId, repeatUnitDataList);
    }

    protected List<RepeatUnitData> buildRepeatUnitData(final Collection<JitterCountsTable> jitterCounts)
    {
        int rows = jitterCounts.stream().mapToInt(x -> x.getRows().size()).sum();
        List<RepeatUnitData> repeatUnitDataList = Lists.newArrayListWithCapacity(rows);

        for(JitterCountsTable jcTable : jitterCounts)
        {
            if(!mModelParams.RepeatTypes.contains(jcTable.RepeatUnit))
                continue;

            boolean isMultiBaseRepeat = isMultiBaseRepeat(jcTable.RepeatUnit);

            for(JitterTableRow jcRow : jcTable.getRows())
            {
                if(jcRow.refNumUnits() < mModelParams.RepeatCountMin || jcRow.refNumUnits() > mModelParams.RepeatCountMax)
                    continue;

                RepeatUnitData repeatUnitData;

                if(isMultiBaseRepeat)
                {
                    // consolidate 3-5 length repeat units
                    repeatUnitData = repeatUnitDataList.stream()
                            .filter(x -> x.RepeatUnit.equals(REPEAT_3_5_GROUP) && x.repeatCount() == jcRow.refNumUnits())
                            .filter(x -> x.JitterRow.getConsensusType() == jcRow.getConsensusType())
                            .findFirst().orElse(null);

                    if(repeatUnitData == null)
                    {
                        repeatUnitData = new RepeatUnitData(REPEAT_3_5_GROUP, jcRow);
                        repeatUnitDataList.add(repeatUnitData);
                    }
                    else
                    {
                        repeatUnitData.mergeJitterRow(jcRow);
                    }
                }
                else
                {
                    repeatUnitData = new RepeatUnitData(jcTable.RepeatUnit, jcRow);
                    repeatUnitDataList.add(repeatUnitData);
                }
            }
        }

        // take the highest count across consensus types for matched units and repeat counts
        Set<String> processedKeys = Sets.newHashSet();

        List<RepeatUnitData> filteredRepeatUnitDataList = Lists.newArrayList();

        for(RepeatUnitData repeatUnitData : repeatUnitDataList)
        {
            if(processedKeys.contains(repeatUnitData.rucKey()))
                continue;

            processedKeys.add(repeatUnitData.rucKey());

            List<RepeatUnitData> matchedData = repeatUnitDataList.stream()
                    .filter(x -> x.rucKey().equals(repeatUnitData.rucKey())).collect(Collectors.toList());

            int maxCount = matchedData.stream().mapToInt(x -> x.JitterRow.totalReadCount()).max().orElse(0);

            if(maxCount < TRAINING_REPEAT_UNIT_MIN_READ_COUNT)
                continue;

            RepeatUnitData maxByConsensusType = matchedData.stream()
                    .filter(x -> x.JitterRow.totalReadCount() == maxCount).findFirst().orElse(null);

            if(maxByConsensusType != null)
            {
                filteredRepeatUnitDataList.add(maxByConsensusType);
            }
        }

        return filteredRepeatUnitDataList;
    }

    private void calculateErrorRates()
    {
        Map<String,double[]> repeatErrorRates = Maps.newHashMap();
        int sampleCount = mSampleRepeatUnits.size();
        int sampleIndex = 0;

        for(List<RepeatUnitData> repeatUnitDataList : mSampleRepeatUnits.values())
        {
            for(RepeatUnitData repeatUnitData : repeatUnitDataList)
            {
                repeatUnitData.computeErrorRates();

                double[] errorRates = repeatErrorRates.get(repeatUnitData.rucKey());

                if(errorRates == null)
                {
                    errorRates = new double[sampleCount];
                    repeatErrorRates.put(repeatUnitData.rucKey(), errorRates);
                }

                errorRates[sampleIndex] = repeatUnitData.errorRate();
            }

            ++sampleIndex;
        }

        for(Map.Entry<String,double[]> entry : repeatErrorRates.entrySet())
        {
            String rucKey = entry.getKey();
            double[] errorRates = entry.getValue();;
            double quantileErrorRate = Quantiles.percentiles().index(mModelParams.NoiseThreshold).compute(errorRates);

            mRepeatUnitCountErrorRates.put(rucKey, quantileErrorRate);
        }

        // and factor this out of the per-sample data

        for(List<RepeatUnitData> repeatUnitDataList : mSampleRepeatUnits.values())
        {
            for(RepeatUnitData repeatUnitData : repeatUnitDataList)
            {
                double combinedErrorRate = mRepeatUnitCountErrorRates.get(repeatUnitData.rucKey());
                repeatUnitData.setAdjustedErrorRate(combinedErrorRate);
            }
        }
    }

    public void calculateCoefficients(final Map<String,PurplePurity> samplePurities)
    {
        RD_LOGGER.info("calculating model coefficients");

        mSamplePurities.putAll(samplePurities);

        runFittingRoutine();
    }

    private void runFittingRoutine()
    {
        // group sample error rates and ms-indels per MB into vectors for the linear regression
        int sampleCountBelowThreshold = (int)mSamplePurities.values().stream()
                .filter(x -> x.Purity <= mModelParams.SamplePurityThreshold).count();

        Map<String,double[]> rucErrorRates = Maps.newHashMap();

        for(String rucKey : mRepeatUnitCountErrorRates.keySet())
        {
            rucErrorRates.put(rucKey, new double[sampleCountBelowThreshold]);
        }

        double[] msIndelsPerMB = new double[sampleCountBelowThreshold];

        int sampleIndex = 0;
        for(Map.Entry<String,PurplePurity> entry : mSamplePurities.entrySet())
        {
            String sampleId = entry.getKey();
            PurplePurity purplePurity = mSamplePurities.get(sampleId);

            if(purplePurity.Purity > mModelParams.SamplePurityThreshold)
                continue;

            msIndelsPerMB[sampleIndex] = purplePurity.MsIndelsPerMb;

            List<RepeatUnitData> repeatUnitDataList = mSampleRepeatUnits.get(sampleId);

            for(RepeatUnitData repeatUnitData : repeatUnitDataList)
            {
                double[] errorRates = rucErrorRates.get(repeatUnitData.rucKey());

                errorRates[sampleIndex] = repeatUnitData.adjustedErrorRate();
            }

            ++sampleIndex;
        }

        for(Map.Entry<String,double[]> entry : rucErrorRates.entrySet())
        {
            String rucKey = entry.getKey();
            double[] errorRates = entry.getValue();

            double[] coefficients = calcCoefficients(errorRates, msIndelsPerMB);

            if(coefficients == null)
            {
                RD_LOGGER.warn("{} coefficient calc failed", rucKey);
                mRepeatUnitCountCoefficients.put(rucKey, new double[] {0, 0});
            }
            else
            {
                RD_LOGGER.debug(format("%s coefficient(b1=%.6f b2=%.6f)", rucKey, coefficients[0], coefficients[1]));
                mRepeatUnitCountCoefficients.put(rucKey, coefficients);
            }
        }
    }

    private static double[] calcCoefficients(final double[] errorRates, final double[] msIndelsPerMB)
    {
        double[][] errorRateArray = new double[errorRates.length][2];

        // model is y = β1 x + β2 x^3, where y = msIndelsPerMB and x is the adjusted error thresholds
        for (int i = 0; i < errorRates.length; i++)
        {
            errorRateArray[i][0] = errorRates[i];
            errorRateArray[i][1] = pow(errorRates[i], 3);
        }

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        regression.setNoIntercept(true);

        regression.newSampleData(msIndelsPerMB, errorRateArray);

        try
        {
            double[] beta = regression.estimateRegressionParameters();
            return beta;
        }
        catch(Exception e)
        {
            RD_LOGGER.warn("linear regression error: {}", e.toString());
            return null;
        }
    }

    protected List<Double> calcMsIndelPerMbValues(final List<RepeatUnitData> repeatUnitDataList)
    {
        List<Double> calcMsIndelsPerMb = Lists.newArrayList();

        for(RepeatUnitData repeatUnitData : repeatUnitDataList)
        {
            repeatUnitData.computeErrorRates();

            double combinedErrorRate = mRepeatUnitCountErrorRates.get(repeatUnitData.rucKey());
            double[] coefficients = mRepeatUnitCountCoefficients.get(repeatUnitData.rucKey());

            if(coefficients == null || coefficients[0] == 0 || coefficients[1] == 0)
            {
                calcMsIndelsPerMb.add(INVALID_VALUE);
                continue;
            }

            repeatUnitData.setAdjustedErrorRate(combinedErrorRate);
            double x = repeatUnitData.adjustedErrorRate();

            double predictedValue = max(coefficients[0] * x + coefficients[1] * pow(x, 3), 0);
            calcMsIndelsPerMb.add(predictedValue);
        }

        return calcMsIndelsPerMb;
    }

    private static final double INVALID_VALUE = -1;

    protected double calcMsIndelPerMb(final List<Double> calcMsIndelsPerMb)
    {
        int validCount = (int)calcMsIndelsPerMb.stream().filter(x -> x != INVALID_VALUE).count();

        if(validCount == 0)
            return INVALID_VALUE;

        double[] validValues = new double[validCount];
        int index = 0;

        for(Double value : calcMsIndelsPerMb)
        {
            if(value != INVALID_VALUE)
            {
                 validValues[index++] = value;
            }
        }

        double avergeMsIndelsPerMb = Doubles.mean(validValues);
        return avergeMsIndelsPerMb;
    }

    public List<MsModelErrorRates> getModelErrorRates()
    {
        List<MsModelErrorRates> modelErrorRates = Lists.newArrayListWithCapacity(mRepeatUnitCountErrorRates.size());

        for(Map.Entry<String,Double> entry : mRepeatUnitCountErrorRates.entrySet())
        {
            String repeatUnit = repeatUnitFromKey(entry.getKey());
            int repeatCount = repeatCountFromKey(entry.getKey());
            double errorRate = entry.getValue();
            modelErrorRates.add(new MsModelErrorRates(repeatUnit, repeatCount, errorRate));
        }

        return modelErrorRates;
    }

    public List<MsModelCoefficients> getModelCoeffcients()
    {
        List<MsModelCoefficients> modelCoefficients = Lists.newArrayListWithCapacity(mRepeatUnitCountCoefficients.size());

        for(Map.Entry<String,double[]> entry : mRepeatUnitCountCoefficients.entrySet())
        {
            String repeatUnit = repeatUnitFromKey(entry.getKey());
            int repeatCount = repeatCountFromKey(entry.getKey());
            double[] coefficients = entry.getValue();
            modelCoefficients.add(new MsModelCoefficients(repeatUnit, repeatCount, coefficients[0], coefficients[1]));
        }

        return modelCoefficients;
    }

    private static String repeatUnitFromKey(final String key)
    {
        String[] items = key.split("_", 2);
        return items[0];
    }

    private static int repeatCountFromKey(final String key)
    {
        String[] items = key.split("_", 2);
        return Integer.parseInt(items[1]);
    }

    private static boolean isMultiBaseRepeat(final String repeatUnit) { return repeatUnit.contains(MULTI_BASE_REPEAT); }
}
