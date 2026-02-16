package com.hartwig.hmftools.redux.ms_model;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.MULTI_BASE_REPEAT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.REPEAT_3_5_GROUP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.math.Quantiles;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterTableRow;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.apache.logging.log4j.Level;

public class MsModelCalculator
{
    private final MsModelParams mModelParams;
    private final Map<String,PurplePurity> mSamplePurities;
    private final Map<String,List<RepeatUnitData>> mSampleRepeatUnits;

    private final Map<String,double[]> mRepeatUnitCountCoefficients;
    private final Map<String,Double> mRepeatUnitCountErrorRates;

    public MsModelCalculator(final MsModelParams modelParams)
    {
        mModelParams = modelParams;

        mSampleRepeatUnits = Maps.newHashMap();
        mSamplePurities = Maps.newHashMap();
        mRepeatUnitCountErrorRates = Maps.newHashMap();
        mRepeatUnitCountCoefficients = Maps.newHashMap();
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

    private List<RepeatUnitData> buildRepeatUnitData(final Collection<JitterCountsTable> jitterCounts)
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


                // consolidate 3-5 length repeat units
                RepeatUnitData repeatUnitData;

                if(isMultiBaseRepeat)
                {
                    repeatUnitData = repeatUnitDataList.stream()
                            .filter(x -> x.RepeatUnit.equals(REPEAT_3_5_GROUP) && x.repeatCount() == jcRow.refNumUnits())
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

        return repeatUnitDataList;
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

                errorRates[sampleIndex] = repeatUnitData.mAdjustedErrorRate;
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

    public void runValidation(
            final Map<String,PurplePurity> samplePurities, final Map<String,Collection<JitterCountsTable>> sampleJitterCounts,
            final MsModelConfig config)
    {
        RD_LOGGER.info("evaluating samples");

        Level logLevel = mSamplePurities.size() > 100 ? Level.TRACE : Level.DEBUG;

        try
        {
            String filename = config.OutputDir + File.separator + "ms_model_evaluation";
            String calcsFilename = config.OutputDir + File.separator + "ms_model_calcs";

            if(config.OutputId != null)
            {
                filename += "." + config.OutputId;
                calcsFilename += "." + config.OutputId;
            }

            filename += TSV_EXTENSION;
            calcsFilename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_SAMPLE_ID).add("Purity").add("MsIndelsPerMb").add("PredictedValue");
            writer.write(sj.toString());
            writer.newLine();

            BufferedWriter calcsWriter = createBufferedWriter(calcsFilename, false);
            sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_SAMPLE_ID).add("RepeatUnit").add("RepeatCount").add("AdjustedErrorRate").add("PredictedValue");
            calcsWriter.write(sj.toString());
            calcsWriter.newLine();

            for(Map.Entry<String,PurplePurity> entry : samplePurities.entrySet())
            {
                String sampleId = entry.getKey();

                PurplePurity purplePurity = entry.getValue();

                Collection<JitterCountsTable> jitterCounts = sampleJitterCounts.get(sampleId);

                List<RepeatUnitData> repeatUnitDataList = buildRepeatUnitData(jitterCounts);

                List<Double> repeatUnitPredictedValues = calcMsIndelPerMbValues(repeatUnitDataList);

                for(int i = 0; i < repeatUnitDataList.size(); ++i)
                {
                    RepeatUnitData repeatUnitData = repeatUnitDataList.get(i);
                    double predictedValue = repeatUnitPredictedValues.get(i);

                    sj = new StringJoiner(TSV_DELIM);
                    sj.add(sampleId);
                    sj.add(repeatUnitData.RepeatUnit);
                    sj.add(String.valueOf(repeatUnitData.repeatCount()));
                    sj.add(format("%4.3e", repeatUnitData.adjustedErrorRate()));
                    sj.add(format("%.4f", predictedValue));
                    calcsWriter.write(sj.toString());
                    calcsWriter.newLine();
                }

                double predictedMsIndelsPerMb = calcMsIndelPerMb(repeatUnitPredictedValues);

                RD_LOGGER.log(logLevel, format("sample(%s) purity(%.3f) msIndelsPerMb(actual=%.4f predicted=%.4f)",
                        sampleId, purplePurity.Purity, purplePurity.MsIndelsPerMb, predictedMsIndelsPerMb));

                sj = new StringJoiner(TSV_DELIM);
                sj.add(sampleId);
                sj.add(format("%.4f", purplePurity.Purity));
                sj.add(format("%.4f", purplePurity.MsIndelsPerMb));
                sj.add(format("%.4f", predictedMsIndelsPerMb));
                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
            calcsWriter.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to create sample evaluation file: {}", e.toString());
            System.exit(1);
        }
    }

    private double calcMsIndelPerMb(final Collection<JitterCountsTable> jitterCounts)
    {
        List<RepeatUnitData> repeatUnitDataList = buildRepeatUnitData(jitterCounts);

        List<Double> repeatUnitPredictedValues = calcMsIndelPerMbValues(repeatUnitDataList);

        return calcMsIndelPerMb(repeatUnitPredictedValues);
    }

    private List<Double> calcMsIndelPerMbValues(final List<RepeatUnitData> repeatUnitDataList)
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

            double predictedValue = coefficients[0] * x + coefficients[1] * pow(x, 3);
            calcMsIndelsPerMb.add(predictedValue);
        }

        return calcMsIndelsPerMb;
    }

    private static final double INVALID_VALUE = -1;

    private double calcMsIndelPerMb(final List<Double> calcMsIndelsPerMb)
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

    private class RepeatUnitData
    {
        public final String RepeatUnit;
        public final JitterTableRow JitterRow;

        private int mWeightedErrorCount;
        private int mWeightedTotal;
        private double mErrorRate;
        private double mAdjustedErrorRate;

        private final String mRucKey;

        public RepeatUnitData(final String repeatUnit, final JitterTableRow jitterRow)
        {
            RepeatUnit = repeatUnit;
            JitterRow = jitterRow;
            mRucKey = repeatUnitAndCountKey(RepeatUnit, jitterRow.refNumUnits());

            mAdjustedErrorRate = 0;
        }

        public String rucKey() { return mRucKey;}
        public int repeatCount() { return JitterRow.refNumUnits(); }

        public void mergeJitterRow(final JitterTableRow jcRow)
        {
            for(Map.Entry<Integer,Integer> entry : jcRow.jitterCounts().entrySet())
            {
                JitterRow.addReads(entry.getKey(), entry.getValue());
            }
        }

        public void computeErrorRates()
        {
            mWeightedErrorCount = 0;
            mWeightedTotal = 0;
            mErrorRate = 0;

            for(Map.Entry<Integer,Integer> entry : JitterRow.jitterCounts().entrySet())
            {
                int count = abs(entry.getKey());
                int reads = entry.getValue();

                if(count == 0)
                    mWeightedTotal += reads;
                else
                    mWeightedErrorCount += count * reads;
            }

            mWeightedTotal += mWeightedErrorCount;
            mErrorRate = mWeightedTotal > 0 ? mWeightedErrorCount / (double) mWeightedTotal : 0;
        }

        public void setAdjustedErrorRate(double noiseAdjustment) { mAdjustedErrorRate = max(mErrorRate - noiseAdjustment, 0); }
        public double errorRate() { return mErrorRate; }
        public double adjustedErrorRate() { return mAdjustedErrorRate; }

        public String toString()
        {
            return format("%s: reads(%d) errorRate(%4.3e)", mRucKey, mWeightedTotal, mErrorRate);
        }
    }

    public void loadModelErrorRates(final List<MsModelErrorRates> modelErrorRates)
    {
        for(MsModelErrorRates errorRates : modelErrorRates)
        {
            String rucKey = repeatUnitAndCountKey(errorRates.RepeatUnit, errorRates.RepeatCount);
            mRepeatUnitCountErrorRates.put(rucKey, errorRates.ErrorRate);
        }
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

    public void loadModelCoeffcients(final List<MsModelCoefficients> modelCoefficients)
    {
        for(MsModelCoefficients coefficients : modelCoefficients)
        {
            String rucKey = repeatUnitAndCountKey(coefficients.RepeatUnit, coefficients.RepeatCount);
            mRepeatUnitCountCoefficients.put(rucKey, new double[] { coefficients.Coefficient1, coefficients.Coefficient2} );
        }
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

    private static String repeatUnitAndCountKey(final String repeatUnit, final int repeatCount)
    {
        // eg A/T with num-units = 4
        return format("%s_%d", repeatUnit, repeatCount);
    }

    private static boolean isMultiBaseRepeat(final String repeatUnit) { return repeatUnit.contains(MULTI_BASE_REPEAT); }
}
