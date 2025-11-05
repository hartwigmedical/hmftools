package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.jetbrains.annotations.Nullable;

import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class PercentileTransformer
{
    private final double[] mPercentiles;
    private final double[] mRefValues;

    private double[] mPercentilesDeduped;
    private double[] mRefValuesDeduped;

    private boolean mFitted;
    private boolean mFittedWithAllNaN;

    private PercentileTransformer(double[] percentiles, double[] refValues, boolean fitted)
    {
        mPercentiles = percentiles;
        mRefValues = refValues;
        mFitted = fitted;
    }

    private PercentileTransformer(double[] percentiles)
    {
        this(percentiles, initialiseRefValues(percentiles.length), false);
    }

    private static double[] initialiseRefValues(int length)
    {
        double[] array = new double[length];
        Arrays.fill(array, Double.NaN);
        return array;
    }

    public static PercentileTransformer withNumPercentiles(int numPercentiles)
    {
        if(numPercentiles <= 0)
            throw new IllegalArgumentException("Number of percentiles should be greater than 0");

        double increment = 100.0 / (numPercentiles - 1);
        double[] percentiles = IntStream.range(0, numPercentiles).mapToDouble(i -> i * increment).toArray();

        return new PercentileTransformer(percentiles);
    }

    public static PercentileTransformer withIncrement(double increment)
    {
        if(increment <= 0 || increment >= 100)
            throw new IllegalArgumentException("Percentile increment should be between 0 and 100");

        int numPercentiles = (int) Math.ceil(100.0 / increment) + 1;

        return withNumPercentiles(numPercentiles);
    }

    public static PercentileTransformer fromPrefitData(double[] percentiles, double[] refValues)
    {
        if(percentiles.length != refValues.length)
            throw new IllegalArgumentException("Percentiles and ref values should be the same length");

        if(notSorted(percentiles) || notSorted(refValues))
            throw new IllegalArgumentException("Percentiles and ref values should be sorted in ascending order");

        if(Arrays.stream(percentiles).anyMatch(x -> x < 0 || x > 100))
            throw new IllegalArgumentException("Percentiles should be between 0 and 100");

        if(percentiles.length < 2 || percentiles[0] != 0 || percentiles[percentiles.length-1] != 100)
            throw new IllegalArgumentException("Data for at least percentiles 0 and 100 should be provided");

        PercentileTransformer transformer = new PercentileTransformer(percentiles, refValues, true);
        transformer.dedupRefValues();

        return transformer;
    }

    private static boolean notSorted(double[] values)
    {
        return IntStream.range(1, values.length).anyMatch(i -> values[i] < values[i - 1]);
    }

    public void fit(double[] cohortValues, @Nullable FeatureKey featureKey)
    {
        if(mFitted)
            throw new IllegalStateException("PercentileTransformer already fitted");

        if(cohortValues.length == 0)
            throw new IllegalArgumentException("Cohort values must not be empty");

        double[] cohortValuesFiltered = Arrays.stream(cohortValues)
                .filter(x -> !Double.isNaN(x) && Double.isFinite(x))
                .toArray();

        int invalidValuesCount = cohortValues.length - cohortValuesFiltered.length;

        if(invalidValuesCount == cohortValues.length)
        {
            String message = String.format("All %s cohort values are NaN/Infinite (transform will return NaN)", cohortValues.length);
            message += (featureKey != null) ? " for feature: " + featureKey : "";
            QC_LOGGER.trace(message);

            mFitted = true;
            mFittedWithAllNaN = true;

            return;
        }
        else if(invalidValuesCount > 0)
        {
            String message = String.format("Removed %s NaN/Infinite value(s) from %s cohort values", invalidValuesCount, cohortValues.length);
            message += (featureKey != null) ? " for feature: " + featureKey : "";
            QC_LOGGER.trace(message);
        }

        double[] cohortValuesSorted = Arrays.copyOf(cohortValuesFiltered, cohortValuesFiltered.length);
        Arrays.sort(cohortValuesSorted);

        for(int i = 0; i < mPercentiles.length; i++)
        {
            mRefValues[i] = calcRefValueAtPercentile(mPercentiles[i], cohortValuesSorted);
        }

        dedupRefValues();

        mFitted = true;
    }

    public void fit(double[] cohortValues) { fit(cohortValues, null); }

    private static double calcRefValueAtPercentile(double percentile, double[] cohortValues)
    {

        double position = (percentile / 100) * (cohortValues.length - 1);

        int lowerIndex = (int) Math.floor(position);
        int upperIndex = (int) Math.ceil(position);

        if (lowerIndex == upperIndex)
        {
            return cohortValues[lowerIndex];
        }
        else
        {
            double fraction = position - lowerIndex;
            return linearInterpolate(cohortValues[lowerIndex], cohortValues[upperIndex], fraction);
        }
    }

    private void dedupRefValues()
    {
        Map<Double, List<Double>> refValuesToPctMap = new LinkedHashMap<>();

        for(int i = 0; i < mRefValues.length; i++)
        {
            double refValue = mRefValues[i];
            double percentile = mPercentiles[i];

            refValuesToPctMap.computeIfAbsent(refValue, k -> new ArrayList<>());
            refValuesToPctMap.get(refValue).add(percentile);
        }

        Map<Double, Double> refValuesToMeanPctMap = new LinkedHashMap<>();
        for(Double refValue : refValuesToPctMap.keySet())
        {
            List<Double> groupPercentiles = refValuesToPctMap.get(refValue);
            double meanPercentile = (groupPercentiles.size() == 1)
                    ? groupPercentiles.get(0)
                    : refValuesToPctMap.get(refValue).stream().mapToDouble(x -> x).average().getAsDouble();

            refValuesToMeanPctMap.put(refValue, meanPercentile);
        }

        mRefValuesDeduped = refValuesToMeanPctMap.keySet().stream().mapToDouble(x -> x).toArray();
        mPercentilesDeduped = refValuesToMeanPctMap.values().stream().mapToDouble(x -> x).toArray();
    }

    public double transform(double inputValue)
    {
        double[] percentiles = mPercentilesDeduped;
        double[] refValues = mRefValuesDeduped;

        if(Double.isNaN(inputValue))
            throw new IllegalArgumentException("Input value cannot be NaN");

        if(!mFitted)
            throw new IllegalStateException("PercentileTransformer not fitted");

        if(mFittedWithAllNaN)
            return Double.NaN;

        if(inputValue < refValues[0])
            return Double.NEGATIVE_INFINITY;

        if(inputValue > refValues[refValues.length - 1])
            return Double.POSITIVE_INFINITY;

        int matchIndex = Arrays.binarySearch(refValues, inputValue);
        if(matchIndex >= 0)
        {
            return percentiles[matchIndex];
        }

        int upperIndex = -matchIndex - 1; // When no exact match is found with binary search, a negative match index is returned encoding the insertion point
        int lowerIndex =  upperIndex - 1;

        double lowerRefValue = refValues[lowerIndex];
        double upperRefValue = refValues[upperIndex];
        double lowerPercentile = percentiles[lowerIndex];
        double upperPercentile = percentiles[upperIndex];

        if(lowerRefValue == inputValue)
        {
            return percentiles[lowerIndex];
        }

        double fraction = (inputValue - lowerRefValue) / (upperRefValue - lowerRefValue);
        return linearInterpolate(lowerPercentile, upperPercentile, fraction);
    }

    private static double linearInterpolate(double lowerValue, double upperValue, double fraction)
    {
        return lowerValue + fraction * (upperValue - lowerValue);
    }

    public double[] getPercentiles() { return mPercentiles; }

    public double[] getRefValues() { return mRefValues; }

    public double[] getPercentilesDeduped() { return mPercentilesDeduped; }

    public double[] getRefValuesDeduped() { return mRefValuesDeduped; }
}
