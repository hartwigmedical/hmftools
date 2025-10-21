package cohort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class PercentileTransformer
{
    private final double[] mPercentiles;
    private final double[] mRefValues;

    private double[] mPercentilesDeduped;
    private double[] mRefValuesDeduped;

    private boolean mFitted = false;

    public PercentileTransformer(double percentileInterval)
    {
        int numPercentiles = (int) Math.ceil(100.0 / percentileInterval) + 1;
        double[] percentiles = new double[numPercentiles];

        double currentPercentile = 0;
        for(int i = 0; i < numPercentiles; i++)
        {
            percentiles[i] = currentPercentile;
            currentPercentile += percentileInterval;
        }

        mPercentiles = percentiles;
        mRefValues = new double[percentiles.length];
    }

    public static PercentileTransformer fromPrefitData(double[] percentiles, double[] refValues)
    {
        return new PercentileTransformer(percentiles, refValues);
    }

    private PercentileTransformer(double[] percentiles, double[] refValues)
    {
        if(percentiles.length != refValues.length)
            throw new IllegalArgumentException("Percentiles and ref values should be the same length");

        checkSorted(percentiles, "Percentiles");
        checkSorted(refValues, "Ref values");

        if(percentiles.length < 2 || percentiles[0] != 0 || percentiles[percentiles.length - 1] != 100)
            throw new IllegalArgumentException("Data for at least percentiles 0 and 100 should be provided");

        mPercentiles = percentiles;
        mRefValues = refValues;

        dedupRefValues();

        mFitted = true;
    }

    private static void checkSorted(double[] values, String errorMessagePrefix)
    {
        for(int i = 1; i < values.length; i++)
        {
            if(values[i] < values[i - 1])
                throw new IllegalArgumentException(errorMessagePrefix + " should be sorted in ascending order");
        }
    }

    public void fit(double[] cohortValues)
    {
        if(mFitted)
            throw new IllegalStateException("PercentileTransformer already fitted");

        double[] cohortValuesSorted = Arrays.copyOf(cohortValues, cohortValues.length);
        Arrays.sort(cohortValuesSorted);

        for(int i = 0; i < mPercentiles.length; i++)
        {
            mRefValues[i] = calcRefValueAtPercentile(mPercentiles[i], cohortValuesSorted);
        }

        dedupRefValues();

        mFitted = true;
    }

    private static double calcRefValueAtPercentile(double percentile, double[] cohortValues) {

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

    private void checkFitted()
    {
        if(!mFitted)
            throw new IllegalStateException("PercentileTransformer not fitted");
    }

    public double transform(double inputValue)
    {
        checkFitted();
        return transformValueToPercentile(inputValue, mRefValuesDeduped, mPercentilesDeduped);
    }

    public double[] transform(double[] inputValues)
    {
        checkFitted();

        double[] percentiles = new double[inputValues.length];
        for(int i = 0; i < inputValues.length; i++)
        {
            percentiles[i] = transformValueToPercentile(inputValues[i], mRefValuesDeduped, mPercentilesDeduped);
        }

        return percentiles;
    }

    private double transformValueToPercentile(double inputValue, double[] refValues, double[] percentiles)
    {
        if(inputValue < refValues[0])
            return Double.NEGATIVE_INFINITY;

        if(inputValue > refValues[refValues.length - 1])
            return Double.POSITIVE_INFINITY;

        int matchIndex = java.util.Arrays.binarySearch(refValues, inputValue);
        if(matchIndex >= 0)
            return percentiles[matchIndex];

        int upperIndex = -matchIndex - 1; // When no exact match is found with binary search, a negative match index is returned encoding the insertion point
        int lowerIndex =  upperIndex - 1;

        double lowerRefValue = refValues[lowerIndex];
        double upperRefValue = refValues[upperIndex];
        double lowerPercentile = percentiles[lowerIndex];
        double upperPercentile = percentiles[upperIndex];

        if(lowerRefValue == inputValue)
            return percentiles[lowerIndex];

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
