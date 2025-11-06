package com.hartwig.hmftools.qsee.cohort;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import com.hartwig.hmftools.qsee.cohort.FeatureMatrix;

import org.junit.Test;

import com.hartwig.hmftools.qsee.feature.Feature;

public class FeatureMatrixTest
{
    @Test
    public void canAddRowsIncrementally()
    {
        List<Feature> sample1Features = List.of(
                createTestFeature("feature1", 1.1),
                createTestFeature("feature2", 1.2)
        );

        List<Feature> sample2Features = List.of(
                createTestFeature("feature1", 2.1),
                createTestFeature("feature2", 2.2)
        );

        List<Feature> sample3Features = List.of(
                createTestFeature("feature1", 3.1),
                createTestFeature("feature2", 3.2),
                createTestFeature("feature3", 3.3)
        );

        List<Feature> sample4Features = List.of(
                createTestFeature("feature1", 4.1),
                createTestFeature("feature2", 4.2),
                createTestFeature("feature3", 4.3)
        );

        FeatureMatrix matrix = new FeatureMatrix(new HashMap<>(), 4);
        matrix.addRow("sample1", sample1Features);
        matrix.addRow("sample2", sample2Features);
        matrix.addRow("sample3", sample3Features);
        matrix.addRow("sample4", sample4Features);

        List<String> expectedSampleIds = List.of("sample1", "sample2", "sample3", "sample4");
        List<String> actualSampleIds = matrix.getRowIds();
        assertEquals(expectedSampleIds, actualSampleIds);

        List<String> actualFeatureNames = matrix.getFeatureKeys().stream().map(x -> x.name()).toList();
        List<String> expectedFeatureNames = List.of("feature1", "feature2", "feature3");
        assertEquals(expectedFeatureNames, actualFeatureNames);

        double[][] expectedValues = {
                { 1.1, 1.2, Double.NaN },
                { 2.1, 2.2, Double.NaN },
                { 3.1, 3.2, 3.3 },
                { 4.1, 4.2, 4.3 },
        };
        double[][] actualValues = matrix.getValues();
        assertArrayEquals(expectedValues, actualValues);

        double[][] expectedValuesTransposed = {
                { 1.1, 2.1, 3.1, 4.1 },
                { 1.2, 2.2, 3.2, 4.2 },
                { Double.NaN, Double.NaN, 3.3, 4.3 },
        };
        double[][] actualValuesTransposed = matrix.getValuesTransposed();

        assertArrayEquals(expectedValuesTransposed, actualValuesTransposed);

    }

    @Test
    public void canAddRowsConcurrently() throws InterruptedException
    {
        int NUM_SAMPLE_THREADS = 1000;
        int NUM_FEATURES = 500;

        List<String> expectedSampleIds = IntStream.range(0, NUM_SAMPLE_THREADS)
                .mapToObj(x -> formTestSampleId(x))
                .toList();

        FeatureMatrix matrix = new FeatureMatrix(new ConcurrentHashMap<>(), expectedSampleIds);

        // Create workers
        List<Thread> threads = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < NUM_SAMPLE_THREADS; sampleIndex++)
        {
            Thread thread = createAddRowThread(matrix, sampleIndex, NUM_FEATURES);
            thread.start();
            threads.add(thread);
        }

        for(Thread thread : threads)
        {
            thread.join();
        }

        // Check row names order
        List<String> actualSampleIds = matrix.getRowIds();

        printDiffs(expectedSampleIds, actualSampleIds);
        assertEquals(expectedSampleIds, actualSampleIds);

        // Check column names order
        List<String> expectedFeatureNames = IntStream.range(0, NUM_FEATURES).mapToObj(x -> formTestFeatureName(x)).toList();
        List<String> actualFeatureNames = matrix.getFeatureKeys().stream().map(x -> x.name()).toList();

        printDiffs(expectedFeatureNames, actualFeatureNames);
        assertEquals(expectedFeatureNames, actualFeatureNames);

        // Check values
        double[][] expectedValues = createExpectedValues(NUM_SAMPLE_THREADS, NUM_FEATURES);
        double[][] actualValues = matrix.getValues();
        assertArrayEquals(expectedValues, actualValues);
    }

    private static Feature createTestFeature(String name, double value)
    {
        return new Feature(name, value, null, null);
    }

    private static String formTestSampleId(int sampleIndex) { return String.format("sample%d", sampleIndex); }
    private static String formTestFeatureName(int featureIndex) { return String.format("feature%d", featureIndex); }
    private static double formTestFeatureValue(int sampleIndex, int featureIndex) { return sampleIndex + 0.1*featureIndex; }

    private Thread createAddRowThread(FeatureMatrix matrix, int sampleIndex, int numFeatures)
    {
        String sampleId = formTestSampleId(sampleIndex);

        List<Feature> features = new ArrayList<>();
        for(int featureIndex = 0; featureIndex < numFeatures; ++featureIndex)
        {
            String featureName = formTestFeatureName(featureIndex);
            double featureValue = formTestFeatureValue(sampleIndex, featureIndex);
            features.add(createTestFeature(featureName, featureValue));
        }

        return new Thread(() -> matrix.addRow(sampleId, features));
    }

    private double[][] createExpectedValues(int numSamples, int numFeatures)
    {
        double[][] expectedValues = new double[numSamples][numFeatures];

        for(int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures; featureIndex++)
            {
                expectedValues[sampleIndex][featureIndex] = formTestFeatureValue(sampleIndex, featureIndex);
            }
        }

        return expectedValues;
    }

    private <T> void printDiffs(List<T> expected, List<T> actual)
    {
         for(int i = 0; i < expected.size(); ++i)
         {
             if(!expected.get(i).equals(actual.get(i)))
             {
                 System.out.printf("index(%d): actual(%s) != expected(%s)\n", i, actual.get(i), expected.get(i));
             }
         }
    }
}
