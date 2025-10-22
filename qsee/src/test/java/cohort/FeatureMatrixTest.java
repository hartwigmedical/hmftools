package cohort;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import org.junit.Test;

import feature.FeatureValue;

public class FeatureMatrixTest
{
    @Test
    public void canAddRowsIncrementally()
    {
        List<FeatureValue> sample1Features = List.of(
                new FeatureValue("feature1", 1.1),
                new FeatureValue("feature2", 1.2)
        );

        List<FeatureValue> sample2Features = List.of(
                new FeatureValue("feature1", 2.1),
                new FeatureValue("feature2", 2.2)
        );

        List<FeatureValue> sample3Features = List.of(
                new FeatureValue("feature1", 3.1),
                new FeatureValue("feature2", 3.2),
                new FeatureValue("feature3", 3.3)
        );

        List<FeatureValue> sample4Features = List.of(
                new FeatureValue("feature1", 4.1),
                new FeatureValue("feature2", 4.2),
                new FeatureValue("feature3", 4.3)
        );

        FeatureMatrix matrix = new FeatureMatrix(new HashMap<>(), 4);
        matrix.addRow("sample1", sample1Features);
        matrix.addRow("sample2", sample2Features);
        matrix.addRow("sample3", sample3Features);
        matrix.addRow("sample4", sample4Features);

        assertEquals(List.of("sample1", "sample2", "sample3", "sample4"), matrix.getRowIds());
        assertEquals(List.of("feature1", "feature2", "feature3"), matrix.getFeatureKeys());

        double[][] expectedValues = {
                { 1.1, 1.2, Double.NaN },
                { 2.1, 2.2, Double.NaN },
                { 3.1, 3.2, 3.3 },
                { 4.1, 4.2, 4.3 },
        };
        double[][] actualValues = matrix.getValues();
        assertArrayEquals(expectedValues, actualValues);
    }

    @Test
    public void canAddRowsConcurrently() throws InterruptedException
    {
        int NUM_SAMPLE_THREADS = 100;
        int NUM_FEATURES = 5;

        FeatureMatrix matrix = new FeatureMatrix(new ConcurrentHashMap<>(), NUM_SAMPLE_THREADS);

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

        // Check result
        List<String> expectedSampleIds = IntStream.range(0, NUM_SAMPLE_THREADS).mapToObj(x -> formTestSampleId(x)).toList();
        assertEquals(expectedSampleIds, matrix.getRowIds());

        List<String> expectedFeatureKeys = IntStream.range(0, NUM_FEATURES).mapToObj(x -> formTestFeatureKey(x)).toList();
        assertEquals(expectedFeatureKeys, matrix.getFeatureKeys());

        double[][] expectedValues = createExpectedValues(NUM_SAMPLE_THREADS, NUM_FEATURES);
        double[][] actualValues = matrix.getValues();
        assertArrayEquals(expectedValues, actualValues);
    }

    private static String formTestSampleId(int sampleIndex) { return String.format("sample%d", sampleIndex); }
    private static String formTestFeatureKey(int featureIndex) { return String.format("feature%d", featureIndex); }
    private static double formTestFeatureValue(int sampleIndex, int featureIndex) { return sampleIndex + 0.1*featureIndex; }

    private Thread createAddRowThread(FeatureMatrix matrix, int sampleIndex, int numFeatures)
    {
        String sampleId = formTestSampleId(sampleIndex);

        List<FeatureValue> features = new ArrayList<>();
        for(int featureIndex = 0; featureIndex < numFeatures; ++featureIndex)
        {
            String featureKey = formTestFeatureKey(featureIndex);
            double featureValue = formTestFeatureValue(sampleIndex, featureIndex);
            features.add(new FeatureValue(featureKey, featureValue));
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
}
