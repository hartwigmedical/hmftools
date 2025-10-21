package cohort;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.List;

import org.junit.Test;

import feature.FeatureValue;

public class SampleFeatureMatrixTest
{
    @Test
    public void canConstructMatrixIncrementally()
    {
        List<FeatureValue> sample1Features = List.of(
                new FeatureValue("feature1", 1.1, null),
                new FeatureValue("feature2", 1.2, null)
        );

        List<FeatureValue> sample2Features = List.of(
                new FeatureValue("feature1", 2.1, null),
                new FeatureValue("feature2", 2.2, null)
        );

        List<FeatureValue> sample3Features = List.of(
                new FeatureValue("feature1", 3.1, null),
                new FeatureValue("feature2", 3.2, null),
                new FeatureValue("feature3", 3.3, null)
        );

        List<FeatureValue> sample4Features = List.of(
                new FeatureValue("feature1", 4.1, null),
                new FeatureValue("feature2", 4.2, null),
                new FeatureValue("feature3", 4.3, null)
        );

        SampleFeatureMatrix matrix = new SampleFeatureMatrix(new HashMap<>(), 4);
        matrix.addSampleFeatures("sample1", sample1Features);
        matrix.addSampleFeatures("sample2", sample2Features);
        matrix.addSampleFeatures("sample3", sample3Features);
        matrix.addSampleFeatures("sample4", sample4Features);

        assertEquals(List.of("sample1", "sample2", "sample3", "sample4"), matrix.getSampleIds());
        assertEquals(List.of("feature1", "feature2", "feature3"), matrix.getFeatureKeys());

        double[][] expectedFeatureValues = {
                {1.1, 1.2, Double.NaN},
                {2.1, 2.2, Double.NaN},
                {3.1, 3.2, 3.3},
                {4.1, 4.2, 4.3},
        };
        double[][] actualFeatureValues = matrix.getValues(Double.NaN);
        assertEquals(expectedFeatureValues, actualFeatureValues);
    }
}
