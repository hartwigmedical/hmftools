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
        List<FeatureValue<Double>> sample1Features = List.of(
                new FeatureValue<>("feature1", 1.1, null),
                new FeatureValue<>("feature2", 1.2, null)
        );

        List<FeatureValue<Double>> sample2Features = List.of(
                new FeatureValue<>("feature1", 2.1, null),
                new FeatureValue<>("feature2", 2.2, null)
        );

        List<FeatureValue<Double>> sample3Features = List.of(
                new FeatureValue<>("feature1", 3.1, null),
                new FeatureValue<>("feature2", 3.2, null),
                new FeatureValue<>("feature3", 3.3, null)
        );

        SampleFeatureMatrix matrix = new SampleFeatureMatrix(new HashMap<>(), List.of("sample1", "sample2", "sample3"));
        matrix.addSampleFeatures("sample1", sample1Features);
        matrix.addSampleFeatures("sample2", sample2Features);
        matrix.addSampleFeatures("sample3", sample3Features);

        assertEquals(List.of("sample1", "sample2", "sample3"), matrix.getSampleIds());
        assertEquals(List.of("feature1", "feature2", "feature3"), matrix.getFeatureKeys());

        Double[][] expectedFeatureValues = {
                {1.1, 1.2, null},
                {2.1, 2.2, null},
                {3.1, 3.2, 3.3},
        };
        Double[][] actualFeatureValues = matrix.getValues(true);
        assertEquals(expectedFeatureValues, actualFeatureValues);
    }
}
