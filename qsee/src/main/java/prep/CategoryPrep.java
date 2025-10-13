package prep;

import java.util.List;

import feature.FeatureValue;

public interface CategoryPrep<T>
{
    List<FeatureValue<T>> extractSampleData(String sampleId);
}
