package prep;

import java.util.List;

import feature.FeatureValue;

public interface CategoryPrep
{
    List<FeatureValue> extractSampleData(String sampleId);
}
