package prep;

import java.util.List;

import feature.Feature;

public interface CategoryPrep
{
    List<Feature> extractSampleData(String sampleId);
}
