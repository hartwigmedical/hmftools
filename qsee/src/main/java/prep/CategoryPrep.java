package prep;

import java.util.List;

import common.SampleType;
import feature.Feature;

public interface CategoryPrep
{
    default List<Feature> extractSampleData(String sampleId)
    {
        return extractSampleData(sampleId, null);
    }

    default List<Feature> extractSampleData(String sampleId, SampleType sampleType)
    {
        // Sample type is used to get the input directory for tools (e.g. BamMetrics) where the directory is different for
        // tumor vs reference samples.
        return extractSampleData(sampleId);
    }
}
