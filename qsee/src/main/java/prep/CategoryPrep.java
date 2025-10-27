package prep;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.Nullable;

import common.SampleType;
import feature.Feature;

public interface CategoryPrep
{
    default @Nullable List<Feature> extractSampleData(String sampleId) throws Exception
    {
        return extractSampleData(sampleId, null);
    }

    default @Nullable List<Feature> extractSampleData(String sampleId, SampleType sampleType) throws Exception
    {
        // Sample type is used to get the input directory for tools (e.g. BamMetrics) where the directory is different for
        // tumor vs reference samples.
        return extractSampleData(sampleId);
    }
}
