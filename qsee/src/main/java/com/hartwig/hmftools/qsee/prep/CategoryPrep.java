package com.hartwig.hmftools.qsee.prep;

import java.util.List;

import org.jetbrains.annotations.Nullable;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.category.PrepCategory;

public interface CategoryPrep
{
    SourceTool sourceTool();

    PrepCategory category();

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
