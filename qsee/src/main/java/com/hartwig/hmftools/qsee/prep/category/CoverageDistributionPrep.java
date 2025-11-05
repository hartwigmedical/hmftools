package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.NotNull;

public class CoverageDistributionPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    public CoverageDistributionPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private BamMetricCoverage loadCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, sampleId);
        return BamMetricCoverage.read(filePath);
    }

    private static List<Feature> calcPropBasesWithCoverage(List<ValueFrequency> coverageBaseCounts)
    {
        long totalCount = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();

        return coverageBaseCounts.stream().map(x -> {
            double propBases = (double) x.Count / totalCount;
            FeatureKey key = new FeatureKey(String.valueOf(x.Value), FeatureType.COVERAGE_DISTRIBUTION, SOURCE_TOOL);
            return new Feature(key, propBases);
        }).toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        BamMetricCoverage coverage = loadCoverage(sampleId, sampleType);
        List<Feature> features = calcPropBasesWithCoverage(coverage.Coverage);
        return features;
    }
}
