package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.BinnedFrequencies;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

import org.jetbrains.annotations.NotNull;

public class CoverageDistributionPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private static final String FIELD_READ_DEPTH = "ReadDepth";

    public CoverageDistributionPrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.COVERAGE_DISTRIBUTION; }

    private BamMetricCoverage loadCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, sampleId);
        return BamMetricCoverage.read(filePath);
    }

    private static List<Feature> calcPropBasesWithCoverage(List<ValueFrequency> coverageBaseCounts)
    {
        BinnedFrequencies coverageBaseFrequencies = BinnedFrequencies.fromValueFrequencies(coverageBaseCounts);

        List<Feature> features = coverageBaseFrequencies.formProportionalDensityFeatures(
                FIELD_READ_DEPTH, FeatureType.COVERAGE_DISTRIBUTION, SOURCE_TOOL);

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        BamMetricCoverage coverage = loadCoverage(sampleId, sampleType);
        List<Feature> features = calcPropBasesWithCoverage(coverage.Coverage);
        return features;
    }
}
