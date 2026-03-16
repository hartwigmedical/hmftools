package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.BinnedFrequencies;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
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

        double[] coverageBinStarts = coverageBaseFrequencies.binStarts();
        double[] propBasesPerBinStart = coverageBaseFrequencies.calcProportionalDensities();

        List<Feature> features = new ArrayList<>();
        for(int i = 0; i < propBasesPerBinStart.length; i++)
        {
            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_READ_DEPTH, String.valueOf(coverageBinStarts[i]));
            FeatureKey key = new FeatureKey(featureName, FeatureType.COVERAGE_DISTRIBUTION, SOURCE_TOOL);
            Feature feature = new Feature(key, propBasesPerBinStart[i]);
            features.add(feature);
        }

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
