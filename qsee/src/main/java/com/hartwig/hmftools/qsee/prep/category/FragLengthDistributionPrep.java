package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricFragmentLength;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.NotNull;

public class FragLengthDistributionPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private final static SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private final static String FIELD_FRAG_LENGTH = "FragLength";

    public FragLengthDistributionPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private BamMetricFragmentLength loadFragmentLengths(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricFragmentLength.generateFilename(baseDir, sampleId);
        return BamMetricFragmentLength.read(filePath);
    }

    private static List<Feature> calcPropFragmentsWithLength(List<ValueFrequency> fragmentLengthCounts)
    {
        long totalFragments = fragmentLengthCounts.stream().mapToLong(x -> x.Count).sum();

        return fragmentLengthCounts.stream().map(x -> {
            double propBases = (double) x.Count / totalFragments;

            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_FRAG_LENGTH, String.valueOf(x.Value));
            FeatureKey key = new FeatureKey(featureName, FeatureType.FRAG_LENGTH_DISTRIBUTION, SOURCE_TOOL);

            return new Feature(key, propBases);
        }).toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        BamMetricFragmentLength fragmentLengths = loadFragmentLengths(sampleId, sampleType);
        List<Feature> features = calcPropFragmentsWithLength(fragmentLengths.FragmentLengths);
        return features;
    }
}
