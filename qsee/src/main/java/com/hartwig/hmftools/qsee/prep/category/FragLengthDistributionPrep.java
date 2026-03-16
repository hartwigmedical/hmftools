package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricFragmentLength;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.BinnedFrequencies;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

import org.jetbrains.annotations.NotNull;

public class FragLengthDistributionPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private final static SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private final static String FIELD_FRAG_LENGTH = "FragLength";

    public FragLengthDistributionPrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.FRAG_LENGTH_DISTRIBUTION; }

    private BamMetricFragmentLength loadFragmentLengths(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricFragmentLength.generateFilename(baseDir, sampleId);
        return BamMetricFragmentLength.read(filePath);
    }

    private static List<Feature> calcPropFragmentsWithLength(List<ValueFrequency> fragmentLengths)
    {
        BinnedFrequencies fragmentLengthFrequencies = BinnedFrequencies.fromValueFrequencies(fragmentLengths);

        List<Feature> features = fragmentLengthFrequencies.formProportionalDensityFeatures(
                FIELD_FRAG_LENGTH, FeatureType.FRAG_LENGTH_DISTRIBUTION, SOURCE_TOOL);

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        BamMetricFragmentLength fragmentLengths = loadFragmentLengths(sampleId, sampleType);
        List<Feature> features = calcPropFragmentsWithLength(fragmentLengths.FragmentLengths);
        return features;
    }
}
