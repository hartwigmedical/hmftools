package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.redux.DuplicateFrequency;
import com.hartwig.hmftools.qsee.common.BinnedFrequencies;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

public class DuplicateFreqPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final String FIELD_READ_COUNT = "ReadCount";

    public DuplicateFreqPrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.DUPLICATE_FREQ; }

    private List<DuplicateFrequency> loadDuplicateFrequencies(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getReduxDir(sampleId, sampleType);
        String filePath = DuplicateFrequency.generateFilename(baseDir, sampleId);
        return DuplicateFrequency.read(filePath);
    }

    private static List<Feature> normaliseAndBinCounts(List<DuplicateFrequency> dupFreqs)
    {
        long[] readGroupSizes = new long[dupFreqs.size()];
        double[] frequencies = new double[dupFreqs.size()];

        for(int i = 0; i < dupFreqs.size(); i++)
        {
            DuplicateFrequency valueFrequency = dupFreqs.get(i);
            readGroupSizes[i] = valueFrequency.ReadCount;
            frequencies[i] = valueFrequency.Count;
        }

        BinnedFrequencies readGroupSizeFrequencies = new BinnedFrequencies(readGroupSizes, frequencies);

        List<Feature> features = readGroupSizeFrequencies.formProportionalDensityFeatures(
                FIELD_READ_COUNT, FeatureType.COVERAGE_DISTRIBUTION, SOURCE_TOOL);

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
        {
            return List.of();
        }

        List<DuplicateFrequency> dupFreqs = loadDuplicateFrequencies(sampleId, sampleType);
        List<Feature> features = normaliseAndBinCounts(dupFreqs);
        return features;
    }
}
